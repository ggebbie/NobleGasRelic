module NobleGasRelic

using PyCall, PyPlot, DrWatson, TMI, TMItransient, Interpolations, LinearAlgebra

export vintages_planview, vintages_section, agedistribution,
    taudeltaresponse, compare_deltaresponses, priorcovariance,
    gaussmarkovsolution, anomalymatrix, magnitude, trendmatrix,propagate

const mpl = PyNULL()
const plt = PyNULL()
const cmocean = PyNULL()
const cartopy = PyNULL()

#Initialize all Python packages - install with conda through Julia
function __init__()

    # following ClimatePlots.jl
    copy!(mpl, pyimport_conda("matplotlib", "matplotlib", "conda-forge"))
    copy!(cartopy, pyimport_conda("cartopy", "cartopy", "conda-forge"))

    #copy!(plt, pyimport_conda("matplotlib.pyplot", "matplotlib", "conda-forge"))
    #copy!(cmocean, pyimport_conda("cmocean", "cmocean", "conda-forge"))

    println("Python libraries installed")
 end

function planviewplotcartopy(c::Field{T}, depth, lims;titlelabel="section plot") where T <: Real

    cmap_seismic = get_cmap("seismic")
    #cmap_hot = get_cmap("hot_r")
    cmap_hot = get_cmap("inferno_r")
    cplan = planview(c,depth)

    fig = figure(202)
    clf()
    cenlon = -160.0
    proj0 = cartopy.crs.PlateCarree()
    proj = cartopy.crs.PlateCarree(central_longitude=cenlon)
    ax = fig.add_subplot(projection = proj)
    ax.set_global()
    #ax.coastlines()
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor="black", facecolor="black")
    
    outdir = plotsdir()
    !isdir(outdir) && mkpath(outdir) 
    outfname = plotsdir("vintage.png")
    xlbl = "longitude "*L"[\degree E]"
    ylbl = "latitude "*L"[\degree N]"
    ax.set_title(titlelabel)
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false, crs=proj0)
    gl.top_labels = false
    gl.right_labels = false

    test = ax.contourf(c.γ.lon,c.γ.lat, cplan', lims, cmap=cmap_hot, transform = proj0)

    colorbar(test,label="[%]",orientation="vertical",ticks=lims, fraction = 0.03)
    CS = ax.contour(c.γ.lon,c.γ.lat, cplan', lims, colors="k", transform = proj0)
    ax.clabel(CS, CS.levels, inline=true, fontsize=10)

    savefig(outfname)

end

function vintages_planview(params)

    @unpack vintage, depth, tinterval, longname = params 

    # doing this every time, not so efficient
    
    Δ,τ = read_stepresponse()
    g = vintagedistribution(tinterval[vintage][1],tinterval[vintage][2],Δ,τ,interp="linear")

    # get the meta-data correct on the output.
    gvintage = Field(g.tracer,g.γ,vintage,longname[vintage],"seawater mass fraction []")
    # save g to file if it hasn't been done before.
    if isapprox(depth,2000) # kludge to not write twice
        println("writing",depth,vintage)
        !isdir(DrWatson.datadir()) && mkdir(DrWatson.datadir())
        writefield(DrWatson.datadir("vintages_TMI_4x4_2012.nc"),gvintage)
    end
    
    froot = plotsdir(savename("TMI_4x4_2012",params,"png",accesses=["vintage","depth"]))
    println(froot)

    tlabel = "Vintage: "* longname[vintage] * ", depth="*string(depth)*"m"
    println(tlabel)
    #fname = "vintage_"*string(v)*"_"*string(depth)*"m.png"

    lims = vcat(collect(0:5:50),100)
    # Plan view plots
    planviewplotcartopy(100g, depth, lims, titlelabel=tlabel)
    mv(plotsdir("vintage.png"),plotsdir(froot),force=true)

end

function vintages_section(params)

    @unpack vintage, lon, tinterval, longname = params 

    lims = vcat(collect(0:5:50),100)
    # doing this every time, not so efficient
    Δ,τ = read_stepresponse()
    local g = vintagedistribution(tinterval[vintage][1],tinterval[vintage][2],Δ,τ)
    
    froot = plotsdir(savename("TMI_4x4_2012",params,"png",accesses=["vintage","lon"]))
    println(froot)

    tlabel = "Vintage: "* longname[vintage] * ", lon="*string(lon)*"E"
    println(tlabel)

    sectionplot(100g, lon, lims; titlelabel=tlabel) 

    savefig(plotsdir(froot))

end

"""
    function compare_deltaresponses(loc)

    Input: two locations

    Side-effect: plot of two delta responses and their difference
"""
function compare_deltaresponses(loc)
    # get age distribution
    g = agedistribution.(loc)
    tg = taudeltaresponse()

    Δg = g[1] - g[2]

    # PyPlot version, not currently showing
    figure(2)
    clf()
    line1, = PyPlot.plot(tg,g[1],"black",label="35°N, 152°W, 3.5 km")
    line2, = PyPlot.plot(tg,g[2],"red",label="20°S, 152°W, 3.5 km")
    line3, = PyPlot.plot(tg,Δg,"green",label="Δ")
    grid("true")
    xlabel("Lag, τ [yr]")
    ylabel("mass fraction per yr [1/yr]")
    legend()
    PyPlot.savefig(plotsdir("deltaresponse_NPACvSPAC.png"))

    # try to use Plots 
    # Plots.plot(tg,g[1],color=:black,label="35°N, 152°W, 3.5 km")
    # Plots.plot!(tg,g[2],color=:red,label="20°S, 152°W, 3.5 km")
    # Plots.plot!(tg,g[1]-g[2],color=:green,label="Δ")
    # plot!(xlabel="Lag, τ [yr]",ylabel="mass fraction per yr [1/yr]")
    # Plots.savefig(plotsdir("deltaresponse_NPACvSPAC.png"))

end

"""
    function priorcovariance(tₚ,σₓ,σlong,Tlong)

    tₚ: estimate times 
    σₓ: year-to-year variations
    σlong: variance on long timescale, T
    T: long timescale, timescale of correlation
"""
function priorcovariance(tₚ,σₓ,σlong,Tlong)

    # permit centennial-scale correlation in SLP
    nt = length(tₚ)
    Clong = Matrix{Float64}(undef,nt,nt)

    for xx = 1:nt
        for yy = 1:nt
            Clong[xx,yy] = σlong^2 * exp(-(tₚ[xx]-tₚ[yy])^2/Tlong^2)
        end
    end

    Cₓₓ = Diagonal((σₓ * ones(length(tₚ))).^2) + Clong

    return Cₓₓ
end

function gaussmarkovsolution(Eᵀ,y,σy,Cₓₓ)

    # Gauss-Markov solution method
    E = transpose(Eᵀ)
    x̃ = Cₓₓ*Eᵀ*((E*Cₓₓ*Eᵀ + σy^2)\y)

    # reduction in uncertainty
    P⁻ = Cₓₓ*Eᵀ*((E*Cₓₓ*Eᵀ + σy^2)\(E*Cₓₓ))
    P  = Cₓₓ - P⁻

    return x̃, P
end

"""
    function anomalymatrix(t₀,tf,tₚ)

    t₀: start of averaging interval
    tf: end of averaging interval
    tₚ: time grid
    M: matrix that takes anomaly relative to this interval 
"""
function anomalymatrix(t₀,tf,tₚ)

    igood = findall(x -> t₀ ≤ x ≤ tf,tₚ)
    nt = length(tₚ)
    M = Matrix{Float64}(undef,nt,nt)

    # for each row, subtract reference value
    for rr = 1:length(tₚ)
        M[rr,rr] = 1
        M[rr,igood] .-= 1 ./ length(igood)
    end
    return M
end

"""
    function trendmatrix(t₀,tf,tₚ)

    t₀: start of averaging interval
    tf: end of averaging interval
    tₚ: time grid
    F: matrix that takes trend over this interval 
"""
function trendmatrix(t₀,tf,tₚ)

    igood = findall(x -> t₀ ≤ x ≤ tf,tₚ)
    nt = length(tₚ)
    ngood = length(igood)

    # least-squares problem
    #G = Matrix{Float64}(undef,ngood,nt)
    G = zeros(ngood,nt)
    for gg = 1:ngood
        G[gg,igood[gg]] = 1.0 
    end

    E = Matrix{Float64}(undef,ngood,2)
    E[:,1] = tₚ[igood]
    E[:,2] = ones(ngood)

    #Ê = (transpose(G)*G)\(transpose(G)*E)
    Ê = transpose(G)*E

    Ffull = (transpose(Ê)*Ê)\(transpose(Ê))

    F = Matrix{Float64}(undef,1,nt)
    F[1,:] = Ffull[1,:]
    return F
end

"""
    function magnitude(xₚ,tₚ)

    make diagnostic: size of Δ between max min, assuming timing doesn't change
     Magnitude diagnostic
     make diagnostic: size of Δ between max min, assuming timing doesn't change
"""
function magnitude(xₚ,tₚ)

    # max, min surface signal
    tmax = tₚ[findmax(xₚ)[2]]
    tmin = tₚ[findmin(xₚ)[2]]

    # make diagnostic: size of Δ between max min, assuming timing doesn't change

    imax = findall(x -> tmax - 100 ≤ x ≤ tmax + 100, tₚ)
    imin = findall(x -> tmin - 100 ≤ x ≤ tmin + 100, tₚ)

    # Magnitude diagnostic
    M = zeros(1,length(tₚ))
    M[imax] .= 1/length(imax)
    M[imin] .= -1/length(imax)

    #Mdacp = zeros(1,length(tₚ))
    #Mdacp[imax] .= 1/length(imax)

    return mag = M*xₚ
    
end

"""
    function propagate(M,x̃,P)
"""
function propagate(M,x,P)
    # diagnostic
    d = dot(M,x)
    σd = √(M*(P*transpose(M)))[1]
    return d,σd
end

end
