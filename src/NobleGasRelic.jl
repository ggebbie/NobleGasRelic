module NobleGasRelic

#using PyCall
#using PyPlot
using DrWatson, TMI, TMItransient, Interpolations
using LinearAlgebra
using GGplot

export vintages_planview, vintages_section,
    agedistribution,
    taudeltaresponse, compare_deltaresponses,
    priorcovariance,
    gaussmarkovsolution, anomalymatrix, magnitude,
    trendmatrix,
    propagate, vintage_atloc, diagnose_deltaresponse,
    define_vintages, vintages_longnames,
    vintages_longnameslabel

t_today = 2022
# each interval is 500 years
define_vintages(t_today) =  OrderedDict(:MOD => (1800, t_today),
                :LIA => (1300,1800),
                :MCA => (800,1300),
                :DACP => (300,800),
                #:DACP2 => (550,650),
                 :RWP => (-200,300),
                 :preRWP => (-Inf,-200))

vintages_longnames() = Dict(:MOD => "Modern Warming",
                :LIA => "Little Ice Age",
                :MCA => "Medieval Climate Anomaly",
                :DACP => "Dark Ages Cold Period",
                #:DACP2 => "Dark Ages Cold Period (strict)",
                :RWP => "Roman Warm Period",
                :preRWP => "Pre-Roman Warm Period")

function vintages_longnameslabel(longname)

    # add dates to longname
    longnamelabel = Dict{Symbol,String}()
    for (kk,vv) in longname
        longnamelabel[kk] = longname[kk]*" "*string(tinterval[kk])*" CE"
    end

    return longnamelabel
end

function vintages_planview(params)

    @unpack vintage, depth, tinterval, longnamelabel = params 

    # doing this every time, not so efficient
    
    Δ,τ = read_stepresponse()
    println(size(Δ))
    println(size(τ))
    println(tinterval[vintage][1],tinterval[vintage][2])
    
    g = vintagedistribution(tinterval[vintage][1],tinterval[vintage][2],Δ,τ,interp="linear")

    # get the meta-data correct on the output.
    gvintage = Field(g.tracer,g.γ,vintage,longnamelabel[vintage],"seawater mass fraction []")
    # save g to file if it hasn't been done before.
    if isapprox(depth,2000) # kludge to not write twice
        println("writing",depth,vintage)
        !isdir(DrWatson.datadir()) && mkdir(DrWatson.datadir())
        writefield(DrWatson.datadir("vintages_TMI_4x4_2012.nc"),gvintage)
    end
    
    froot = plotsdir(savename("TMI_4x4_2012",params,"png",accesses=["vintage","depth"]))
    println(froot)

    tlabel = "Vintage: "* longnamelabel[vintage] * ", depth="*string(depth)*"m"
    println(tlabel)
    #fname = "vintage_"*string(v)*"_"*string(depth)*"m.png"

    lims = vcat(collect(0:5:50),100)
    # Plan view plots
    plotfname = plotsdir("vintage.png")

    #GGplot.planviewplotcartopy(100g, depth, lims, titlelabel=tlabel)
    GGplot.planviewplotcartopy(100g, depth, lims, titlelabel=tlabel,fname=plotfname,cenlon=-160.0) 
    mv(plotsdir("vintage.png"),plotsdir(froot),force=true)

end

function vintages_section(params)

    @unpack vintage, lon, tinterval, longnamelabel = params 

    lims = vcat(collect(0:5:50),100)
    # doing this every time, not so efficient
    Δ,τ = read_stepresponse()
    local g = vintagedistribution(tinterval[vintage][1],tinterval[vintage][2],Δ,τ)
    
    froot = plotsdir(savename("TMI_4x4_2012",params,"png",accesses=["vintage","lon"]))
    println(froot)

    tlabel = "Vintage: "* longnamelabel[vintage] * ", lon="*string(lon)*"E"
    println(tlabel)

    sectionplot(100g, lon, lims, titlelabel=tlabel,fname=plotsdir(froot)) 

end

"""
    function diagnose_deltaresponse(loc)

    Input: one location

    Side-effect: plot of delta response and heaviside response
"""
function diagnose_deltaresponse(loc)
    g = agedistribution(loc)
    tg = taudeltaresponse()

    if loc[2] > 0
        leglabel = string((loc[2]))*"°N, "*string(360-loc[1])*"°W, "*string(round(loc[3]/1000,sigdigits=2))*" km"
    else
        leglabel = string((-loc[2]))*"°S, "*string(360-loc[1])*"°W, "*string(round(loc[3]/1000,sigdigits=2))*" km"
    end

    # UPDATE THIS TO PLOTS.JL
    # PyPlot version, not currently showing
    # figure(2)
    # clf()
    # line1, = PyPlot.plot(tg,g,"black",label=leglabel[1])
    # grid("true")
    # xlabel("Lag, τ [yr]")
    # ylabel("mass fraction per yr [1/yr]")
    # legend()
    # PyPlot.savefig(plotsdir("deltaresponse_at_loc.png"))

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

    leglabel = Vector{String}(undef,2)
    for ii = 1:2
        if loc[ii][2] > 0
            leglabel[ii] = string((loc[ii][2]))*"°N, "*string(360-loc[ii][1])*"°W, "*string(round(loc[ii][3]/1000,sigdigits=2))*" km"
        else
            leglabel[ii] = string((-loc[ii][2]))*"°S, "*string(360-loc[ii][1])*"°W, "*string(round(loc[ii][3]/1000,sigdigits=2))*" km"
        end           
    end

    # UPDATE TO USE PLOTS.JL
    # PyPlot version, not currently showing
    # figure(2)
    # clf()
    # line1, = PyPlot.plot(tg,g[1],"black",label=leglabel[1])
    # line2, = PyPlot.plot(tg,g[2],"red",label=leglabel[2])
    # line3, = PyPlot.plot(tg,Δg,"green",label="Δ")
    # grid("true")
    # xlabel("Lag, τ [yr]")
    # ylabel("mass fraction per yr [1/yr]")
    # legend()
    # PyPlot.savefig(plotsdir("deltaresponse_NPACvSPAC.png"))

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

# UPDATE TO USE BLUES.JL
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

"""
function vintage_atloc(vname,loc)

Compute fractional contribution of vintage `vname`
at location/s `loc`
"""
function vintage_atloc(vname,loc)

    gname = datadir("vintages_TMI_4x4_2012.nc")
    γ = Grid(TMI.pkgdatadir("TMI_modern_90x45x33_GH10_GH12.nc"))
    gvintage = readfield(gname,vname,γ)

    # get weighted interpolation indices
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,2)

    for (i,v) in enumerate(loc)
        wis[i] = interpindex(v,γ)
    end

    return gloc = observe(gvintage,wis,γ) 
end
    
end
