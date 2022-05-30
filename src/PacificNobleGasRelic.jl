module PacificNobleGasRelic

using PyCall, PyPlot, DrWatson, TMI, TMItransient, Interpolations

export vintages_planview, vintages_section

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
    proj0 = PacificNobleGasRelic.cartopy.crs.PlateCarree()
    proj = PacificNobleGasRelic.cartopy.crs.PlateCarree(central_longitude=cenlon)
    ax = fig.add_subplot(projection = proj)
    ax.set_global()
    #ax.coastlines()

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

    froot = plotsdir(savename("TMI_4x4_2012",params,"png",accesses=["vintage","depth"]))
    println(froot)

    tlabel = "Vintage: "* longname[vintage] * " " * string(tinterval[vintage]) * " CE, depth="*string(depth)*"m"
    println(tlabel)
    #fname = "vintage_"*string(v)*"_"*string(depth)*"m.png"

    local g = vintagedistribution(tinterval[vintage][1],tinterval[vintage][2],Δ,τ)

    lims = vcat(collect(0:5:50),100)
    # Plan view plots
    PacificNobleGasRelic.planviewplotcartopy(100g, depth, lims, titlelabel=tlabel)
    mv(plotsdir("vintage.png"),plotsdir(froot),force=true)

end

function vintages_section(params)

    @unpack vintage, lon, tinterval, longname = params 

    lims = vcat(collect(0:5:50),100)
    # doing this every time, not so efficient
    Δ,τ = read_stepresponse()

    
    froot = plotsdir(savename("TMI_4x4_2012",params,"png",accesses=["vintage","lon"]))
    println(froot)

    tlabel = "Vintage: "* longname[vintage] * " " * string(tinterval[vintage]) * " CE, lon="*string(lon)*"E"
    println(tlabel)
    #fname = "vintage_"*string(v)*"_"*string(depth)*"m.png"

    local g = vintagedistribution(tinterval[vintage][1],tinterval[vintage][2],Δ,τ)

    #println(size(g.tracer))
    sectionplot(100g, lon, lims; titlelabel=tlabel) 

    savefig(plotsdir(froot))

end

"""
    function agedistribution

    age distribution refers to distribution over lags, not space
    sometimes called a transit time distribution
"""
function agedistribution(Δ,τΔ,loc)
    Δ,τ = read_stepresponse()

    # get weighted interpolation indices
    wis1= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,2)
    [wis[i] = interpindex(locs[i],Δ[1].γ) for i in 1]

    Δloc = Vector(undef,length(Δ))
    for tt in 1:length(Δ)
        Δloc[tt] = observe(Δ[tt],wis,Δ[tt].γ)
    end

    g = Vector{Any}(undef,2)
    for nn = 1:2
        temp = Vector{Float64}(undef,length(Δloc))
        for mm in 1:length(Δloc)
            temp[mm] = Δloc[mm][nn]
        end
        tg, g[nn] = deltaresponse(tmp,τ)
    end
    
    return nothing
    
end

"""
    function deltaresponse

    Take CDF and turn it into PDF
"""
function deltaresponse(Δ,τΔ)

    # hard coded annual resolution, easy because don't have to normalize
    tgedge = 0:2000
    tg = 0.5:2000

    interp_linear = LinearInterpolation(τΔ, Δ)
    Δhires = interp_linear(tgedge)
    g = diff(Δhires)
        
    return tg, g
end


end
