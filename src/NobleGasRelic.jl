module NobleGasRelic

using PyCall, PyPlot, DrWatson, TMI, TMItransient, Interpolations

export vintages_planview, vintages_section, agedistribution, taudeltaresponse

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
    local g = vintagedistribution(tinterval[vintage][1],tinterval[vintage][2],Δ,τ)

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
    #fname = "vintage_"*string(v)*"_"*string(depth)*"m.png"

    #println(size(g.tracer))
    sectionplot(100g, lon, lims; titlelabel=tlabel) 

    savefig(plotsdir(froot))

end

"""
    function agedistribution

    age distribution refers to distribution over lags, not space
    sometimes called a transit time distribution
"""
function agedistribution(loc)

    Δ,τ = read_stepresponse()
    ncdf = length(τ)
    npdf = ncdf - 1

    # get weighted interpolation indices
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,1)
    #[wis[xx] = interpindex(loc[xx],Δ[xx].γ) for xx in 1]
    wis[1] = interpindex(loc,Δ[1].γ)
    
    Δloc = Vector{Float64}(undef,ncdf)
    for tt in 1:ncdf
        Δloc[tt] = observe(Δ[tt],wis,Δ[tt].γ)[1] # kludge to convert to scalar
    end
    tg, g = deltaresponse(Δloc,τ)
    
    return g
    
end

"""
    function deltaresponse

    Take CDF and turn it into PDF
"""
function deltaresponse(Δ,τΔ)

    # hard coded annual resolution, easy because don't have to normalize
    τmax = maximum(τΔ)
    tgedge = 0:floor(τmax)
    tg = 0.5:floor(τmax)

    #    interp_linear = LinearInterpolation(τΔ, Δ)
    # Δhires = interp_linear(tgedge)

    itp = interpolate(τΔ, Δ, FritschCarlsonMonotonicInterpolation())
    Δhires = itp.(tgedge)
    g = diff(Δhires)
        
    return tg, g
end

"""
    function taudeltaresponse

    Take CDF and turn it into PDF, get lag timescale
"""
function taudeltaresponse()

    Δ,τ = read_stepresponse()

    # hard coded annual resolution, easy because don't have to normalize
    τmax = maximum(τ)
    tgedge = 0:floor(τmax)
    tg = 0.5:floor(τmax)
    return tg
end


end
