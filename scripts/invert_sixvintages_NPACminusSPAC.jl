#using DrWatson
#@quickactivate "NobleGasRelic"

using Revise
using NobleGasRelic
using DrWatson
using DataFrames
using OrderedCollections
using CSV
using BLUEs
using DimensionalData
using DimensionalData:@dim
using Unitful
using UnitfulLinearAlgebra
using Measurements

include(srcdir("config_vintages.jl"));

# HERE DEEP NORTH PACIFIC VS. DEEP SOUTH PACIFIC 
n = 2
loc = Vector{Tuple}(undef,n)
loc[1] = (360-152,35,3500m) # North Pacific
loc[2] = (360-152,-10,3500m) # South Pacific

# read output of vintages_table_NPACvSPAC into DataFrame
csvinput = datadir("sixvintages_"*TMIversion*".csv")

# Do computations if necessary
!isfile(csvinput) && include(scriptsdir("vintages_table_NPACvSPAC.jl"))

df = DataFrame(CSV.File(csvinput))

# extract list of vintages
iloc = [4,5] # indices of vintages in DataFrame
vintages = df[:,:Vintage]
locnames = names(df)[iloc]
    
# associate vintages with a Dimension
#@dim YearCE "years Common Era"
@dim Vintage "vintage"
@dim InteriorLocation "interior location"

# observations have units of mbar
mbar = u"mbar"
M = 1 # just deal with the difference, otherwise length(iloc) # number of obs
urange1 = fill(mbar,M)

# solution also has units of mbar
N = length(vintages)
udomain = fill(mbar,N)

# manually determine that units on matrix are percent
# please encode in CSV, probably in the header for portability
percent = u"percent"
𝐌 = uconvert.(NoUnits,(Matrix(df)[:,iloc])percent)
𝐦 = Matrix(transpose(𝐌[:,2]-𝐌[:,1]))

E = UnitfulDimMatrix(ustrip.(𝐦),urange1,udomain,dims=(InteriorLocation([:NPACminusSPAC]),Vintage(vintages)))

iszero(sum(E)) && println("not normalized correctly")

Δp★ = (2.8 ± 0.4)mbar # observed value with 1σ uncertainty
yr = u"yr"

# other fixed parameters (could be `const`)
scentury = 4mbar/100yr
σref = (0.0001)mbar # error in pre-industrial reference
σSLP₀ = (10.0)mbar  # from existing SLP gradients and historical model simulations

# Solve three cases.
cases = ("min_trend","min_variance","min_trend_variance")
for case in cases

    df = DataFrame(CSV.File(csvinput)) # reload: seems to be having a problem without this.
    
    # make a covariance matrix that penalizes differences
    # greater than 1 mbar/century
    if case == "min_trend"
        # strip units from functions
        S⁻ = invcovariance_temporalsmoothness(tinterval,scentury)
        S⁻ += invcovariance_preindustrialmean(vintage,σref)
    elseif case == "min_variance"
        S⁻ = invcovariance_minenergy(vintage,σSLP₀) # 20 dbar magnitude of SLP changes
        S⁻ += invcovariance_preindustrialmean(vintage,σref)
    elseif case == "min_trend_variance"
        S⁻ = invcovariance_temporalsmoothness(tinterval,scentury)
        S⁻ += invcovariance_minenergy(vintage,σSLP₀)
        S⁻ += invcovariance_preindustrialmean(vintage,σref)
    else
        error("no case chosen")
    end
    # add units at the end
    Cxxdims = (last(dims(E)),last(dims(E)))
    Cₓₓ = UnitfulDimMatrix(inv(S⁻),unitdomain(E),unitdomain(E).^-1,dims=Cxxdims)

    Cnndims = (first(dims(E)),first(dims(E)))
    Cnn = UnitfulDimMatrix([ustrip(Measurements.uncertainty(Δp★).^2);;],unitrange(E),unitrange(E).^-1,dims=Cnndims)

    y = UnitfulDimMatrix([ustrip(Measurements.value(Δp★));],unitrange(E),dims=(first(dims(E))))
    x₀ = UnitfulDimMatrix(zeros(N),unitdomain(E),dims=(last(dims(E))))
    
    uproblem = UnderdeterminedProblem(y,E,Cnn,Cₓₓ,x₀)
    x̃ = solve(uproblem)
    
    col5 = "SLP Anomaly [mbar]"
    col6 = "SLP Error [mbar]"

    insertcols!(df, col5 => [round(x̃.v[At(string(vv))],digits=1) for vv in vintage])

    #upstream bug in BLUEs.jl: dimensions are lost
    # assume order is correct
    #insertcols!(df, col6 => [round(x̃.σ[At(string(vv))],digits=1) for vv in vintage])
    insertcols!(df, col6 => [round(ustrip(x̃.σ[vv]),digits=1) for vv in eachindex(vintage)])
    CSV.write(datadir("sixvintages_"*case*".csv"),df)

    # # make a plot
    # t̄ = midtime(tinterval)
    # if case == "min_trend"
    #     plot(collect(values(t̄)),df[:,6],ribbon=df[:,7],label=case,xlabel="calendar years",ylabel="SLP anomaly [dbar]")
    # else
    #     plot!(collect(values(t̄)),df[:,6],ribbon=df[:,7],label=case,xlabel="years [CE]",ylabel="SLP anomaly [dbar]",xticks=(-450:250:2022))
        
    # end
    # savefig(plotsdir("SLP_CommonEra_combined.pdf"))
end
