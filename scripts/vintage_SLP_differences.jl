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
using JLD2

include(srcdir("config_vintages.jl"));

# Solve three cases.
cases = ["min_trend","min_variance","min_trend_variance"]
vintage1 = vintage
vintage2 = vintage
    
params = @strdict cases vintage1 vintage2
dicts = dict_list(params) # 108 cases to test

## Stopped here
for case in cases

jldinput = datadir("sixvintages_"*case*".jld2")

load(jldinput)

# Do computations if necessary UPDATE THIS
!isfile(csvinput) && include(scriptsdir("vintages_table_NPACvSPAC.jl"))

    # NOT SURE IF THIS IS NECESSARy
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


    df = DataFrame(CSV.File(csvinput)) # reload: seems to be having a problem without this.
    
    # make a covariance matrix that penalizes differences
    # greater than 1 mbar/century
    if case == "min_trend"

    elseif case == "min_variance"

    elseif case == "min_trend_variance"

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

    # Save full estimate in JLD2 format.
    jldsave(datadir("sixvintages_"*case*".jld2"),x̃)

end
    jldsave(datadir("sixvintages_"*case*".jld2"),x̃)
