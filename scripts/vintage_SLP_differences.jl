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

# associate vintages with a Dimension
#@dim Vintage "vintage"
#@dim Diagnostic "vintage differences"

# Solve three cases.
cases = ["min_trend","min_variance","min_trend_variance"]
#vintage1 = vintage
#vintage2 = vintage
    
#params = @strdict vintage1 vintage2
#dicts = dict_list(params) # 36 permutations to test for each case

# Could use DrWatson, special naming, etc. below here.
for case in cases

    jldinput = datadir("sixvintages_"*case*".jld2")

    # Do computations if necessary UPDATE THIS
    !isfile(jldinput) && include(scriptsdir("invert_sixvintages_NPACminusSPAC.jl"))

    datain = load(jldinput)
    x̃ = datain["x̃"]  # seems awkward

    # make a Diagnostic or Difference matrix

    # diagnostic output has units of mbar
    mbar = u"mbar"
    M = 1 # just deal with the difference, otherwise length(iloc) # number of obs
    urange = fill(mbar,M)

    # solution also has units of mbar
    N = length(vintages)
    udomain = fill(mbar,N)
    # make a DimensionalData Array to store data?
    # or use CSV?
    Δp★ = DimArray(Matrix{Quantity}(undef,N,N),(Vintage(vintages),Vintage(vintages))); #no show

    for v1 in vintages
        for v2 in vintages
            D = UnitfulDimMatrix(zeros(M,N),urange,udomain,dims=(Diagnostic([:Δp★]),Vintage(vintages)))

            D[1,At(v1)] += 1
            D[1,At(v2)] -= 1

            Δp★[At(v1),At(v2)] = (D*x̃).x[1]
        end
    end

#    CSV.write(datadir("SLP_differences_"*case*".csv"),DimTable(DA))
    jldsave(datadir("SLP_differences_"*case*".jld2");Δp★)

    # save latexified table to textfile
    caseroot = "SLP_differences_"*case
    redirect_stdio(stdout=datadir(caseroot*".txt"), stderr=datadir(caseroot*"_err.txt")) do
        println(latexify(load(datadir(caseroot*".jld2"))["Δp★"]))
        #println(stderr, "hello error")
    end
end
# send output to REPL; loses header info unfortunately
latexify(load(datadir(caseroot*".jld2"))["Δp★"])

