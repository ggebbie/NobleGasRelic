#using DrWatson
#@quickactivate "NobleGasRelic"

using Revise
using NobleGasRelic
using DrWatson
using DataFrames
using CSV
using BLUEs
using Unitful
using Plots

include(srcdir("config_vintages.jl"));

# should actually be computed from CSV file
t̄ = midtime(tinterval)

# Solve three cases.
cases = ("min_trend","min_variance","min_trend_variance")
for case in cases

    csvfile = datadir("sixvintages_"*case*".csv")

    # Do computations if necessary
    !isfile(csvfile) && include(scriptsdir("invert_sixvintages_NPACminusSPAC.jl"))

    df = DataFrame(CSV.File(csvfile))
    # make a plot

    if case == "min_trend"
        # should read units from CSV
        plot(ustrip.(collect(values(t̄))),ustrip.(df[:,6]mbar),ribbon=ustrip.(df[:,7]mbar),label=case,xlabel="calendar years [CE]",ylabel="SLP anomaly [mbar]",xticks=(-450:250:2022))
        #plot(collect(values(t̄)),df[:,6]mbar,ribbon=df[:,7]mbar,label=case,xlabel="calendar years",ylabel="SLP anomaly")
    else
        plot!(ustrip.(collect(values(t̄))),df[:,6],ribbon=df[:,7],label=case)

        #plot!(ustrip.(collect(values(t̄))),df[:,6],ribbon=df[:,7],label=case,xlabel="years [CE]",ylabel="SLP anomaly [mbar]",xticks=(-450:250:2022))
        
               # plot!(collect(values(t̄)),df[:,6],ribbon=df[:,7],label=case,xlabel="years [CE]",ylabel="SLP anomaly",xticks=(-450:250:2022))
        
    end
    savefig(plotsdir("SLP_CommonEra_3cases.pdf"))
end
