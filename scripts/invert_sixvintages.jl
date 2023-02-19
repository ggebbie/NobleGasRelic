using Revise
using NobleGasRelic
using DrWatson
using LinearAlgebra
using TMI
using DataFrames
using Interpolations
using OrderedCollections
using CSV

TMIversion = "TMI_4x4x33"
t_today = 2022
# each interval is 500 years
tinterval = define_vintages(t_today)
longname = vintages_longnames()
# add dates to longname
longnamelabel = vintages_longnameslabel(longname,tinterval)
vintage = collect(keys(tinterval))
depth = collect(2000:500:4000)

# HERE DEEP NORTH PACIFIC VS. DEEP SOUTH PACIFIC 
n = 2
loc = Vector{Tuple}(undef,n)
loc[1] = (360-152,35,3500) # North Pacific
loc[2] = (360-152,-10,3500) # South Pacific

# read output of vintages_table_NPACvSPAC into DataFrame
csvinput = datadir("sixvintages_"*TMIversion*".csv")

E = Matrix(df)[:,4:5]
ΔE = transpose(E[:,2]-E[:,1])/100

# Solve it.
ΔNe = 2.8#mbar # mbar
σΔNe = 0.4#mbar # mbar

cases = ("min_trend","min_variance","min_trend_variance")
scentury = 4
referror = 0.0001
σSLP₀ = 10.0 #dbar
for case in cases
    df = DataFrame(CSV.File(csvinput))
    
    # make a covariance matrix that penalizes differences
    # greater than 1 mbar/century
    if case == "min_trend"
        S⁻ = invcovariance_temporalsmoothness(tinterval,scentury)
        S⁻ += invcovariance_preindustrialmean(vintage,referror)
    elseif case == "min_variance"
        S⁻ = invcovariance_minenergy(vintage,σSLP₀) # 20 dbar magnitude of SLP changes
        S⁻ += invcovariance_preindustrialmean(vintage,referror)
    elseif case == "min_trend_variance"
        S⁻ = invcovariance_temporalsmoothness(tinterval,scentury)
        S⁻ += invcovariance_minenergy(vintage,σSLP₀)
        S⁻ += invcovariance_preindustrialmean(vintage,referror)
    else
        error("no case chosen")
    end
    Cₓₓ = inv(S⁻)

    xtmp,P = gaussmarkovsolution(transpose(ΔE),ΔNe,σΔNe,Cₓₓ)
    
    #xtmp = (transpose(ΔE)*W⁻*ΔE + S⁻) \ (transpose(ΔE)*W⁻*y)
    ỹ = (transpose(E)*xtmp)/100
    ñ = ΔE*xtmp - ΔNe

    println("noise = ",ñ)
    
    x̃ = OrderedDict{Symbol,Float64}()
    σx̃ = OrderedDict{Symbol,Float64}()
    for (mm,ii) in enumerate(vintage)
        x̃[ii] = xtmp[mm]
        σx̃[ii] = √P[mm,mm]
    end

    col5 = "SLP Anomaly [mbar]"
    col6 = "SLP Error [mbar]"

    insertcols!(df, col5 => [round(x̃[vv],digits=1) for vv in vintage])
    insertcols!(df, col6 => [round(σx̃[vv],digits=1) for vv in vintage])
    #CSV.write(datadir("sixvintages_"*case*".csv"),df)

    # make a plot
    t̄ = midtime(tinterval)
    #plot(collect(values(t̄)),df[:,6],ribbon=df[:,7])
    plot(collect(values(t̄)),df[:,6],ribbon=df[:,7],label=case,xlabel="years [CE]",ylabel="SLP anomaly [dbar]",xticks=(-450:250:2022))
    savefig(plotsdir("SLP_CommonEra_"*case*".pdf"))

end
