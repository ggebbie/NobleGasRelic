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
df = DataFrame(CSV.File(csvinput))

E = Matrix(df)[:,4:5]
ΔE = transpose(E[:,2]-E[:,1])/100

# find time difference between vintages
#t̄ = midtime(tinterval)

# make a covariance matrix that penalizes differences
# greater than 1 mbar/century
S⁻ = covariance_temporalsmoothness(tinterval)
Cₓₓ = inv(S⁻)

# Solve it.
ΔNe = 2.8#mbar # mbar
σΔNe = 0.4#mbar # mbar

xtmp,P = gaussmarkovsolution(transpose(ΔE),ΔNe,σΔNe,Cₓₓ)

#xtmp = (transpose(ΔE)*W⁻*ΔE + S⁻) \ (transpose(ΔE)*W⁻*y)
ỹ = (transpose(E)*xtmp)/100
ñ = ΔE*xtmp - ΔNe

x̃ = Dict{Symbol,Float64}()
σx̃ = Dict{Symbol,Float64}()
for (mm,ii) in enumerate(vintage)
    x̃[ii] = xtmp[mm]
    σx̃[ii] = √P[mm,mm]
end

col5 = "SLP Anomaly [mbar]"
col6 = "SLP Error [mbar]"

insertcols!(df, col5 => [round(x̃[vv],digits=1) for vv in vintage])
insertcols!(df, col6 => [round(σx̃[vv],digits=1) for vv in vintage])
CSV.write(datadir("sixvintages_mintrend.csv"),df)
