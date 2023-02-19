using Revise
using NobleGasRelic
using DrWatson
#using PyPlot
using LinearAlgebra
using TMI
using DataFrames
using Interpolations
#using PlotlyJS
#using Plots
using OrderedCollections
using CSV

t_today = 2022
# each interval is 500 years
tinterval = define_vintages(t_today)
longname = vintages_longnames()
# add dates to longname
longnamelabel = vintages_longnameslabel(longname,tinterval)
vintage = collect(keys(tinterval))
depth = collect(2000:500:4000)

# HERE DEEP NORTH PACIFIC VS. DEEP SOUTH PACIFIC 
# get 2 age distributions (TTDs), ultimately we have info about their difference
n = 2
loc = Vector{Tuple}(undef,n)
loc[1] = (360-152,35,3500) # North Pacific
loc[2] = (360-152,-10,3500) # South Pacific

col1 = "Vintage"
col2 = "Vintage Name"
col2b = "Years CE"
col3 = "Southern Region"
col4 = "Northern Region"
col5 = "SLP Anomaly [mbar]"
col6 = "SLP Error [mbar]"

defs = Dict(col1 => vintage,
            col2 => [longnamelabel[vv] for vv in vintage],
            col2b => [tinterval[vv] for vv in vintage])
df = DataFrame(defs)

gnorth = Dict{Symbol,Float64}()
gsouth = Dict{Symbol,Float64}()
for vv in vintage
    println(vv)
    gtmp = vintage_atloc(vv,loc)
    gnorth[vv] = gtmp[1]
    gsouth[vv] = gtmp[2]
end

insertcols!(df, col3 => [round(100gsouth[vv],digits=1) for vv in vintage])
insertcols!(df, col4 => [round(100gnorth[vv],digits=1) for vv in vintage])
