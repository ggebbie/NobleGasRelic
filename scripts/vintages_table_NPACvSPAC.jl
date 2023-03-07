#=
`vintages_table_NPAC_SPAC.jl`: compute water-mass percentages for six vintages
=#
using Revise
using NobleGasRelic
using DrWatson
using LinearAlgebra
using TMI
using DataFrames
using Interpolations
#using PlotlyJS
#using Plots
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
# get 2 age distributions (TTDs), ultimately we have info about their difference
n = 2
loc = Vector{Tuple}(undef,n)
loc[1] = (360-152,35,3500) # North Pacific
loc[2] = (360-152,-10,3500) # South Pacific

df = vintages_table(loc,vintage,tinterval,longnamelabel)
CSV.write(datadir("sixvintages_"*TMIversion*".csv"),df)
