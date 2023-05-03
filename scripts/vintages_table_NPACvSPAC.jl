#=
`vintages_table_NPAC_SPAC.jl`: compute water-mass percentages for six vintages at 2 locations
=#
using Revise
using NobleGasRelic
using DrWatson
using OrderedCollections
using CSV

include(srcdir("config_vintages.jl"));

# HERE DEEP NORTH PACIFIC VS. DEEP SOUTH PACIFIC 
# get 2 age distributions (TTDs), ultimately we have info about their difference
n = 2
loc = Vector{Tuple}(undef,n)
loc[1] = (360-152,35,3500) # North Pacific
loc[2] = (360-152,-10,3500) # South Pacific

!isfile(datadir("vintages_TMI_4x4_2012.nc")) && include("vintages_calculate.jl")

df = vintages_table(loc,vintage,tinterval,longnamelabel)
CSV.write(datadir("sixvintages_"*TMIversion*".csv"),df)
