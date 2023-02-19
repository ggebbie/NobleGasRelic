using Revise
using NobleGasRelic
using DrWatson
#using PyPlot
using LinearAlgebra
using TMI
using DataFrames
using Interpolations
#using PlotlyJS
using Plots
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

# a script to understand TTDs etc.

# get 2 age distributions (TTDs), ultimately we have info about their difference
n = 2
loc = Vector{Tuple}(undef,n)
loc[1] = (360-152,35,3500) # North Pacific
loc[2] = (360-152,-10,3500) # South Pacific
compare_deltaresponses(loc)

fracs = Dict{Symbol,Float64}()
fracs[:MOD] = vintage_atloc(:MOD,loc)[1]

# N.Pac: why more RWP than DACP water?
#diagnose_deltaresponse(loc,vintage,tinterval) # seems to repeat compare_deltaresponses
