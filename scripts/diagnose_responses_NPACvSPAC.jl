# a script to understand and compare two TTDs
using Revise
using NobleGasRelic
using DrWatson
using LinearAlgebra
using TMI
using DataFrames
using Interpolations
#using PlotlyJS
using Plots
using OrderedCollections
using CSV

# get 2 age distributions (TTDs), ultimately we have info about their difference
n = 2
loc = Vector{Tuple}(undef,n)
loc[1] = (360-152,35,3500) # North Pacific
loc[2] = (360-152,-10,3500) # South Pacific
compare_deltaresponses(loc)

## extra lost code

#fracs = Dict{Symbol,Float64}()
#fracs[:MOD] = vintage_atloc(:MOD,loc)[1]

# N.Pac: why more RWP than DACP water?
#diagnose_deltaresponse(loc,vintage,tinterval) # seems to repeat compare_deltaresponses
