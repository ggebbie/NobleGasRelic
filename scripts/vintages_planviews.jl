#=
Calculate vintage percentages as seen in planviews and Pacific sections
=# 
using Revise
using NobleGasRelic
using DrWatson
#using LinearAlgebra
#using TMI
#using DataFrames
#using Interpolations
#using OrderedCollections
#using CSV

include(srcdir("config_vintages.jl"));

# if NetCDF file doesn't exist, use this as a proxy that calculations haven't been done.
# in reality, only jld2 files are actually needed
!isfile(datadir("vintages_TMI_4x4_2012.nc")) && include(scriptsdir("vintages_calculate.jl"))

# planviews
params = @strdict vintage depth tinterval longnamelabel
dicts = dict_list(params)

map(vintages_planview,dicts)
