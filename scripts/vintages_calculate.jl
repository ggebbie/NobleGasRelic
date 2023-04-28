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
using OrderedCollections
using CSV

include(srcdir("config_vintages.jl"))

params = @strdict vintage tinterval longnamelabel
#params = @strdict vintage depth tinterval longnamelabel
dicts = dict_list(params)

# avoid error with NetCDF key already existing by
# moving existing file to SAVE version
isfile(datadir("vintages_TMI_4x4_2012.nc")) &&
    mv(datadir("vintages_TMI_4x4_2012.nc"),
       datadir("vintages_TMI_4x4_2012_SAVE.nc"),force=true)

map(vintages_calculate,dicts)
