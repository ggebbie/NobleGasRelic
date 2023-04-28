#=
Calculate vintage percentages as seen in planviews and Pacific sections
=# 
using Revise
using NobleGasRelic
using DrWatson
using OrderedCollections
using CSV

include(srcdir("config_vintages.jl"))

# if NetCDF file doesn't exist, use this as a proxy that calculations haven't been done.
# in reality, only jld2 files are actually needed
!isfile(datadir("vintages_TMI_4x4_2012.nc")) && include(scriptsdir("vintages_calculate.jl"))

## sections
# lon must match grid exactly or else "dropped dims" error
lon = [162, 210]

params = @strdict vintage lon tinterval longnamelabel
dicts = dict_list(params)
map(vintages_section,dicts)
