#=
Calculate vintage percentages as seen in planviews and Pacific sections
=# 
using Revise
using NobleGasRelic
using DrWatson
using LinearAlgebra
using TMI
using DataFrames
using Interpolations
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

# lon must match grid exactly or else "dropped dims" error
lon = [162, 210]

params = @strdict vintage depth tinterval longnamelabel
dicts = dict_list(params)

# avoid error with NetCDF key already existing by
# moving existing file to SAVE version
isfile(datadir("vintages_TMI_4x4_2012.nc")) &&
    mv(datadir("vintages_TMI_4x4_2012.nc"),
       datadir("vintages_TMI_4x4_2012_SAVE.nc"),force=true)

# planviews
map(vintages_planview,dicts)

## sections
params = @strdict vintage lon tinterval longnamelabel
dicts = dict_list(params)
map(vintages_section,dicts)
