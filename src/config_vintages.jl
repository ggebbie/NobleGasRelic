#= define the vintages

This information is used in multiple scripts
=#

const TMIversion = "TMI_4x4x33" # relax if run with other versions

# related to units
using Unitful
ENV["UNITFUL_FANCY_EXPONENTS"] = true
const yr = u"yr"
const m = u"m"
const mbar = u"mbar"

const t_today = 2022yr # calendar year

# each interval is 500 years
tinterval = define_vintages(t_today)

longname = vintages_longnames()

# add dates to longname
longnamelabel = vintages_longnameslabel(longname,tinterval)

vintage = collect(keys(tinterval))

depth = collect(2000:500:4000)m
