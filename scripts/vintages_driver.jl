#= Overall driver to organize scripts

Gives a sense of the order of each script
=#

#=
Activate project environment
=#
using Revise
using DrWatson # must be available in default environment
@quickactivate "NobleGasRelic"

#=
Numerical calculations section
=#

using Revise
using DrWatson

# calculations for global distribution of vintages
include(scriptsdir("vintages_calculate.jl"))

# Gives vintage information needed for inversion
# calls vintages_calculate if necessary
include(scriptsdir("vintages_table_NPACvSPAC.jl"))

# invert for sea level pressure
include(scriptsdir("invert_sixvintages_NPACminusSPAC.jl"))

# are the SLP differences significant?
include(scriptsdir("vintage_SLP_differences.jl"))

#=
Diagostics Section
=#

# make planview figures of vintages
include(scriptsdir("vintages_planviews.jl"))

# make section figures of vintages
include(scriptsdir("vintages_sections.jl"))

# Look at age distributions
include(scriptsdir("diagnose_responses_NPACvSPAC.jl"))

# plot output from inversion
include(scriptsdir("plot_sixvintages_NPACminusSPAC.jl"))
