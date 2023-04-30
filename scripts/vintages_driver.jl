#= Overall driver to organize scripts

Gives a sense of the order of each script
=#

#=
Numerical calculations section

=#

# calculations for global distribution of vintages
include(scriptsdir("vintages_calculate.jl"))

# Gives vintage information needed for inversion
# calls vintages_calculate if necessary
include(scriptsdir("vintages_table_NPACvSPAC.jl"))

## HERE: add inversion routines


#=
Diagostics Section
=#

# make planview figures of vintages
include(scriptsdir("vintages_planviews.jl"))

# make section figures of vintages
include(scriptsdir("vintages_sections.jl"))

# Look at age distributions
include(scriptsdir("diagnose_responses_NPACvSPAC.jl"))
