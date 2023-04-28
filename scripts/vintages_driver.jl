#= Overall driver to organize scripts

Gives a sense of the order of each script
=#

# calls vintages_calculate if necessary
include(scriptsdir("vintages_table_NPACvSPAC.jl"))

# make figures of vintages
include(scriptsdir("vintages_plot_all.jl"))

include(scriptsdir("diagnose_responses_NPACvSPAC.jl"))
