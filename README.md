# PacificNobleGasRelic

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> PacificNobleGasRelic

It is authored by G Jake Gebbie.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

# Scripts

- `vintages_calculate_and_plot.jl`: Calculate vintage percentages and create planviews and Pacific sections
- `vintages_table_NPACvSPAC.jl`: compute water-mass percentages for six vintages at two locations (NPAC, SPAC)
- `diagnose_responses_NPACvSPAC.jl`: understand, plots, and compare two TTDs, one in NPAC, one in SPAC
- `invert_sixvintages.jl`: invert for time history of each vintage, plot 3 different cases
- `invert_sixvintages_oneplot.jl`: invert for time history of each vintage, plot 3 different cases on one plot

Deprecated
- `vintages_table.jl`
- `vintages_driver.jl` original driver, deprecated
