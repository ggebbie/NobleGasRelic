# NobleGasRelic

## Algorithms to invert modern oceanic noble gas anomalies for the history of sea level pressure

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://github.com/ggebbie/NobleGasRelic/papers/Methods.pdf)

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> NobleGasRelic

It is authored by G Jake Gebbie <ggebbie@whoi.edu>.

# Scientific context

The scientific context for the algorithms is given in `papers/Methods.pdf`. 

# Reproducibility

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently. You could clone from GitHub:
   ```sh
   git clone https://github.com/ggebbie/NobleGasRelic
   ```

1. Download julia. I recommend the `juliaup` program. On Linux or MacOSX, try the following to install Julia 1.9.0 release candidate 3:
```sh
curl -fsSL https://install.julialang.org | sh
juliaup add 1.9.0-rc3
juliaup default 1.9.0-rc3
```

2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

# Scripts

You may start from any script of interest. Necessary prerequisite calculations will automatically be done. After activating the project environment you may run a script with the following command at the repl: `include("sample_script_name.jl")`

- `vintages_driver.jl` : An overview of the scripts with an optimal run order.

## Numerical calculations

- `vintages_calculate_and_plot.jl`: Calculate vintage percentages and distributions

- `vintages_table_NPACvSPAC.jl`: compute water-mass percentages for six vintages at two locations (NPAC, SPAC)

- `invert_sixvintages_NPACminusSPAC.jl`: invert for time history of sea-level pressure of each vintage, plot 3 different cases

- `vintages_SLP_differences.jl` : are the SLP differences significant?

## Diagnostics and figures

- `diagnose_responses_NPACvSPAC.jl`: understand, plots, and compare two TTDs, one in NPAC, one in SPAC, look at age distributions

- `vintages_planviews.jl`: make planview figures of vintages

- `vintages_sections.jl`: make section figures of vintages

- `plot_sixvintages_NPACminusSPAC.jl` : plot output from inversion

