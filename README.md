# PacificNobleGasRelic

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> PacificNobleGasRelic

It is authored by G Jake Gebbie.

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


# How this scientific project was constructed

Make bare repository on GitHub (no readme, .gitignore, or anything else). 

In your default Julia environment (i.e, @v1.7), `add DrWatson` and `using DrWatson`.

Change to a place in your file system where you can make a new directory, i.e., `cd("GOODPLACE")` at the Julia REPL

Initialize DrWatson scientific project in Julia REPL \
`initialize_project("PROJECTNAME"; authors="G Jake Gebbie)` 
PROJECTNAME doesn't need .jl at the end because this is not a Julia software package. DrWatson creates a README and makes a git repository. Sadly it doesn't make a license file. 

Add GitHub as the remote repository: 

Via command line with GitHub SSH keys and default branch name = main (not master): \
`git remote add origin git@github.com:USERNAME/PROJECTNAME.git` \
`git branch --set-upstream-to=origin/main main`

Via Magit \
`M a` ;; remote add origin\
then Magit asks if the upstream should be set to `origin/main`, answer yes.

Think about adding a License at this point.
