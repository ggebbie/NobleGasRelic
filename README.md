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

- `vintages_table_NPAC_SPAC.jl`: compute water-mass percentages for six vintages
- `vintages_calculate_and_plot.jl`: Calculate vintage percentages and create planviews and Pacific sections
- `diagnose_responses.jl`: understand, plots, and compare two TTDs
- `invert_sixvintages.jl`: invert for time history of each vintage, plot 3 different cases
- `invert_sixvintages_oneplot.jl`: invert for time history of each vintage, plot 3 different cases on one plot

Deprecated
- `vintages_table.jl`
- `vintages_driver.jl` original driver, deprecated

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
