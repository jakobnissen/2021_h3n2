# H3N2 project

Date of project initialization: 2021-01-25

## Setup and installation
* Julia version 1.6.0
* Python 3.9.1
* This was run on a 16-inch MacBook Pro running MacOS Catalina 10.15.7

Activate the Conda env specified by the `environment.yml` file. Activate it in the Julia package `Conda.jl` by doing the following:
```julia
julia> ENV["CONDA_JL_HOME"] = "/path/to/miniconda/envs/h3n2"  # change this to your path

pkg> build Conda
```

##  How to run:
You need to have a Unix system or a shell with common Unix commands like `cp` and `gunzip`. Using the environment above, and being in this directory, execute the `main.jl` file.
