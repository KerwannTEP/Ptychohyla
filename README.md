# Ptychohyla
Leapfrog integrator for (stellar) stream

## Aim of the code

This repository aims to evolve a set of self-gravitating Plummer globular clusters under the influence of a host potential, in order to study the creation of stellar streams.

## Installation

Install Julia by following the instructions at [https://julialang.org/downloads/platform/](https://julialang.org/downloads/platform/).

To invoke Julia in the Terminal, you need to make sure that the julia command-line program is in your `PATH`. 

On MacOS, we must create a link in `/usr/local/bin` (here for Julia 1.8):

```
$ sudo ln -s /Applications/Julia-1.8.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia
```

## Julia packages

Open the terminal in the folder `packages` and type

```
$ julia Install_pkg.jl
```

to install the necessary packages.

## Local Julia environment

One may create a local Julia environment by following the instructions at [https://jkrumbiegel.com/pages/2022-08-26-pkg-introduction/](https://jkrumbiegel.com/pages/2022-08-26-pkg-introduction/).

One can create a local environment called `LocalEnv` by following the next instructions. First, open `Julia` on the console.

```
$ module load julia/1.8.0
$ julia
```

Then, open the Pkg mode and create the local environment

```
julia> ]
(@v1.8) pkg> generate LocalEnv
```


One may import packages within that local environment by following the instructions at [https://pkgdocs.julialang.org/v1/environments/](https://pkgdocs.julialang.org/v1/environments/).

First, go to the folder containing `LocalEnv/` and open Julia. Then, add the packages (here, `SpecialFunctions`) in Pkg mode after activating the local environment.

```
julia> ]
(@v1.8) pkg> activate LocalEnv
(LocalEnv) pkg> add SpecialFunctions
```

One may also deactivate the local environment by writing in Pkg mode

```
(LocalEnv) pkg> activate
```

Alternatively, one may activate the local environment without having to switch to Pkg mode by writing 

```
julia> using Pkg; Pkg.activate("LocalEnv")
```

Similarly, deactivating the local environment can by done through the line

```
julia> using Pkg; Pkg.activate()
```



## Integration

One may perform a CPU run for a King sphere by going to the `job` folder and launching the bash script

```
$ bash job_stream_king_run.sh
```

One may modify the run's parameters within the bash script file.
Snapshots are saved in the `data/snapshot` folder.

## Running the GPU version

A GPU-accelerated version of `Ptychohyla`, using the `CUDA` library, can be found in the folder `code/src_gpu`. Should one have access to a computing cluster using SLURM, one may adapt the script given in the folder `job_slurm` to perform a run

```
$ sbatch job_stream_king.sh
```