# Ptychohyla
Leapfrog integrator for (stellar) stream

## Aim of the code

This repository aims to evolve a set of self-gravitating Plummer globular clusters under the influence of a host potential, in order to study the creation of stellar streams.

## Generation of initial conditions

TODO: Characteristics of the clusters

The generation of initial conditions is done using the [PlummerPlus.py](https://github.com/michael-petersen/PlummerPlus) python package, in `code/IC_generator`. One can generate a set of initial conditions by modifying the bash script, `generate_IC_Plummer.sh`, which creates a text datafile in the folder `data_IC` (TODO).

## Integration

One may launch the script by using the command

```
$ julia -t 12 Main.jl
```

where `-t 12` launches a multi-threaded execution over 12 CPU cores. One may find the list of options in the `Args.jl` file. 

One is advised to change the default option of the `--path_dir` argument to the user's directory path, which one can obtain by typing `pwd` in the console after moving to the `Ptychohyla` folder.

As an exemple, the command

```
$ julia -t 12 Main.jl --dt 0.001 --N_dt 100 --eps 0.001
```

launches an execution with timestep 0.001 HU, snapshots saved every 100 timesteps and softening length of 0.001 HU. Snapshots are saved in the `data/snapshot` folder.