# Ptychohyla
Leapfrog integrator for (stellar) stream

## Aim of the code

This repository aims to evolve a set of self-gravitating Plummer globular clusters under the influence of a host potential, in order to study the creation of stellar streams.

## Generation of initial conditions

TODO: Characteristics of the clusters

The generation of initial conditions is done using the [PlummerPlus.py](https://github.com/michael-petersen/PlummerPlus) python package, in `code/IC_generator`. One can generate a set of initial conditions by modifying the bash script, `generate_IC_Plummer.sh`, which creates a text datafile in the folder `data_IC` (TODO).

## Integration

One may both generate the IC and perform a run by going to the `job` folder and launching the bash script

```
$ bash job_stream_run.sh
```

One may modify the run's parameters within the bash script file.
Snapshots are saved in the `data/snapshot` folder.