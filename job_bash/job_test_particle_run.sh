#! /bin/bash

# Test star IC
RV=0.02 # Virial radius in kpc
DIST=4.0 # Distance of the cluster to the host potential's centre (in kpc)

# Host
HAS_HOST=true # Is there a host potential ? (true or false)
HOST_TYPE=MW2022

# Run parameters
TEND=500.0 # Final time, in Henon units
DT=0.01 # Timestep, in Henon units
NDT=100 # Frequency of snapshot save
EPS=0.01 # Softening length of the gravitational interaction

# Perform the run
RUN=Main.jl
FOLDER_OUTPUT=/path/to/data/output/

cd ../code/src/test_star

julia -t 12 ${RUN} --d_cluster ${DIST} --Rv_cluster ${RV} \
                --t_end ${TEND} --dt ${DT} --N_dt ${NDT} --eps ${EPS} \
                --host ${HAS_HOST} --host_type ${HOST_TYPE} \
                # --folder_output ${FOLDER_OUTPUT}
                
