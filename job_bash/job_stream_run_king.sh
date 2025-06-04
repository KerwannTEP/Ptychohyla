#! /bin/bash

# King IC
N=10000
MCLUSTER=100000 # Total mass in solar masses
RH=0.020 # Half-mass radius in kpc
DIST=10.0 # Distance of the cluster to the host potential's centre (in kpc)

# Host
HAS_HOST=true # Is there a host potential ? (true or false)
HOST_TYPE=MW2022

# Run parameters
TEND=200.0 # Final time, in Henon units
DT=0.01 # Timestep, in Henon units
NDT=100 # Frequency of snapshot save
EPS=0.01 # Softening length of the gravitational interaction, in Henon units

# Perform the run
RUN=Main.jl
FOLDER_OUTPUT=/path/to/data/output/

# Restart?
RESTART=false # Is this run a restart ? (true or false)
ID=-1 # Set to -1 for automatic id creation. So to the run's id you want to restart if needed

cd ../code/src/king

julia -t 12 ${RUN} --Npart ${N} --M_cluster ${MCLUSTER} --Rh_cluster ${RH} \
                --d_cluster ${DIST} \
                --t_end ${TEND} --dt ${DT} --N_dt ${NDT} --eps ${EPS} \
                --host ${HAS_HOST} --host_type ${HOST_TYPE} \
#                --folder_output ${FOLDER_OUTPUT} \
#                --restart ${RESTART} --id ${ID}
