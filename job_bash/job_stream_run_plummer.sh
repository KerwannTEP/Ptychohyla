#! /bin/bash

# Plummer IC
N=10000
Q=0.0
SEED=0
MCLUSTER=100000 # Total mass in solar masses
RV=0.02 # Virial radius in kpc
DIST=4.0 # Distance of the cluster to the host potential's centre (in kpc)

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

RESTART=false # Is this run a restart ? (true or false)
ID=-1 # Set to -1 for automatic id creation. So to the run's id you want to restart if needed

cd ../code/src/plummer

julia -t 12 ${RUN} --Npart ${N} --q ${Q} --M_cluster ${MCLUSTER} --Rv_cluster ${RV} \
                --d_cluster ${DIST} \
                --t_end ${TEND} --dt ${DT} --N_dt ${NDT} --eps ${EPS} \
                --host ${HAS_HOST} --host_type ${HOST_TYPE} \
                # --folder_output ${FOLDER_OUTPUT} \
                # --restart ${RESTART} --id ${ID}
                
