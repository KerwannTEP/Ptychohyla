#!/bin/bash
#SBATCH --job-name=stream
#SBATCH --partition=pscomp
#SBATCH --gres=gpu
#SBATCH --exclusive
#SBATCH --time=2-00:00:00
#SBATCH --output=./stream_king.out
#SBATCH --mail-type=TIME_LIMIT,END,FAIL
#SBATCH --mail-user=<your email>

# Use whichever versions are compatible with your GPU driver
# Incompatible and not up-to-date versions may cause GPU performance to randomly drop to 0%, and effectively bringing the code to a stop
module purge
module load julia/1.11.3
module load cuda/12.8 
module load hdf5/1.14.6-mpi-intel
module load inteloneapi/2024.0

NBTHREADSGPU=128

# King IC
N=10000
MCLUSTER=100000 # Total mass in solar masses
RH=0.020 # Half-mass radius in kpc
DIST=20.0 # Distance of the cluster to the host potential's centre (in kpc)

# Host
HAS_HOST=true # Is there a host potential ? (true or false)
HOST_TYPE=MW2022

# Use custom IC file ?
CUSTOM_IC=false

# Run parameters
TEND=2000.0 # Final time, in Henon units
DT=0.0001 # Timestep, in Henon units
NDT=10000 # Frequency of snapshot save
EPS=0.001 # Softening length of the gravitational interaction, in Henon units

# Perform the run
RUN=Main.jl
FOLDER_OUTPUT=/path/to/data/output/

RESTART=false # Is this run a restart ? (true or false)
ID=-1 # Set to -1 for automatic id creation. So to the run's id you want to restart if needed

cd ../code/src_gpu/king

julia -t auto ${RUN} --Npart ${N} --M_cluster ${MCLUSTER} --Rh_cluster ${RH} \
                --d_cluster ${DIST} --custom_IC ${CUSTOM_IC} \
                --t_end ${TEND} --dt ${DT} --N_dt ${NDT} --eps ${EPS} \
                --host ${HAS_HOST} --host_type ${HOST_TYPE} \
                --nbThreadsPerBlocks ${NBTHREADSGPU} \
                --folder_output ${FOLDER_OUTPUT} \
                --restart ${RESTART} --id ${ID}