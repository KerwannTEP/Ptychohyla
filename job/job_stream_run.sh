#! /bin/bash

# Plummer IC
N=1000
Q=0.0
SEED=0
MCLUSTER=100000 # Total mass in solar masses
RV=0.02 # Virial radius in kpc
DIST=8.5 # Distance of the cluster to the host potential's centre (in kpc)

# Host IC
C=9.4 # Concentration parameter of the NFW halo profile
RS=23.8 # Scale radius of the NFW halo profile
MHALO=970000000000 # "Mass of the dark halo, M200, in solar masses

# Run parameters
TEND=120.0 # Final time, in Henon units
DT=0.001 # Timestep, in Henon units
NDT=100 # Frequency of snapshot save
EPS=0.001 # Softening length of the gravitational interaction

# Generate initial Plummer cluster

cd ../code/IC_generator

source PlummerPlus/venv/bin/activate
./PlummerPlus/PlummerPlus.py -n ${N} -q ${Q} -rs ${SEED}
julia ConvertToHenon.jl --N ${N} --q ${Q} --seed ${SEED}
rm output.txt
deactivate

# Perform the run

cd ../src

RUN=Main.jl

julia -t 8 ${RUN} --Npart ${N} --c ${C} --q ${Q} --t_end ${TEND} --M_cluster ${MCLUSTER} --Rv_cluster ${RV} --Rs_host ${RS} --d_cluster ${DIST} --M_halo_200 ${MHALO} --dt ${DT} --N_dt ${NDT} --eps ${EPS}


