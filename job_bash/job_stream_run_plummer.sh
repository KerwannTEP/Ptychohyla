#! /bin/bash

# Plummer IC
N=10000
Q=0.0
SEED=0
MCLUSTER=100000 # Total mass in solar masses
RV=0.02 # Virial radius in kpc
DIST=4.0 # Distance of the cluster to the host potential's centre (in kpc)

# Dark halo
C=15.3 # Concentration parameter of the NFW halo profile
RS=16.0 # Scale radius of the NFW halo profile
MHALO=800000000000 # "Mass of the dark halo, M200, in solar masses

# Bulge
MBULGE=5000000000
ALPHA=-1.8
RC=1.9

# Disk
MDISK=68000000000
ADISK=3.0
BDISK=0.280

# Run parameters
TEND=1000.0 # Final time, in Henon units
DT=0.001 # Timestep, in Henon units
NDT=1000 # Frequency of snapshot save
EPS=0.001 # Softening length of the gravitational interaction
HAS_HOST=false # Is there a host potential ? (true or false)

# Perform the run

RUN=Main.jl
FOLDER_OUTPUT=/path/to/data/output/

RESTART=false # Is this run a restart ? (true or false)
ID=-1 # Set to -1 for automatic id creation. So to the run's id you want to restart if needed

cd ../code/src/plummer

julia -t 12 ${RUN} --Npart ${N} --q ${Q} --M_cluster ${MCLUSTER} --Rv_cluster ${RV} \
                --d_cluster ${DIST} --c ${C} --Rs_host ${RS} --Mvir_halo ${MHALO} \
                --M_bulge ${MBULGE} --alpha_bulge ${ALPHA} --rc_bulge ${RC} \
                --M_disk ${MDISK} --a_disk ${ADISK} --b_disk ${BDISK} \
                --t_end ${TEND} --dt ${DT} --N_dt ${NDT} --eps ${EPS} \
                --host ${HAS_HOST} \
                --folder_output ${FOLDER_OUTPUT} \
                --restart ${RESTART} --id ${ID}
                
