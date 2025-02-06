using ArgParse

# https://github.com/carlrodriguez/CHC/blob/mw2014/src/potential/host_potential/static_analytic_potential/host_potential.h

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--Npart"
    help = "Number of stars in the cluster"
    arg_type = Int64
    default = 10000
    "--M_cluster"
    help = "Mass of the King cluster (in solar masses)"
    arg_type = Float64
    default = 1.0e+5
    "--Rh_cluster"
    help = "Half-mass radius of the King cluster (in kpc)"
    arg_type = Float64
    default = 5.00e-2
    "--d_cluster"
    help = "Distance of the cluster to the host potential's centre (in kpc)"
    arg_type = Float64
    default = 2.0e+0



    "--c"
    help = "Halo (NFW) concentration"
    arg_type = Float64
    default = 15.3
    "--Rs_host"
    help = "Halo (NFW) scale radius (in kpc)"
    arg_type = Float64
    default = 16.0
    "--Mvir_halo"
    help = "Virial massass of the dark halo (NFW), (in solar masses)"
    arg_type = Float64
    default = 0.8e+12

    "--M_bulge"
    help = "Mass of the bulge (in solar masses)"
    arg_type = Float64
    default = 0.5e+10
    "--alpha_bulge"
    help = "Power law index of the bulge"
    arg_type = Float64
    default = -1.8
    "--rc_bulge"
    help = "Exponential cutoff radius of the bulge (in kpc)"
    arg_type = Float64
    default = 1.9e+0

    "--M_disk"
    help = "Mass of the Miyamoto-Nagai disk (in solar masses)"
    arg_type = Float64
    default = 6.8e+10
    "--a_disk"
    help = "Scale length of the Miyamoto-Nagai disk (in kpc)"
    arg_type = Float64
    default = 3.0e+0
    "--b_disk"
    help = "Scale height of the Miyamoto-Nagai disk (in kpc)"
    arg_type = Float64
    default = 0.280e+0


    "--t_end"
    help = "Stopping time (in Henon units)"
    arg_type = Float64
    default = 600.0
    "--dt"
    help = "Timestep"
    arg_type = Float64
    default = 1.0e-3
    "--N_dt"
    help = "Save snapshots every N_dt timesteps"
    arg_type = Int64
    default = 100
    "--eps"
    help = "Softening length"
    arg_type = Float64
    default = 0.001
end
parsed_args = parse_args(tabargs)


const Npart = parsed_args["Npart"]
const Mtot_Msun = parsed_args["M_cluster"]
const Rh_kpc = parsed_args["Rh_cluster"]
const d_kpc =  parsed_args["d_cluster"]


const Mvir_Msun = parsed_args["Mvir_halo"]
const Rs_kpc = parsed_args["Rs_host"]
const c = parsed_args["c"]

const M_bulge_Msun = parsed_args["M_bulge"]
const alpha_bulge = parsed_args["alpha_bulge"]
const rc_bulge_kpc = parsed_args["rc_bulge"]

const M_disk_Msun = parsed_args["M_disk"]
const a_disk_kpc = parsed_args["a_disk"]
const b_disk_kpc = parsed_args["b_disk"]


const time_end = parsed_args["t_end"]
const dt = parsed_args["dt"]
const N_dt = parsed_args["N_dt"]
const eps = parsed_args["eps"]