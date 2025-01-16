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
    default = 1000
    "--c"
    help = "Concentration parameter of the NFW dark halo profile"
    arg_type = Float64
    default = 9.4
    "--q"
    help = "Anisotropy parameter"
    arg_type = Float64
    default = 0.0
    "--t_end"
    help = "Stopping time (in Henon units)"
    arg_type = Float64
    default = 120.0
    "--M_cluster"
    help = "Mass of the Plummer cluster (in solar masses)"
    arg_type = Float64
    default = 1.0e+5
    "--Rv_cluster"
    help = "Virial radius of the Plummer cluster (in kpc)"
    arg_type = Float64
    default = 2.00e-2
    "--Rs_host"
    help = "Scale radius of the NFW dark halo (in kpc)"
    arg_type = Float64
    default = 23.8
    "--d_cluster"
    help = "Distance of the cluster to the host potential's centre (in kpc)"
    arg_type = Float64
    default = 8.50e+0
    "--M_halo_200"
    help = "Mass of the dark halo, M200, within a sphere whose mean density is 200 times the critical density of the universe (in solar masses)"
    arg_type = Float64
    default = 0.97e+12
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
const c = parsed_args["c"]
const q = parsed_args["q"]
const time_end = parsed_args["t_end"]
const dt = parsed_args["dt"]
const N_dt = parsed_args["N_dt"]
const eps = parsed_args["eps"]

const Mtot_Msun = parsed_args["M_cluster"]
const Rv_kpc = parsed_args["Rv_cluster"]
const d_kpc =  parsed_args["d_cluster"]

const M_DH_200_Msun = parsed_args["M_halo_200"]
const Rs_kpc = parsed_args["Rs_host"]