using ArgParse

# https://github.com/carlrodriguez/CHC/blob/mw2014/src/potential/host_potential/static_analytic_potential/host_potential.h

const folder_dir = @__DIR__ 

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--M_cluster"
    help = "Mass of the cluster (in solar masses)"
    arg_type = Float64
    default = 1.0e+5
    "--Rv_cluster"
    help = "Virial radius of the cluster (in kpc)"
    arg_type = Float64
    default = 2.00e-2
    "--d_cluster"
    help = "Distance of the cluster to the host potential's centre (in kpc)"
    arg_type = Float64
    default = 4.0e+0

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

    "--host"
    help = "Do we include a host potential (true or false) ?"
    arg_type = Bool
    default = true
    "--host_type"
    help = "Type of host potential"
    arg_type = String
    default = "MW2014"

    "--folder_output"
    help = "Output folder of the data"
    arg_type = String
    default = folder_dir  * "/../../../data/"

end
parsed_args = parse_args(tabargs)


const Mtot_Msun = parsed_args["M_cluster"]
const Rv_kpc = parsed_args["Rv_cluster"]
const d_kpc =  parsed_args["d_cluster"]

const time_end = parsed_args["t_end"]
const dt = parsed_args["dt"]
const N_dt = parsed_args["N_dt"]
const eps = parsed_args["eps"]

const HAS_HOST = parsed_args["host"]
const HOST_TYPE = parsed_args["host_type"]

const folder_output = parsed_args["folder_output"]
