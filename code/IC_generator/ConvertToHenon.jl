using ArgParse
using DelimitedFiles

tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--N"
    help = "Number of particles"
    arg_type = Int64
    default = 10
   "--q"     
    help = "Anisotropy parameter"              
    arg_type = Float64 
    default = 0.0
    "--seed"     
    help = "Random seed"              
    arg_type = Int64 
    default = 0
end

parsed_args = parse_args(tabargs)

const N = parsed_args["N"]
const q = parsed_args["q"]
const rs = parsed_args["seed"]

# https://www.geeksforgeeks.org/storing-output-on-a-file-in-julia/



data = readdlm("output.txt")
posvel = data[:,2:7]

const G = 1.0
const M = 1.0
const b = 3*pi/16
const vh = sqrt(G*M/b)

# Convert to Henon units
posvel[:,1:3] = posvel[:,1:3] .* b
posvel[:,4:6] = posvel[:,4:6] .* vh   

mkpath("data_IC")
writedlm("data_IC/plummer_sphere_n_"*string(N)*"_q_"*string(q)*".txt", posvel)

