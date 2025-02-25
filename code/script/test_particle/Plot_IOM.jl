using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse
using NearestNeighbors
using LaTeXStrings

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--run"
    help = "Run id"
    arg_type = Int64
    default = 63876166088701
end
parsed_args = parse_args(tabargs)


const run = parsed_args["run"]


const path_to_script = @__DIR__
const path_data = path_to_script * "/../../../data/"

const srun = string(run)

function plot_data()

    data=readdlm(path_data*"snapshots_"*srun*"/iom_snapshots_"*srun*".txt")

    dataE = data[:,4]
    nsnap = length(dataE)


    datafE = abs.(dataE[2:nsnap] ./ dataE[1] .- 1.0)

    plt = plot(data[2:nsnap,1], [datafE], 
        label=:false, 
        xlabel="Time [HU]", 
        ylabel="Fractional energy [HU]", 
        frame=:box)

    display(plt)
    readline()

    datab=readdlm(path_data*"snapshots_"*srun*"/pos_snapshots_"*srun*".txt")

    pltb = plot(datab[:,2], [datab[:,3]], 
            # xlims=(-10,10), ylims=(-10,10), 
            aspect_ratio=1, 
            # xticks=-10:2:10, yticks=-10:2:10, 
            frame=:box, label=:false, 
            xlabel="x [HU]", ylabel="y [HU]")

    display(pltb)
    readline()


    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/IOM_test_particle_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    mkpath(path_data*"plot/")
    namefilep_pdf = path_data*"plot/pos_test_particle_"*srun*".pdf"
    savefig(pltb, namefilep_pdf)




end

@time plot_data()