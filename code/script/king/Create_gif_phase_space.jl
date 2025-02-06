using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse


##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--M_cluster"
    help = "Mass of the Plummer cluster (in solar masses)"
    arg_type = Float64
    default = 1.0e+5
    "--Rv_cluster"
    help = "Virial radius of the Plummer cluster (in kpc)"
    arg_type = Float64
    default = 0.011839088782478026
    "--framerate"
    help = "Number of frames per second"
    arg_type = Int64
    default = 30
    "--run"
    help = "Run id"
    arg_type = Int64
    default = 63874531748065

end
parsed_args = parse_args(tabargs)

const Mtot_Msun = parsed_args["M_cluster"]
const Rv_kpc = parsed_args["Rv_cluster"]
const framepersec = parsed_args["framerate"]
const run = parsed_args["run"]

const path_to_script = @__DIR__
const path_data = path_to_script * "/../../../data/"


# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_MSun_Myr = 4.49851e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_MSun_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 
const V_HU_in_kpc_Myr = sqrt((G_in_kpc_MSun_Myr*M_HU_in_Msun)/R_HU_in_kpc)
const V_HU_in_km_s = V_HU_in_kpc_Myr * 977.792
const srun = string(run)

function plot_data()

    if (isfile(path_data*"snapshots_"*srun*"/.DS_Store"))
        rm(path_data*"snapshots_"*srun*"/.DS_Store")
    end

    listFile = readdir(path_data*"snapshots_"*srun*"/";join=true)
    nsnap = length(listFile)

    tab_time = zeros(Float64, nsnap)

    for i=1:nsnap
        interm = split(split(listFile[i],"_")[end],".")
        index = parse(Float64, interm[1]*"."*interm[2])
        tab_time[i] = index
    end

    p = sortperm(tab_time)

    # nsnap = 10
    # (x,z)
    anim = @animate for i=1:1:nsnap

        println("Progress = ", i/nsnap)
        namefile = listFile[p[i]]
        data = readdlm(namefile, header=false)
        interm = split(split(namefile,"_")[end],".")
        time = parse(Float64, interm[1]*"."*interm[2]) * T_HU_in_Myr
        time = round(time, digits=1)

        rmax = 5 # kpc
        zmax = 0.5
        vzmax = 10.0/1000
        s = 1.0

        # https://docs.juliaplots.org/latest/generated/attributes_plot/
        # https://stackoverflow.com/questions/71992758/size-and-colour-in-julia-scatter-plot

        scatter(data[:,3] .* R_HU_in_kpc, data[:, 6] .* V_HU_in_kpc_Myr , 
                xlabel=L"z"*" [kpc]", ylabel=L"v_z"*" [kpc/Myr]", 
                framestyle=:box, labels=:false,
                xlims=(-zmax, zmax), ylims=(-vzmax,vzmax), 
                # aspect_ratio=1, 
                size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                title="t = "*string(time)*" Myr")

    end

    mkpath(path_data*"gif/")
    namefile_gif = path_data*"gif/king_"*srun*"_z_vz.gif"
    gif(anim, namefile_gif, fps = framepersec)





end

@time plot_data()

