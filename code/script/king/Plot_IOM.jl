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
    "--M_cluster"
    help = "Mass of the Plummer cluster (in solar masses)"
    arg_type = Float64
    default = 1.0e+5
    "--Rv_cluster"
    help = "Virial radius of the Plummer cluster (in kpc)"
    arg_type = Float64
    default = 0.011839088782478026
    "--run"
    help = "Run id"
    arg_type = Int64
    default = 63874531748065
    "--N"
    help = "Number for particles"
    arg_type = Int64
    default = 10^4

end
parsed_args = parse_args(tabargs)

const Mtot_Msun = parsed_args["M_cluster"]
const Rv_kpc = parsed_args["Rv_cluster"]
const run = parsed_args["run"]
const Npart = parsed_args["N"]

const path_to_script = @__DIR__
const path_data = path_to_script * "/../../../data/"

const srun = string(run)


# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_MSun_Myr = 4.49851e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_MSun_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 

const E_HU_in_MSun_kpc_Myr = M_HU_in_Msun^2 * G_in_kpc_MSun_Myr/R_HU_in_kpc
# MSun kpc/Myr

# data=readdlm("bary_snapshots_63874397777699.txt")
# plot(data[:,2], [data[:,3]], xlims=(-100,100), ylims=(-100,100), aspect_ratio=1, xticks=-100:20:100, yticks=-100:20:100, frame=:box, label=:false, xlabel="x [HU]", ylabel="y [HU]")
# data=readdlm("iom_snapshots_63874397777699.txt")
# plot(data[:,1], [data[:,4]], label=:false, xlabel="Time [HU]", ylabel="Energy [HU]", frame=:box)

function plot_data()

    data=readdlm(path_data*"iom_snapshots_"*srun*".txt")

    plt = plot(data[:,1] .* T_HU_in_Myr, [data[:,4] .* E_HU_in_MSun_kpc_Myr], 
        label=:false, 
        xlabel="Time [Myr]", 
        ylabel="Energy "*L"[M_{\odot}\,\mathrm{ kpc}^2 \, \mathrm{Myr}^{-2}]", 
        frame=:box)

    display(plt)
    readline()

    datab=readdlm(path_data*"bary_snapshots_"*srun*".txt")

    pltb = plot(datab[:,2] .* R_HU_in_kpc, [datab[:,3]] .* R_HU_in_kpc, 
            xlims=(-5,5), ylims=(-5,5), 
            aspect_ratio=1, 
            xticks=-5:1:5, yticks=-5:1:5, 
            frame=:box, label=:false, 
            xlabel="x [kpc]", ylabel="y [kpc]")

    display(pltb)
    readline()


    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/IOM_king_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    mkpath(path_data*"plot/")
    namefilep_pdf = path_data*"plot/bary_king_"*srun*".pdf"
    savefig(pltb, namefilep_pdf)


end

@time plot_data()