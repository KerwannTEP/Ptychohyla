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
    default = 2.00e-2
    "--framerate"
    help = "Number of frames per second"
    arg_type = Int64
    default = 60
    "--path_dir"
    help = "IC data file in the cluster's frame"
    arg_type = String
    default = "/Users/tep/Documents/Astro_projects/UNC/Julia_fewbody_integrator/"
end
parsed_args = parse_args(tabargs)

const Mtot_Msun = parsed_args["M_cluster"]
const Rv_kpc = parsed_args["Rv_cluster"]
const framepersec = parsed_args["framerate"]
const path_dir = parsed_args["path_dir"]

# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_Mpc_Myr = 4.49e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_Mpc_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 



function plot_data()

    if (isfile(path_dir*"data/snapshots/.DS_Store"))
        rm(path_dir*"data/snapshots/.DS_Store")
    end

    listFile = readdir(path_dir*"data/snapshots/";join=true)
    nsnap = length(listFile)

    tab_time = zeros(Float64, nsnap)

    for i=1:nsnap
        interm = split(split(listFile[i],"_")[end],".")
        index = parse(Float64, interm[1]*"."*interm[2])
        tab_time[i] = index
    end

    p = sortperm(tab_time)

   
    namefile = listFile[p[nsnap]]
    data = readdlm(namefile, header=false)
    interm = split(split(namefile,"_")[end],".")
    time = parse(Float64, interm[1]*"."*interm[2]) * T_HU_in_Myr
    time = round(time, digits=1)

    rmax = 500
    p = scatter(data[:,1], data[:, 2], xlabel=L"x"*" [HU]", ylabel=L"y"*" [HU]", framestyle=:box, labels=:false, xlims=(-rmax, rmax), ylims=(-rmax,rmax), aspect_ratio=1, size=(800,800), left_margin = [2mm 0mm], right_margin = [2mm 0mm], title="t = "*string(time)*" Myr")

    mkpath(path_dir*"data/plot/")
    namefile_pdf = path_dir*"data/plot/plummer_last_snapshot.pdf"
    savefig(p, namefile_pdf)


end

@time plot_data()

