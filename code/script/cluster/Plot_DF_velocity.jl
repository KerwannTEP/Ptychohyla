using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse

# TODO


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
    default = 0.01228105689696044
    "--run"
    help = "Run id"
    arg_type = Int64
    default = 63875434463958

end
parsed_args = parse_args(tabargs)

const Mtot_Msun = parsed_args["M_cluster"]
const Rv_kpc = parsed_args["Rv_cluster"]
const run = parsed_args["run"]

const path_to_script = @__DIR__
const path_data = path_to_script * "/../../../data/"

# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_MSun_Myr = 4.49851e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_MSun_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 

const srun = string(run)


function read_data()

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

   
    namefile = listFile[p[nsnap]]
    data = readdlm(namefile, header=false)
    interm = split(split(namefile,"_")[end],".")
    time = parse(Float64, interm[1]*"."*interm[2]) * T_HU_in_Myr
    time = round(time, digits=1)

    rmax = 5.0 # kpc
    s = 1.0

    # Plot circular orbit 
    d_c = 4.0 # kpc

    nba = 200
    tab_ang = range(0, 2*pi, length=200)
    tab_orbit = zeros(Float64, nba, 2)

    for i=1:nba 
        phi = tab_ang[i]
        tab_orbit[i, 1] = d_c * cos(phi)
        tab_orbit[i, 2] = d_c * sin(phi)
    end

    p = plot(tab_orbit[:,1] , tab_orbit[:, 2], 
                xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]", 
                framestyle=:box, labels="Circular orbit", 
                legend=:topright,
                xlims=(-rmax, rmax), ylims=(-rmax,rmax), 
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                linestyle=:dash,
                markersize=s, color=:red)

    # https://docs.juliaplots.org/latest/generated/attributes_plot/
    # https://stackoverflow.com/questions/71992758/size-and-colour-in-julia-scatter-plot

    scatter!(p, data[:,1] .* R_HU_in_kpc, data[:, 2] .* R_HU_in_kpc , 
                xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]", 
                framestyle=:box, labels=:false, 
                xlims=(-rmax, rmax), ylims=(-rmax,rmax), 
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                title="t = "*string(time)*" Myr")

    display(p)
    readline()


    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/cluster_"*srun*"_last_snapshot.pdf"
    savefig(p, namefile_pdf)


end

@time plot_data()

