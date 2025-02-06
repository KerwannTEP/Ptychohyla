using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse
using NearestNeighbors

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
    default = 0.06478848136966434
    "--framerate"
    help = "Number of frames per second"
    arg_type = Int64
    default = 30
    "--run"
    help = "Run id"
    arg_type = Int64
    default = 63874516769458
    "--N"
    help = "Number for particles"
    arg_type = Int64
    default = 10^4

end
parsed_args = parse_args(tabargs)

const Mtot_Msun = parsed_args["M_cluster"]
const Rv_kpc = parsed_args["Rv_cluster"]
const framepersec = parsed_args["framerate"]
const run = parsed_args["run"]
const Npart = parsed_args["N"]

const path_to_script = @__DIR__
const path_data = path_to_script * "/../../../data/"

const srun = string(run)

# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_MSun_Myr = 4.49e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_MSun_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 


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

    # Center 

    tab_pos = zeros(Float64, 3, Npart)
    tab_dens = zeros(Float64, Npart)

    nb_neigh = 6

    # nsnap = 1000

    anim = @animate for i=1:1:nsnap

        println("Progress = ", i/nsnap)
        namefile = listFile[p[i]]
        data = readdlm(namefile, header=false)
        interm = split(split(namefile,"_")[end],".")
        time = parse(Float64, interm[1]*"."*interm[2]) * T_HU_in_Myr
        time = round(time, digits=1)

        tab_pos[1, :] = data[:, 1] .* R_HU_in_kpc * 1.0 # In parsecs
        tab_pos[2, :] = data[:, 2] .* R_HU_in_kpc * 1.0 # In parsecs
        tab_pos[3, :] = data[:, 3] .* R_HU_in_kpc * 1.0
        # Compute density centre
        # Recentre

        tree_neigh = KDTree(tab_pos)

        Threads.@threads for i=1:Npart # Loop over all the particles in the cluster. ATTENTION, we go over all particles, not just the one within the core
            #####
            tab_indneigh, tab_distneigh = knn(tree_neigh, tab_pos[:,i], nb_neigh+1, true) # Finding the (nb_neigh+1) neighbors. ATTENTION, to the `+1', because we are accounting for the distance between a particle and itself
            #####
            # Getting the index and distance of the nb_neigh neighbor
            # ATTENTION, to the `+1', because we are accounting for the distance between a particle and itself
            ind_neigh = tab_indneigh[nb_neigh+1] # Index of the nb_neigh neighbor
            dist_neigh = tab_distneigh[nb_neigh+1] # Distance between the current particle and the nb_neigh neighbor
            #####
            # We can finally have an estimate of the local density,
            # see Eq. (II.2) of (Casertano & Hut, Ap.J. 298, 80)
            # ATTENTION, we set the individual mass to 1.0, as for single-mass system, it cancels out by proportionality
            vol_sphere = 4.0*pi/(3.0)*dist_neigh^(3) # Volume of the sphere between particle i and its nb_neigh neighbor
            dens_loc = (nb_neigh - 1.0)/(vol_sphere) # Local estimation of the density
            #####
            tab_dens[i] = dens_loc # Filling in the array of densities
        end

        Xc = 0.0
        Yc = 0.0
        rho_tot = 0.0

        for i=1:Npart 
            rho_i = tab_dens[i] # Estimation of the local density
            x_i, y_i, z_i = tab_pos[1,i], tab_pos[2,i], tab_pos[3,i] # Position of the current particle
    
            rho_tot  += rho_i # Updating the total density
            Xc += x_i*rho_i # Updating the centre's position
            Yc += y_i*rho_i # Updating the centre's position
        end

        Xc /= rho_tot # Rescaling the centre's position
        Yc /= rho_tot # Rescaling the centre's position
        
        datax = tab_pos[1, :] .- Xc
        datay = tab_pos[2, :] .- Yc

        theta = atan(Yc, Xc)

        dataxc = cos(-theta) .* datax - sin(-theta) .* datay
        datayc = sin(-theta) .* datax + cos(-theta) .* datay


        rmax = 0.5 # kpc
        s = 1.0

        # https://docs.juliaplots.org/latest/generated/attributes_plot/
        # https://stackoverflow.com/questions/71992758/size-and-colour-in-julia-scatter-plot

        plot([0, -rmax/2.0],[0.0,0.0 ],
                xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]", 
                framestyle=:box, labels=:false,
                xlims=(-2*rmax, 2*rmax), ylims=(-3*rmax,3*rmax), 
                aspect_ratio=1, size=(400,600), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                #markersize=s, color=:white, 
                title="t = "*string(time)*" Myr",
                arrow=true,color=:red,linewidth=2,label=:false)

        scatter!(dataxc, datayc, 
                #xlabel=L"x"*" [pc]", ylabel=L"y"*" [pc]", 
                #framestyle=:box, 
                labels=:false,
                xlims=(-rmax, rmax), ylims=(-rmax,rmax), 
                aspect_ratio=1, size=(800,800), 
                #left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                #background_color = :black,
                markersize=s, color=:white)#, 
                #title="t = "*string(time)*" Myr")

    end

    mkpath(path_data*"gif/")
    namefile_gif = path_data*"gif/king_centered_"*srun*".gif"
    gif(anim, namefile_gif, fps = framepersec)


end

@time plot_data()

