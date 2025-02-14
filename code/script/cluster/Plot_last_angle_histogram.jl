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
    default = 0.01228105689696044
    "--run"
    help = "Run id"
    arg_type = Int64
    default = 63875213898039
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

# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_MSun_Myr = 4.49851e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_MSun_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 

const srun = string(run)

const nb_neigh = 6

function get_data()

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

    tab_pos = zeros(Float64, 3, Npart)

    nb_neigh = 6


    tab_pos[1, :] = data[:, 1] .* R_HU_in_kpc * 1.0
    tab_pos[2, :] = data[:, 2] .* R_HU_in_kpc * 1.0
    tab_pos[3, :] = data[:, 3] .* R_HU_in_kpc * 1.0

    return tab_pos, time
end



function compute_tab_delta_phi()

    tab_pos, time = get_data()
    
    # Compute density centre
    # Recentre
    tab_dens = zeros(Float64, Npart)


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

    println("(Xc,Yc) = ",(Xc, Yc))

    phi_avg = atan(Yc, Xc)

    tab_delta_phi = zeros(Float64, Npart)

    Threads.@threads for i=1:Npart
        x = tab_pos[1,i]
        y = tab_pos[2,i]

        phi = atan(y, x)
        tab_delta_phi[i] = ((phi - phi_avg) * 180.0/pi)

        if (180 <= tab_delta_phi[i] <= 360.0)
            tab_delta_phi[i] = -(360.0 - tab_delta_phi[i])
        end
        # tab_delta_phi[i] = (phi ) * 180.0/pi
    end

    return tab_delta_phi, time

end

function compute_histo(phi_max::Float64, dphi::Float64=1.0)

    tab_delta_phi, time = compute_tab_delta_phi()
    nb_phi_pos = floor(Int64, phi_max/dphi)
    nb_phi = 2 * nb_phi_pos

    tab_dphi = [-phi_max + dphi * (i-0.5) for i=1:nb_phi]
    tab_count = zeros(Float64, nb_phi)

    # phi_i <= phi < phi_{i+1}

    for k=1:Npart 
        phi = tab_delta_phi[k]

        # -phi_max + dphi * (i-0.5) <= phi < -phi_max + dphi * (i+0.5)
        # dphi * (i-0.5) <= phi + phi_max < dphi * (i+0.5)
        # i <= (phi + phi_max)/dphi + 0.5 < i+1

        i = floor(Int64, (phi + phi_max)/ dphi + 0.5)

        if (1 <= i <= nb_phi) 
            tab_count[i] += 1.0
        end
    end

    for i=1:nb_phi

        if(tab_count[i] == 0.0)

            tab_count[i] = 0.001

        end

    end


    return tab_dphi, tab_count, time

end

function plot_delta_phi(delta_phi_max::Float64, dphi::Float64=1.0)

    tab_dphi, tab_count, time = compute_histo(delta_phi_max, dphi)
    time = round(time, digits=1)

    p = plot(tab_dphi, [tab_count],
            yaxis=:log10, 
            xlims=(-delta_phi_max, delta_phi_max),
            ylims=(0.1, 10^4),
            frame=:box, label=:false, 
            linetype=:steppost,
            xticks=-180:20:180,
            xlabel=L"\Delta \phi"*" [deg]", ylabel="Count",
            title="t = "*string(time)*" Myr")
    
    # p = histogram(tab_delta_phi, 
    #             yaxis=:log10, 
    #             xlims=(-delta_phi_max, delta_phi_max),
    #             ylims=(0.1, 10^3),
    #             frame=:box, label=:false, 
    #             # normalize=:pdf,
    #             xlabel=L"\Delta \phi"*" [deg]", ylabel="Count",
    #             title="t = "*string(time)*" Myr")

    display(p)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/angles_cluster_"*srun*".pdf"
    savefig(p, namefile_pdf)

end

plot_delta_phi(180.0, 1.0)