
using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse
using NearestNeighbors
using LaTeXStrings
using HDF5

const path_to_script = @__DIR__

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
    default = 63875434463958 #63875411207673
    "--N"
    help = "Number for particles"
    arg_type = Int64
    default = 10^4
    "--data_output"
    help = "Output folder of the data"
    arg_type = String
    default = path_to_script * "/../../../data/"

end
parsed_args = parse_args(tabargs)

const Mtot_Msun = parsed_args["M_cluster"]
const Rv_kpc = parsed_args["Rv_cluster"]
const run = parsed_args["run"]
const Npart = parsed_args["N"]
const path_data = parsed_args["data_output"]

const srun = string(run)

# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_MSun_Myr = 4.49851e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_MSun_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 

const mass = 1.0/Npart
const nb_neigh = 10

function get_data()

    if (isfile(path_data*"snapshots_"*srun*"/.DS_Store"))
        rm(path_data*"snapshots_"*srun*"/.DS_Store")
    end

    listFiles = readdir(path_data*"snapshots_"*srun*"/";join=true)
    nsnap = length(listFiles)
    tabt = zeros(Float64, nsnap)

    for i=1:nsnap 
        interm = split(split(listFiles[i],"_")[end],".")
        interm = interm[1]*"."*interm[2]
        time = parse(Float64, interm)

        tabt[i] = time 
    end

    p = sortperm(tabt)

    sortedFiles = Array{String}(undef, nsnap)
    tabtsort = zeros(Float64, nsnap)

    for i=1:nsnap 

        tabtsort[i] = tabt[p[i]]
        sortedFiles[i] = listFiles[p[i]]

    end

    tab_IOM = zeros(Float64, nsnap, 6) # time, E_tot, Lx, Ly, Lz, nb_unbound
    tab_lag = zeros(Float64, nsnap, 5) # 1% 10% 20% 50% 90% 
    tab_nc = zeros(Float64, nsnap) # Central density

    for isnap=1:nsnap

        println("Progress : ", round(isnap/nsnap, digits=4))

        namefile = sortedFiles[isnap]
        time = tabtsort[isnap]

        data_stars = readdlm(namefile) # x, y, z, vx, vy, vz, Uint, Uc


        # tab_vb_t = zeros(Float64, Threads.nthreads(), 3)
        K_t =  zeros(Float64, Threads.nthreads())
        Uc_t =  zeros(Float64, Threads.nthreads())
        Uh_t =  zeros(Float64, Threads.nthreads())
        L_t =  zeros(Float64, Threads.nthreads(), 3)
        

        Threads.@threads for i=1:Npart 

            tid = Threads.threadid()
            x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]
            

            v2 = vx^2 + vy^2 + vz^2

            # # Barycenter 
            # tab_vb_t[tid, 1] += vx
            # tab_vb_t[tid, 2] += vy
            # tab_vb_t[tid, 3] += vz

            # Kinetic energy 
            K_t[tid] += 0.5 * mass * v2 

            # Host potential energy
            Uh_t[tid] += Uc

            # Cluster potential energy 
            Uc_t[tid] += 0.5*Uint

                
            # Angular momenta 
            Lx = y*vz - z*vy 
            Ly = z*vx - x*vz 
            Lz = x*vy - y*vx 

            L_t[tid, 1] += mass * Lx
            L_t[tid, 2] += mass * Ly
            L_t[tid, 3] += mass * Lz

        end
       
        # Vbx = 0.0
        # Vby = 0.0
        # Vbz = 0.0

        K = 0.0
        U = 0.0
        L = zeros(Float64, 6)


        for tid=1:Threads.nthreads()
            # Vbx += tab_vb_t[tid, 1]/Npart
            # Vby += tab_vb_t[tid, 2]/Npart
            # Vbz += tab_vb_t[tid, 3]/Npart

            K += K_t[tid]
            U += Uh_t[tid] + Uc_t[tid]
            L[1] += L_t[tid, 1]
            L[2] += L_t[tid, 2]
            L[3] += L_t[tid, 3]
    
        end

        Etot = K + U

        # Compute number of unbound stars
        n_unbound_t = zeros(Int64, Threads.nthreads())

        tab_pos = zeros(Float64, 3, Npart)

        tab_pos[1, :] = data_stars[:, 1]
        tab_pos[2, :] = data_stars[:, 2]
        tab_pos[3, :] = data_stars[:, 3]

        tab_dens = zeros(Float64, Npart)

        tree_neigh = KDTree(tab_pos)

        Threads.@threads for k=1:Npart # Loop over all the particles in the cluster. ATTENTION, we go over all particles, not just the one within the core
            #####
            tab_indneigh, tab_distneigh = knn(tree_neigh, tab_pos[:,k], nb_neigh+1, true) # Finding the (nb_neigh+1) neighbors. ATTENTION, to the `+1', because we are accounting for the distance between a particle and itself
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
            tab_dens[k] = dens_loc # Filling in the array of densities
        end

        tab_dens_pos_t = zeros(Float64, Threads.nthreads(), 3)
        tab_dens_vel_t = zeros(Float64, Threads.nthreads(), 3)
        tab_dens_t =  zeros(Float64, Threads.nthreads())
        

        Threads.@threads for i=1:Npart 

            tid = Threads.threadid()
            x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]
            rho = tab_dens[i]
            
            tab_dens_pos_t[tid, 1] += x * rho
            tab_dens_pos_t[tid, 2] += y * rho
            tab_dens_pos_t[tid, 3] += z * rho

            tab_dens_vel_t[tid, 1] += vx * rho
            tab_dens_vel_t[tid, 2] += vy * rho
            tab_dens_vel_t[tid, 3] += vz * rho

            tab_dens_t[tid] += rho

            
        end

        rho_tot = 0.0
        Xc = 0.0
        Yc = 0.0
        Zc = 0.0
        Vcx = 0.0
        Vcy = 0.0
        Vcz = 0.0

        for tid=1:Threads.nthreads()

            Xc += tab_dens_pos_t[tid, 1]
            Yc += tab_dens_pos_t[tid, 2]
            Zc += tab_dens_pos_t[tid, 3]

            Vcx += tab_dens_vel_t[tid, 1]
            Vcy += tab_dens_vel_t[tid, 2]
            Vcz += tab_dens_vel_t[tid, 3]

            rho_tot += tab_dens_t[tid]


        end

        Vcx /= rho_tot
        Vcy /= rho_tot
        Vcz /= rho_tot

        Xc /= rho_tot
        Yc /= rho_tot
        Zc /= rho_tot

        

        Threads.@threads for i=1:Npart

            tid = Threads.threadid()
    
            # Unbound particles 
            x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]

       
            # Let x_c is the density center of the cluster.
            # It is a proxy for the cluster's center 
            # What about the velocity vc of that center ?
            # Calculation show that vc = dxc/dt + fluctuations 1/N

            vc_x = vx - Vcx
            vc_y = vy - Vcy
            vc_z = vz - Vcz
    
            Ec = 0.5 * mass * (vc_x^2 + vc_y^2 + vc_z^2) + Uint
    
            if (Ec >= 0.0)
                n_unbound_t[tid] += 1
            end
    
        end

        n_unbound = 0

        for tid=1:Threads.nthreads()
            n_unbound += n_unbound_t[tid]
        end
    
        tab_IOM[isnap,1] = time
        tab_IOM[isnap,2] = Etot
        tab_IOM[isnap,3] = L[1]
        tab_IOM[isnap,4] = L[2]
        tab_IOM[isnap,5] = L[3]
        tab_IOM[isnap,6] = n_unbound

        # Lagrange radii 
        n_bound = Npart -n_unbound
        tabr = zeros(Float64, n_bound)

        index = 1
        for i=1:Npart

            
    
            # Unbound particles 
            x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]

            

            # Let x_c is the density center of the cluster.
            # It is a proxy for the cluster's center 
            # What about the velocity vc of that center ?
            # Calculation show that vc = dxc/dt + fluctuations 1/N

            vc_x = vx - Vcx
            vc_y = vy - Vcy
            vc_z = vz - Vcz
    
            Ec = 0.5 * mass * (vc_x^2 + vc_y^2 + vc_z^2) + Uint
    
            if (Ec < 0.0)
                r = sqrt((x-Xc)^2 + (y-Yc)^2 + (z-Zc)^2) 
                tabr[index] = r
                index += 1
            end
    
        end

        tabr = sort(tabr)

        index_1 = floor(Int64, n_bound * 0.01)
        index_10 = floor(Int64, n_bound * 0.10)
        index_20 = floor(Int64, n_bound * 0.20)
        index_50 = floor(Int64, n_bound * 0.50)
        index_90 = floor(Int64, n_bound * 0.90)

        tab_lag[isnap, 1] = tabr[index_1]
        tab_lag[isnap, 2] = tabr[index_10]
        tab_lag[isnap, 3] = tabr[index_20]
        tab_lag[isnap, 4] = tabr[index_50]
        tab_lag[isnap, 5] = tabr[index_90]

        index_05 = floor(Int64, n_bound * 0.005)
        tab_nc[isnap] = index_05/(4.0/3.0*pi*tabr[index_05]^3) # Number/Volume
        

    end


    return tab_IOM, tab_lag, tab_nc

end

function plot_data!()

    tab_IOM, tab_lag, tab_nc = get_data()
    datat = tab_IOM[:, 1]  .* T_HU_in_Myr
    dataE = tab_IOM[:, 2]
    dataLx = tab_IOM[:, 3]
    dataLy = tab_IOM[:, 4]
    dataLz = tab_IOM[:, 5]
    dataUnbound = tab_IOM[:, 6]

    namefile_hf5 = path_data*"iom_cluster_"*srun*".hf5"
    file = h5open(namefile_hf5, "w")

    write(file, "data_time", datat)
    write(file, "data_Etot", dataE)
    write(file, "data_Lx", dataLx)
    write(file, "data_Ly", dataLy)
    write(file, "data_Lz", dataLz)
    write(file, "data_unbound_frac", dataUnbound ./ Npart)

    write(file, "Npart", Npart)
    write(file, "kpc_per_HU", R_HU_in_kpc)
    write(file, "G_in_kpc_MSun_Myr", G_in_kpc_MSun_Myr)
    write(file, "Msun_per_HU", M_HU_in_Msun)

    write(file, "Myr_per_HU", T_HU_in_Myr)

    close(file)


    # Energy
    n = length(dataE)
    dataFracE = abs.(1.0 .- dataE[2:n] ./ dataE[1])

    plt = plot(datat[2:n], [dataFracE], 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel="Fractional energy", 
        yaxis=:log10,
        yticks=10.0 .^ (-15:1:0),
        xlims=(0, datat[n]),
        frame=:box)

    display(plt)
    readline()
    

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/E_frac_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    # Angular momentum

    # x y
    plt = plot(datat, [dataLx dataLy], 
        labels=[L"L_x" L"L_y"], 
        xlabel="Time [Myr]", 
        ylabel="Angular momentum", 
        # yaxis=:log10,
        # yticks=10.0 .^ (-15:1:0),
        xlims=(0, datat[n]),
        frame=:box)

    display(plt)
    readline()
 
    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/Lxy_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    # z 
    plt = plot(datat, [dataLz], 
        label=L"L_z", 
        xlabel="Time [Myr]", 
        ylabel="Angular momentum", 
        # yaxis=:log10,
        # yticks=10.0 .^ (-15:1:0),
        xlims=(0, datat[n]),
        frame=:box)

    display(plt)
    readline()
 
    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/Lz_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    # Unbound particles

    plt = plot(datat[1:n], [dataUnbound[1:n] ./ Npart .* 100], 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel="Fraction of unbound stars [%]", 
        xlims=(0, datat[n]),
        # aspect_ratio=1,
        xticks=0:250:5000,
        yticks=0:5:100,
        frame=:box)

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/count_unbound_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    

    # Lagrange radii
    rh = tab_lag[1, 4] # Half-mass radius at t=0
    Trh = 0.138 * rh^(3/2) * Npart/log(0.11*Npart)


    plt = plot(datat[1:n] ./ Trh, [tab_lag[:, 1] tab_lag[:, 2] tab_lag[:, 3] tab_lag[:, 4] tab_lag[:, 5]], 
        labels=[L"r_{0.01}" L"r_{0.10}" L"r_{0.20}" L"r_{0.50}" L"r_{0.90}"], 
        xlabel="Time "*L"[T_{\mathrm{rh}}]", 
        ylabel="Lagrange radii", 
        yaxis=:log10,
        ylims=(0.01, 10.0),
        xlims=(0, datat[n]/ Trh),
        color=[:blue :orange :green :purple :red],
        # aspect_ratio=1,
        xticks=0:5:25,
        yticks=10.0 .^ (-2:1:1),
        frame=:box)

    # Add the softening radius
    eps_soft = 0.001
    plot!(plt, [0, datat[n]/ Trh], [eps_soft, eps_soft], linestyle=:dash, color=:black)#, label="Softening length")

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/lagrange_radii_"*srun*".pdf"
    savefig(plt, namefile_pdf)


    # Central density n_c(0)
    plt = plot(datat[1:n] ./ Trh, [tab_nc], 
        labels=:false, 
        xlabel="Time "*L"[T_{\mathrm{rh}}]", 
        ylabel="Central density "*L"n_c"*" [HU]", 
        yaxis=:log10,
        # ylims=(0.01, 10.0),
        xlims=(0, datat[n]./ Trh),
        color=:black,
        # aspect_ratio=1,
        xticks=0:5:25,
        # yticks=10.0 .^ (-2:1:1),
        frame=:box)

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/central_density_"*srun*".pdf"
    savefig(plt, namefile_pdf)
end 


plot_data!()