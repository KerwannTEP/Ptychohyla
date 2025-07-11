
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
    default = 0.012280754430794917
    "--run"
    help = "Run id"
    arg_type = Int64
    default = 63876778223743
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

const mass_avg = 1.0/Npart
const nb_neigh = 10

const R_HU_in_pc = R_HU_in_kpc * 1000.0

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

    tab_IOM = zeros(Float64, nsnap, 10) # time, E_tot, Lx, Ly, Lz, nb_unbound, E_wrt_cluster, Lx_cluster, Ly_cluster, Lz_cluster
    tab_lag = zeros(Float64, nsnap, 5) # 1% 10% 20% 50% 90% 
    tab_rhoc = zeros(Float64, nsnap) # Central density

    for isnap=1:nsnap

        println("Progress : ", round(isnap/nsnap, digits=4))

        namefile = sortedFiles[isnap]
        time = tabtsort[isnap]

        data_stars = readdlm(namefile) # x, y, z, vx, vy, vz, m, Uint, Uc


        # tab_vb_t = zeros(Float64, Threads.nthreads(), 3)
        K_t =  zeros(Float64, Threads.nthreads())
        Uc_t =  zeros(Float64, Threads.nthreads())
        Uh_t =  zeros(Float64, Threads.nthreads())
        L_t =  zeros(Float64, Threads.nthreads(), 3)
        

        Threads.@threads for i=1:Npart 

            tid = Threads.threadid()
            x, y, z, vx, vy, vz, m, Uint, Uc = data_stars[i, :]
            

            v2 = vx^2 + vy^2 + vz^2

            # # Barycenter 
            # tab_vb_t[tid, 1] += vx
            # tab_vb_t[tid, 2] += vy
            # tab_vb_t[tid, 3] += vz

            # Kinetic energy 
            K_t[tid] += 0.5 * m * v2 

            # Host potential energy
            Uh_t[tid] += Uc

            # Cluster potential energy 
            Uc_t[tid] += 0.5*Uint

                
            # Angular momenta 
            Lx = y*vz - z*vy 
            Ly = z*vx - x*vz 
            Lz = x*vy - y*vx 

            L_t[tid, 1] += m * Lx
            L_t[tid, 2] += m * Ly
            L_t[tid, 3] += m * Lz

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
            x, y, z, vx, vy, vz, m, Uint, Uc = data_stars[i, :]
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

        Etot_wrt_cluster_t = zeros(Float64, Threads.nthreads())
        L_wrt_cluster_t = zeros(Float64, Threads.nthreads(), 3)


        Threads.@threads for i=1:Npart

            tid = Threads.threadid()
    
            # Unbound particles 
            x, y, z, vx, vy, vz, m, Uint, Uc = data_stars[i, :]

       
            # Let x_c is the density center of the cluster.
            # It is a proxy for the cluster's center 
            # What about the velocity vc of that center ?
            # Calculation show that vc = dxc/dt + fluctuations 1/N

            vc_x = vx - Vcx
            vc_y = vy - Vcy
            vc_z = vz - Vcz
    
            Ec = 0.5 * m * (vc_x^2 + vc_y^2 + vc_z^2) + Uint

            x_c = x - Xc
            y_c = y - Yc
            z_c = z - Zc

            Etot_wrt_cluster_t[tid] += 0.5 * m * (vc_x^2 + vc_y^2 + vc_z^2) + 0.5*Uint + Uc 
            L_wrt_cluster_t[tid, 1] += m * (y_c*vc_z - z_c*vc_y)
            L_wrt_cluster_t[tid, 2] += m * (x_c*vc_z - z_c*vc_x)
            L_wrt_cluster_t[tid, 3] += m * (x_c*vc_y - y_c*vc_x)

            if (Ec >= 0.0)
                n_unbound_t[tid] += 1
            end
    
        end

        n_unbound = 0

        for tid=1:Threads.nthreads()
            n_unbound += n_unbound_t[tid]
            tab_IOM[isnap,7] += Etot_wrt_cluster_t[tid]
            tab_IOM[isnap,8] += L_wrt_cluster_t[tid, 1]
            tab_IOM[isnap,9] += L_wrt_cluster_t[tid, 2]
            tab_IOM[isnap,10] += L_wrt_cluster_t[tid, 3]
        end
    
        tab_IOM[isnap,1] = time
        tab_IOM[isnap,2] = Etot
        tab_IOM[isnap,3] = L[1]
        tab_IOM[isnap,4] = L[2]
        tab_IOM[isnap,5] = L[3]
        tab_IOM[isnap,6] = n_unbound

        # Lagrange radii 
        n_bound = Npart #-n_unbound
        tabr = zeros(Float64, n_bound)
        tabM = zeros(Float64, n_bound)
        index = 1

        Mbound = 0.0

        for i=1:Npart

            # Unbound particles 
            x, y, z, vx, vy, vz, m, Uint, Uc = data_stars[i, :]

            # Let x_c is the density center of the cluster.
            # It is a proxy for the cluster's center 
            # What about the velocity vc of that center ?
            # Calculation show that vc = dxc/dt + fluctuations 1/N

            vc_x = vx - Vcx
            vc_y = vy - Vcy
            vc_z = vz - Vcz
    
            Ec = 0.5 * m * (vc_x^2 + vc_y^2 + vc_z^2) + Uint
    
            # if (Ec < 0.0)
                r = sqrt((x-Xc)^2 + (y-Yc)^2 + (z-Zc)^2) 
                tabr[index] = r

                tabM[index] = m 
                Mbound += m

                index += 1
            # end
    
        end

        p = sortperm(tabr)
        tabr = tabr[p]
        tabM = tabM[p]

        m_enc = 0.0

        for index=1:Npart#n_bound

            r = tabr[index]
            m = tabM[index]

            m_enc += m

            tabM[index] = m_enc

            if (index == 1)
                
                # Initialize lagrange radii
                tab_lag[isnap, 1] = r 
                tab_lag[isnap, 2] = r 
                tab_lag[isnap, 3] = r 
                tab_lag[isnap, 4] = r 
                tab_lag[isnap, 5] = r 
            else

                if (tabM[index-1] < 0.01 * Mbound <= tabM[index])
                    tab_lag[isnap, 1] = r 
                end

                if (tabM[index-1] < 0.10 * Mbound <= tabM[index])
                    tab_lag[isnap, 2] = r 
                end

                if (tabM[index-1] < 0.20 * Mbound <= tabM[index])
                    tab_lag[isnap, 3] = r 
                end

                if (tabM[index-1] < 0.50 * Mbound <= tabM[index])
                    tab_lag[isnap, 4] = r 
                end

                if (tabM[index-1] < 0.90 * Mbound <= tabM[index])
                    tab_lag[isnap, 5] = r 
                end

            end


        end
    
        

        # tabr = sort(tabr)

        # index_1 = floor(Int64, n_bound * 0.01)
        # index_10 = floor(Int64, n_bound * 0.10)
        # index_20 = floor(Int64, n_bound * 0.20)
        # index_50 = floor(Int64, n_bound * 0.50)
        # index_90 = floor(Int64, n_bound * 0.90)

        # tab_lag[isnap, 1] = tabr[index_1]
        # tab_lag[isnap, 2] = tabr[index_10]
        # tab_lag[isnap, 3] = tabr[index_20]
        # tab_lag[isnap, 4] = tabr[index_50]
        # tab_lag[isnap, 5] = tabr[index_90]

        # index_01 = floor(Int64, n_bound * 0.01)
        M_01 = 0.01 * Mbound
        r_01 = tab_lag[isnap, 1]
        tab_rhoc[isnap] = M_01/(4.0/3.0*pi*r_01^3) # Number/Volume
        

    end


    return tab_IOM, tab_lag, tab_rhoc

end

function plot_data!()

    tab_IOM, tab_lag, tab_rhoc = get_data()
    datat = tab_IOM[:, 1]
    dataE = tab_IOM[:, 2]
    dataLx = tab_IOM[:, 3]
    dataLy = tab_IOM[:, 4]
    dataLz = tab_IOM[:, 5]
    dataUnbound = tab_IOM[:, 6]

    rh = tab_lag[1, 4] # Half-mass radius at t=0 
    Trh = 0.138 * rh^(3/2) * Npart/log(0.11*Npart)

    println("Trh [HU] = ", Trh)

    namefile_hf5 = path_data*"iom_cluster_"*srun*".hf5"
    file = h5open(namefile_hf5, "w")

    write(file, "data_time_HU", datat)
    write(file, "data_time_Myr", datat .* T_HU_in_Myr)
    write(file, "data_Etot", dataE)
    write(file, "data_Lx", dataLx)
    write(file, "data_Ly", dataLy)
    write(file, "data_Lz", dataLz)
    write(file, "data_unbound_frac", dataUnbound ./ Npart)

    write(file, "data_E_wrt_cluster_HU", tab_IOM[:, 7])
    write(file, "data_Lx_wrt_cluster_HU", tab_IOM[:, 8])
    write(file, "data_Ly_wrt_cluster_HU", tab_IOM[:, 9])
    write(file, "data_Lz_wrt_cluster_HU", tab_IOM[:, 10])

    write(file, "data_lagrange_rad_01_pc", tab_lag[:, 1].*R_HU_in_pc)
    write(file, "data_lagrange_rad_10_pc", tab_lag[:, 2].*R_HU_in_pc)
    write(file, "data_lagrange_rad_20_pc", tab_lag[:, 3].*R_HU_in_pc)
    write(file, "data_lagrange_rad_50_pc", tab_lag[:, 4].*R_HU_in_pc)
    write(file, "data_lagrange_rad_90_pc", tab_lag[:, 5].*R_HU_in_pc)

    write(file, "Npart", Npart)
    write(file, "kpc_per_HU", R_HU_in_kpc)
    write(file, "G_in_kpc_MSun_Myr", G_in_kpc_MSun_Myr)
    write(file, "Msun_per_HU", M_HU_in_Msun)

    write(file, "Myr_per_HU", T_HU_in_Myr)
    write(file, "Trh_HU", Trh)

    close(file)


    # https://stackoverflow.com/questions/68511668/how-to-set-number-of-minor-tick-marks-for-only-one-but-not-both-the-axes
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8

    # Energy
    n = length(dataE)
    dataFracE = abs.(1.0 .- dataE[2:n] ./ dataE[1])
    minE = 10.0^floor(Int64, max(log10(minimum(dataFracE)), -16))
    maxE = 10.0^(floor(Int64, log10(maximum(dataFracE)))+1)

    plt = plot(datat[2:n] ./ Trh, [dataFracE], 
        labels=:false, 
        xlabel="Time "*L"[T_{\mathrm{rh}}]", 
        ylabel="Fractional energy", 
        yaxis=:log10,
        # xticks=0:2:50,
        # xminorticks=4,
        yticks=10.0 .^ (-20:1:2),
        yminorticks=10,
        xlims=(0, datat[n]/Trh),
        ylims=(minE, maxE),
        frame=:box)

    # Workaround for x ticks on top
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twinx(plt),
        xlims=(0, datat[n]/Trh),
        # xticks=0:2:50,
        # xminorticks=4,

        yaxis=:log10,
        ylims=(minE, maxE),
        yticks=10.0 .^ (-20:1:2),
        yminorticks=10)

    # Workaround for y ticks on the right
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twiny(plt), 
        xlims=(0, datat[n]/Trh),
        # xticks=0:2:50,
        # xminorticks=4,

        yaxis=:log10,
        ylims=(minE, maxE),
        yticks=10.0 .^ (-20:1:2),
        yminorticks=10)

    display(plt)
    readline()
    

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/E_frac_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    # Angular momentum

    # x y
    plt = plot(datat ./Trh, [dataLx dataLy], 
        labels=[L"L_x" L"L_y"], 
        xlabel="Time "*L"[T_{\mathrm{rh}}]", 
        ylabel="Angular momentum", 
        # xticks=0:2:50,
        # xminorticks=4,
        # yaxis=:log10,
        # yticks=10.0 .^ (-15:1:0),
        xlims=(0, datat[n]/Trh),
        frame=:box)

    # # Workaround for x ticks on top
    # # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    # plot!(twinx(plt),
    #     xlims=(0, datat[n] .* T_HU_in_Myr),
    #     xticks=0:500:5000,
    #     xminorticks=5)

    # # Workaround for y ticks on the right
    # # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    # plot!(twiny(plt), 
    #     xlims=(0, datat[n] .* T_HU_in_Myr),
    #     xticks=0:500:5000,
    #     xminorticks=5)

    display(plt)
    readline()
 
    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/Lxy_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    # z 
    plt = plot(datat ./ Trh, [dataLz], 
        label=L"L_z", 
        xlabel="Time "*L"[T_{\mathrm{rh}}]", 
        ylabel="Angular momentum", 
        # xticks=0:2:50,
        # xminorticks=4,
        # yaxis=:log10,
        # yticks=10.0 .^ (-15:1:0),
        xlims=(0, datat[n]/Trh),
        frame=:box)

    # # Workaround for x ticks on top
    # # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    # plot!(twinx(plt),
    #     xlims=(0, datat[n] .* T_HU_in_Myr),
    #     xticks=0:500:5000,
    #     xminorticks=5)

    # # Workaround for y ticks on the right
    # # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    # plot!(twiny(plt), 
    #     xlims=(0, datat[n] .* T_HU_in_Myr),
    #     xticks=0:500:5000,
    #     xminorticks=5)

    display(plt)
    readline()
 
    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/Lz_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    # Unbound particles

    minu = 0.0

    maxu = min(100.0, 10.0*(1.0+floor(Int64, maximum(dataUnbound./ Npart .* 100)/10)))
    

    # plt = plot(datat[1:n] ./Trh, [dataUnbound[1:n] ./ Npart .* 100], 
    plt = plot(datat[1:n] .* T_HU_in_Myr, [dataUnbound[1:n] ./ Npart .* 100], 
        labels=:false, 
        # xlabel="Time "*L"[T_{\mathrm{rh}}]", 
        xlabel="Time [Myr]", 
        ylabel="Fraction of unbound stars [%]", 
        # xlims=(0, datat[n]/Trh),
        xlims=(0, datat[n].* T_HU_in_Myr),
        # aspect_ratio=1,
        # xticks=0:500:15000,
        # xticks=0:2:50,
        # xminorticks=4,
        # xminorticks=5,
        ylims=(minu, maxu),
        # ylims=(0, 32),
        yticks=0:10:100,
        yminorticks=5,
        gridline=:true,
        frame=:box)

        

    # # Workaround for x ticks on top
    # # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    # plot!(twinx(plt),
    #     xlims=(0, datat[n]/Trh),
    #     xticks=0:2:50,
    #     xminorticks=4,

    #     ylims=(minu, maxu),
    #     yticks=0:10:100,
    #     yminorticks=5)

    # # Workaround for y ticks on the right
    # # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    # plot!(twiny(plt), 
    #     xlims=(0, datat[n] /Trh),
    #     xticks=0:2:50,
    #     xminorticks=4,

    #     ylims=(minu, maxu),
    #     yticks=0:10:100,
    #     yminorticks=5)

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/count_unbound_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    

    # Lagrange radii
    
    minl = 10.0^floor(Int64, max(log10(minimum(tab_lag[:, 1].*R_HU_in_pc)), -16))
    maxl = 10.0^(floor(Int64, log10(maximum(tab_lag[:, 5].*R_HU_in_pc)))+1)

    plt = plot(datat[1:n] ./ Trh, [tab_lag[:, 1].*R_HU_in_pc tab_lag[:, 2].*R_HU_in_pc tab_lag[:, 3].*R_HU_in_pc tab_lag[:, 4].*R_HU_in_pc tab_lag[:, 5].*R_HU_in_pc], 
    # plt = plot(datat[1:n] .* T_HU_in_Myr, [tab_lag[:, 1].*R_HU_in_pc tab_lag[:, 2].*R_HU_in_pc tab_lag[:, 3].*R_HU_in_pc tab_lag[:, 4].*R_HU_in_pc tab_lag[:, 5].*R_HU_in_pc], 
        labels=[L"r_{0.01}" L"r_{0.10}" L"r_{0.20}" L"r_{0.50}" L"r_{0.90}"], 
        xlabel="Time "*L"[T_{\mathrm{rh}}]", 
        # xlabel="Time [Myr]", 
        ylabel="Lagrange radii [pc]", 
        yaxis=:log10,
        # ylims=(0.001, 100.0),
        # ylims=(10.0^(-2.0), 10.0^(2.0)),
        # ylims=(minl, maxl),
        # xlims=(0, datat[n] .* T_HU_in_Myr),
        xlims=(0, datat[n]/ Trh),
        color=[:blue :orange :green :purple :red],
        # aspect_ratio=1,
        # xticks=0:2:50,
        # xminorticks=4,
        # xticks=0:250:5000,
        yticks=10.0 .^ (-7:1:7),
        yminorticks=10,
        legend=:bottomleft,
        frame=:box)

    # Workaround for x ticks on top
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twinx(plt),
        xlims=(0, datat[n]/ Trh),
        # xticks=0:2:50,
        # xminorticks=4,

        yaxis=:log10,
        # ylims=(minl, maxl),
        yticks=10.0 .^ (-7:1:7),
        yminorticks=10)

    # Workaround for y ticks on the right
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twiny(plt), 
        xlims=(0, datat[n]/ Trh),
        # xticks=0:2:50,
        # xminorticks=4,

        yaxis=:log10,
        # ylims=(minl, maxl),
        yticks=10.0 .^ (-7:1:7),
        yminorticks=10)


    # Add the softening radius
    eps_soft = 0.001 
    plot!(plt, [0, datat[n]/ Trh], [eps_soft .*R_HU_in_pc, eps_soft .*R_HU_in_pc], linestyle=:dash, color=:black, label=:false)#, label="Softening length")

    # Add CC time
    tcc_trh = 16.45 # For an isolated King sphere with W0=5.0
    plot!(plt, [tcc_trh, tcc_trh], [minl, maxl], linestyle=:dot, color=:black, label=:false)

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/lagrange_radii_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    minrho = 10.0^floor(Int64, max(log10(minimum(tab_rhoc)), -16))
    maxrho = 10.0^(floor(Int64, log10(maximum(tab_rhoc)))+1)

    # Central density n_c(0)
    plt = plot(datat[1:n] ./ Trh, [tab_rhoc], 
    # plt = plot(datat[1:n] .* T_HU_in_Myr, [tab_rhoc], 
        labels=:false, 
        xlabel="Time "*L"[T_{\mathrm{rh}}]", 
        # xlabel="Time [Myr]", 
        ylabel="Central density "*L"\rho_c"*" [HU]", 
        yaxis=:log10,
        # xlims=(0, datat[n] .* T_HU_in_Myr),
        xlims=(0, datat[n]/ Trh),
        color=:black,
        # aspect_ratio=1,
        # xticks=0:2:50,
        # xminorticks=4,
        # xticks=0:250:5000,
        yticks=10.0 .^ (-2:1:10),
        yminorticks=10,
        ylims=(minrho, maxrho),
        frame=:box)

    # Workaround for x ticks on top
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twinx(plt),
        xlims=(0, datat[n]/ Trh),
        # xticks=0:2:50,
        # xminorticks=4,

        yaxis=:log10,
        ylims=(minrho, maxrho),
        yticks=10.0 .^ (-2:1:10),
        yminorticks=10)

    # Workaround for y ticks on the right
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twiny(plt), 
        xlims=(0, datat[n]/ Trh),
        # xticks=0:2:50,
        # xminorticks=4,

        yaxis=:log10,
        ylims=(minrho, maxrho),
        yticks=10.0 .^ (-2:1:10),
        yminorticks=10)

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/central_density_"*srun*".pdf"
    savefig(plt, namefile_pdf)
end 


plot_data!()