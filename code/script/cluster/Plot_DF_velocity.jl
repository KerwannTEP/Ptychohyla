using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse
using NearestNeighbors
using Statistics

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

# Conversion kpc/Myr to km/s (using astropy)
# x kpc/Myr = x kpc/km s/Myr km/s = y km/s ; y = x kpc/km s/Myr
# >>> l=1*u.kpc
# >>> l.to(u.km)
# <Quantity 3.08567758e+16 km>
# >>> t=u.Myr
# >>> t.to(u.s)
# <Quantity 3.15576e+13 s>
# lkm=l.to(u.km)
# ts=t.to(u.s)
# >>> lkm/ts
# <Quantity 977.79222168 km / s>

const V_HU_in_kpc_Myr = sqrt((G_in_kpc_MSun_Myr*M_HU_in_Msun)/R_HU_in_kpc)
const V_HU_in_km_s = V_HU_in_kpc_Myr * 977.79222168

const srun = string(run)
const nb_neigh = 6

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

    # nsnap = 1500
   
    namefile = listFile[p[end]]
    data = readdlm(namefile, header=false)
    interm = split(split(namefile,"_")[end],".")
    time = parse(Float64, interm[1]*"."*interm[2]) * T_HU_in_Myr
    time = round(time, digits=1)

    return data, time 
end

# Velocity norm in host frame
function get_unbound_velocities()

    data_stars, time = read_data()

    # Determine cluster's center (position and velocity)

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
        # x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]
        # m = 1.0/Npart 

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

    Vt_mean = sqrt(Vcx^2 + Vcy^2)

    println("Cluster's tangential velocity [km/s] = ", Vt_mean  * V_HU_in_km_s)

    tab_dens = zeros(Float64, Npart)

    # Compute number of unbound stars (i.e. stars in the stream)
    n_unbound_t = zeros(Int64, Threads.nthreads())

    Threads.@threads for i=1:Npart

        tid = Threads.threadid()

        # Unbound particles 
        x, y, z, vx, vy, vz, m, Uint, Uc = data_stars[i, :]
        # x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]
        # m = 1.0/Npart 

        # Let x_c is the density center of the cluster.
        # It is a proxy for the cluster's center 
        # What about the velocity vc of that center ?
        # Calculation show that vc = dxc/dt + fluctuations 1/N

        vc_x = vx - Vcx
        vc_y = vy - Vcy
        vc_z = vz - Vcz

        Ec = 0.5 * m * (vc_x^2 + vc_y^2 + vc_z^2) + Uint

        if (Ec >= 0.0)
            n_unbound_t[tid] += 1
        end

    end

    n_unbound = 0

    for tid=1:Threads.nthreads()
        n_unbound += n_unbound_t[tid]
    end

    tab_vel_stream = zeros(Float64, 2, n_unbound) # (v_r, v_t) in (x, y) plane
    
    # n = x/R e_x + y/R e_y
    # t = -y/R e_x + x/R e_y

    # v_X = <v,n> = x*vx + y*vy
    # v_Y = <v,t> = -y*vx + x*vy

    index = 1

    for i=1:Npart 
        x, y, z, vx, vy, vz, m, Uint, Uc = data_stars[i, :]
        # x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]
        # m = 1.0/Npart

        R = sqrt(x^2 + y^2)

        vc_x = vx - Vcx
        vc_y = vy - Vcy
        vc_z = vz - Vcz

        Ec = 0.5 * m * (vc_x^2 + vc_y^2 + vc_z^2) + Uint

        if (Ec >= 0.0)

            vr = x/R*vx + y/R*vy
            vt = -y/R*vx + x/R*vy   

            tab_vel_stream[1, index] = vr
            tab_vel_stream[2, index] = vt - Vt_mean

            index += 1
        end

    end


    tab_vel_stream *= V_HU_in_km_s

    return tab_vel_stream, time
end

function plot_data()

    tab_vel_stream, time = get_unbound_velocities()
    mean_vr = mean(tab_vel_stream[1, :])
    mean_vt = mean(tab_vel_stream[2, :])
    sigma_vr = sqrt(var(tab_vel_stream[1, :]))
    sigma_vt = sqrt(var(tab_vel_stream[2, :]))

    ns = 5
    ds = 2*ns/30

    println("Mean velocity = ", (mean_vr,mean_vt))
    println("Disp velocity = ", (sigma_vr,sigma_vt))

    p = histogram(tab_vel_stream[2, :], bins=mean_vt-ns*sigma_vt:ds:mean_vt+ns*sigma_vt,
                normalize=:density,
                xlabel=L"v_{\mathrm{t}} - v_{\mathrm{t}}^{\mathrm{cluster}}"*" [km/s]", ylabel="Number DF",
                framestyle=:box, label=:false,
                title="t = "*string(time)*" Myr")
                

    display(p)
    readline()


    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/cluster_"*srun*"_DF_vel_last_snapshot.pdf"
    savefig(p, namefile_pdf)


end

@time plot_data()

