
using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse
using NearestNeighbors
using LaTeXStrings
using Statistics

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
    default = 63875411207673
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

    namefile = sortedFiles[end]

    time = round(tabtsort[end] * T_HU_in_Myr, digits=1)
    data_stars = readdlm(namefile) # x, y, z, vx, vy, vz, Uint, Uc

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

    tab_dens_vel_t = zeros(Float64, Threads.nthreads(), 3)
    tab_dens_t =  zeros(Float64, Threads.nthreads())
    

    Threads.@threads for i=1:Npart 

        tid = Threads.threadid()
        x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]
        rho = tab_dens[i]
        
        tab_dens_vel_t[tid, 1] += vx * rho
        tab_dens_vel_t[tid, 2] += vy * rho
        tab_dens_vel_t[tid, 3] += vz * rho

        tab_dens_t[tid] += rho

        
    end

    rho_tot = 0.0
    Vcx = 0.0
    Vcy = 0.0
    Vcz = 0.0

    for tid=1:Threads.nthreads()

        Vcx += tab_dens_vel_t[tid, 1]
        Vcy += tab_dens_vel_t[tid, 2]
        Vcz += tab_dens_vel_t[tid, 3]

        rho_tot += tab_dens_t[tid]


    end

    Vcx /= rho_tot
    Vcy /= rho_tot
    Vcz /= rho_tot

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

    tab_IOM_unbound = zeros(Float64, n_unbound, 2) # E_wrt_host, Lz

    index = 1
    for i=1:Npart 
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
            E =  0.5 * mass * (vx^2 + vy^2 + vz^2) + Uint + Uc
            Lz = x*vy - y*vx 
            tab_IOM_unbound[index, 1] = E
            tab_IOM_unbound[index, 2] = Lz 

            index += 1
        end
    end

    mean_E = mean(tab_IOM_unbound[:, 1])
    var_E = var(tab_IOM_unbound[:, 1], corrected=true, mean=mean_E)
    sigma_E = sqrt(var_E)

    mean_Lz = mean(tab_IOM_unbound[:, 2])
    var_Lz = var(tab_IOM_unbound[:, 2], corrected=true, mean=mean_Lz)
    sigma_Lz = sqrt(var_Lz)

    tab_IOM_unbound[:, 1] = (tab_IOM_unbound[:, 1] .- mean_E) ./ sigma_E
    tab_IOM_unbound[:, 2] = (tab_IOM_unbound[:, 2] .- mean_Lz) ./ sigma_Lz 

    s = 1.0
    rmax = 2
    plt = scatter(tab_IOM_unbound[:, 2], tab_IOM_unbound[:, 1],
                xlabel=L"\Delta L_z", ylabel=L"\Delta E", 
                framestyle=:box, labels=:false, 
                xlims=(-rmax, rmax), ylims=(-rmax,rmax), 
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                xticks=-3:0.5:3, yticks=-3:0.5:3,
                title="t = "*string(time)*" Myr")

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/cluster_"*srun*"_last_IOM_spread.pdf"
    savefig(plt, namefile_pdf)
end

get_data()