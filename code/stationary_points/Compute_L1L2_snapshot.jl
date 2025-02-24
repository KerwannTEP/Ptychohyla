include("../src/king/Args.jl")
include("../src/king/Constants.jl")
include("../src/Host.jl")
include("../src/king/Cluster.jl")

using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse
using NearestNeighbors
using LaTeXStrings

# z = 0

const run_data = 63875411207673
const srun_data = string(run_data)

const path_to_script = @__DIR__
const path_data = path_to_script * "/../../data/"

const nb_neigh = 10

##########################################################################################
# Host potential (at all times)
##########################################################################################

function psi_host(x::Float64, y::Float64)
    R = sqrt(x^2 + y^2)

    return psi_halo(R) + psi_bulge(R) + psi_disk(R, 0.0)

end


##########################################################################################
# Compute potential for each stars of the bound cluster (at given snapshot)
##########################################################################################

function get_rt_tabr(namefile::String)

    data_stars = readdlm(namefile, header=false) # x, y, z, vx, vy, vz, Uint, Uc

    # Compute density center (position and velocity)

    tab_pos = zeros(Float64, 3, Npart)
    tab_dens = zeros(Float64, Npart)

    tab_pos[1, :] = data_stars[:, 1]
    tab_pos[2, :] = data_stars[:, 2]
    tab_pos[3, :] = data_stars[:, 3]

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

    
    tab_Rc_t = zeros(Float64, Threads.nthreads(), 3)
    tab_Vc_t = zeros(Float64, Threads.nthreads(), 3)
    tab_dens_t =  zeros(Float64, Threads.nthreads())
    
    Threads.@threads for i=1:Npart 

        tid = Threads.threadid()
        x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]
        rho = tab_dens[i]
        
        tab_Rc_t[tid, 1] += x * rho
        tab_Rc_t[tid, 2] += y * rho
        tab_Rc_t[tid, 3] += z * rho

        tab_Vc_t[tid, 1] += vx * rho
        tab_Vc_t[tid, 2] += vy * rho
        tab_Vc_t[tid, 3] += vz * rho

        tab_dens_t[tid] += rho
        
    end

    rho_tot = 0.0
    tab_Rc = zeros(Float64, 3)
    tab_Vc = zeros(Float64, 3)

    for tid=1:Threads.nthreads()

        tab_Rc[1] += tab_Rc_t[tid, 1]
        tab_Rc[2] += tab_Rc_t[tid, 2]
        tab_Rc[3] += tab_Rc_t[tid, 3]

        tab_Vc[1] += tab_Vc_t[tid, 1]
        tab_Vc[2] += tab_Vc_t[tid, 2]
        tab_Vc[3] += tab_Vc_t[tid, 3]

        rho_tot += tab_dens_t[tid]

    end

    tab_Rc = tab_Rc ./ rho_tot
    tab_Vc = tab_Vc ./ rho_tot


    # Compute number of bound stars
    n_bound_t = zeros(Int64, Threads.nthreads())

    Threads.@threads for i=1:Npart

        tid = Threads.threadid()

        # Unbound particles 
        x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]

        # Let x_c is the density center of the cluster.
        # It is a proxy for the cluster's center 
        # What about the velocity vc of that center ?
        # Calculation show that vc = dxc/dt + fluctuations 1/N

        vc_x = vx - tab_Vc[1]
        vc_y = vy - tab_Vc[2]
        vc_z = vz - tab_Vc[3]

        Ec = 0.5 * mass * (vc_x^2 + vc_y^2 + vc_z^2) + Uint

        if (Ec < 0.0)
            n_bound_t[tid] += 1
        end

    end

    n_bound = 0

    for tid=1:Threads.nthreads()
        n_bound += n_bound_t[tid]
    end

    tabr = zeros(Float64, n_bound)

    index = 1
    for i=1:Npart 
        x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]

        vc_x = vx - tab_Vc[1]
        vc_y = vy - tab_Vc[2]
        vc_z = vz - tab_Vc[3]

        Ec = 0.5 * mass * (vc_x^2 + vc_y^2 + vc_z^2) + Uint

        if (Ec < 0.0)
            x_c = x - tab_Rc[1]
            y_c = y - tab_Rc[2]
            z_c = z - tab_Rc[3]

            tabr[index] = sqrt(x_c^2 + y_c^2 + z_c^2)

            index += 1
        end

    end

    tabr = sort(tabr)

    rt = tabr[n_bound]

    return rt, tabr

end

function get_tabpsi_tabM(tabr::Array{Float64})

    nbound = length(tabr)
    
    tabM = zeros(Float64, nbound)
    tabpsi = zeros(Float64, nbound)

    tabM[nbound] = _Mtot * nbound/Npart
    tabpsi[nbound] = -_G*tabM[nbound]/tabr[nbound]

    for i=nbound-1:-1:1

        tabM[i] = tabM[i+1] - mass 
        tabpsi[i] = tabpsi[i+1] - _G * tabM[i] * (1.0/tabr[i] - 1.0/tabr[i+1])

    end

    return tabpsi, tabM

end



##########################################################################################
# Cluster potential
##########################################################################################

# r = distance to center
# Henon 1971 (https://link.springer.com/article/10.1007/BF00649201)
# psi(r) = psi_{k+1} + (1/r_{k} - 1/r)/(1/r_{k} - 1/r_{k+1}) * (psi_{k+1} - psi_{k}) for 
function psi_cluster_centered(r::Float64, tabr::Array{Float64}, tabpsi::Array{Float64})

    nbound = length(tabr)

    if (r <= tabr[1])
        return tabpsi[1]

    elseif (r >= tabr[nbound])
        return tabpsi[nbound]*tabr[nbound]/r

    else # r_{1} < r < r_{nbound}

        # Find k by bisection 
        k_left = 1
        k_right = nbound
        k_mid = div(k_left + k_right, 2)

        r_left = tabr[k_left]
        r_right = tabr[k_right]
        r_mid = tabr[k_mid]

        

        while (k_right - k_left != 1)

            if (r_left <= r < r_mid)
                k_right = k_mid
                k_mid = div(k_left + k_right, 2)
                r_right = tabr[k_right]
                r_mid = tabr[k_mid]
            else
                k_left = k_mid
                k_mid = div(k_left + k_right, 2)
                r_left = tabr[k_left]
                r_mid = tabr[k_mid]
            end

        end

        psi_left = tabpsi[k_left]
        psi_right = tabpsi[k_right]

        psi = psi_left + (1.0/r_left-1.0/r)/(1.0/r_left-1.0/r_right) * (psi_right-psi_left)
        return psi

    end

end

# Valid for particles outside of cluster 
function psi_cluster(x::Float64, y::Float64, tabr::Array{Float64}, tabpsi::Array{Float64})

    return psi_cluster_centered(sqrt((x-d_host)^2 + y^2), tabr, tabpsi) #-_G * _Mtot/sqrt((x-d_host)^2 + y^2)

end

function psi_total(x::Float64, y::Float64, tabr::Array{Float64}, tabpsi::Array{Float64})

    return psi_host(x, y) + psi_cluster(x, y, tabr, tabpsi)

end

function psi_eff(x::Float64, y::Float64, tabr::Array{Float64}, tabpsi::Array{Float64})

    vc = circular_velocity(d_host)

    Omegac = vc/d_host

    return psi_total(x, y, tabr, tabpsi) - 0.5*Omegac^2*(x^2+y^2)

end


##########################################################################################
# Plot data
##########################################################################################

function plot_data!()

    namefile = path_data*"snapshots_"*srun_data*"/time_1000.0.txt"

    rt, tabr = get_rt_tabr(namefile)
    tabpsi, tabM = get_tabpsi_tabM(tabr)

    nbx = 1000
    rmax = 0.5 # kpc

    tabx = range(d_host-rmax/R_HU_in_kpc, d_host+rmax/R_HU_in_kpc, length=nbx)
    tabpsix = [psi_eff(tabx[i], 0.0, tabr, tabpsi) for i=1:nbx]

    ixL1 = findmax(tabpsix[1:div(nbx,2)])[2]
    ixL2 = findmax(tabpsix[div(nbx,2):nbx])[2] + div(nbx,2) - 1

    minpsi = minimum(tabpsix)
    maxpsi = maximum(tabpsix) + 2

    # Plot potential
    p = plot((tabx .- d_host) .* R_HU_in_kpc, tabpsix, 
        xlabel=L"x-x_{\mathrm{c}}"*" [kpc]", ylabel=L"\psi_{\mathrm{eff}} (X,0)",
        ylims=(minpsi, maxpsi),
        # xticks=-20:1:20,
        label=:false,
        color=:black,
        size=((600, 300)),
        frame=:box)

    # Plot Lagrange radii
    scatter!(p, ([tabx[ixL1]] .- d_host) .* R_HU_in_kpc,  [tabpsix[ixL1]], color=:red, label=L"\mathrm{L}_1", ylims=(minpsi, maxpsi))
    scatter!(p, ([tabx[ixL2]] .- d_host) .* R_HU_in_kpc,  [tabpsix[ixL2]], color=:blue, label=L"\mathrm{L}_2", ylims=(minpsi, maxpsi))

    # Plot rt
    plot!(p, [-rt, -rt].* R_HU_in_kpc, [-10000, 10000], linestyle=:dash, color=:blue, ylims=(minpsi, maxpsi), label=L"r_t")
    plot!(p, [rt, rt].* R_HU_in_kpc, [-10000, 10000], linestyle=:dash, color=:blue, label=:false, ylims=(minpsi, maxpsi))

    display(p)
    readline()

    p = plot(tabr, [tabpsi],
        xlabel=L"r"*" [HU]", ylabel=L"\psi(r)",
        # ylims=(minpsi, maxpsi),
        # xticks=-20:1:20,
        label=:false,
        color=:black,
        xaxis=:log10,
        size=((600, 300)),
        frame=:box)

    eps_soft = 0.01
    plot!(p, [eps_soft, eps_soft], [tabpsi[1], tabpsi[end]], linestyle=:dash, color=:blue, label="Softening length")


    display(p)
    readline()

    p = plot(tabr, [tabM],
        xlabel=L"r"*" [HU]", ylabel=L"M(r)",
        # ylims=(minpsi, maxpsi),
        # xticks=-20:1:20,
        label=:false,
        color=:black,
        xaxis=:log10,
        size=((600, 300)),
        frame=:box)

    plot!(p, [eps_soft, eps_soft], [0, 1], linestyle=:dash, color=:blue, label="Softening length")


    display(p)
    readline()

end

plot_data!()
