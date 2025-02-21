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

const Mh = M_bulge + M_disk + Mvir
const M = Mh + _Mtot

const mu_h = Mh/M
const mu_c = _Mtot/M

const xh = -mu_c * d_host
const xc = mu_h * d_host



##########################################################################################
# Host potential (at all times)
# Coordinate w.r.t. barycenter
##########################################################################################

function psi_host(x::Float64, y::Float64)
    R = sqrt((x-xh)^2 + y^2)

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

    return rt, tabr, tab_Rc, tab_Vc

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

# Coordinate w.r.t. barycenter
function psi_cluster(x::Float64, y::Float64, tabr::Array{Float64}, tabpsi::Array{Float64})

    return psi_cluster_centered(sqrt((x-xc)^2 + y^2), tabr, tabpsi) 

end

# Coordinate w.r.t. barycenter
function psi_total(x::Float64, y::Float64, tabr::Array{Float64}, tabpsi::Array{Float64})

    return psi_host(x, y) + psi_cluster(x, y, tabr, tabpsi)

end

# Coordinate w.r.t. barycenter
function psi_eff(x::Float64, y::Float64, tabr::Array{Float64}, tabpsi::Array{Float64})

    vc = circular_velocity(d_host)

    Omegac = vc/d_host

    return psi_total(x, y, tabr, tabpsi) - 0.5*Omegac^2*(x^2+y^2)

end

function grad_phi_eff(x::Float64, y::Float64, tabr::Array{Float64}, tabpsi::Array{Float64}, eps::Float64=0.001)

    # x 
    dphidx = (psi_eff(x+eps, y, tabr, tabpsi) - psi_eff(x-eps, y, tabr, tabpsi))/(2.0*eps)

    # y
    dphidy = (psi_eff(x, y+eps, tabr, tabpsi) - psi_eff(x, y-eps, tabr, tabpsi))/(2.0*eps)

    return [dphidx, dphidy]

end

function hessian_phi_eff(x::Float64, y::Float64, tabr::Array{Float64}, tabpsi::Array{Float64}, eps::Float64=0.001)

    dPhidx_xp, _ =  grad_phi_ef(x+eps, y, tabr, tabpsi, eps)
    dPhidx_xm, _ =  grad_phi_ef(x-eps, y, tabr, tabpsi, eps)
    d2Phidx2 = (dPhidx_xp-dPhidx_xm)/(2.0*eps)

    dPhidx_yp, dPhidy_yp =  grad_phi_ef(x, y+eps, tabr, tabpsi, eps)
    dPhidx_ym, dPhidy_ym =  grad_phi_ef(x, y-eps, tabr, tabpsi, eps)
    d2Phidxy = (dPhidx_yp-dPhidx_ym)/(2.0*eps)
    d2Phidy2 = (dPhidy_yp-dPhidy_ym)/(2.0*eps)

    return [d2Phidx2 d2Phidxy; d2Phidxy d2Phidy2]

end

function norm2(v)
    return v[1]^2 + v[2]^2
end

function extremize_phi_eff(x0, y0, tabr::Array{Float64}, tabpsi::Array{Float64}, err=1.0*10^(-5), iterMAX=100, eps=0.001)

    r_guess = [x0, y0]
    gradPhi = grad_phi_eff(r_guess[1], r_guess[2], tabr, tabpsi, eps)

    iter = 0

    while (( norm2(gradPhi) > err^2) && (iter < iterMAX))

        Hess = hessian_phi_eff(r_guess[1], r_guess[2], tabr, tabpsi, eps)

        r_guess = r_guess - inv(Hess) * gradPhi
        gradPhi = grad_phi_eff(r_guess[1], r_guess[2], tabr, tabpsi, eps)

    end

    return r_guess
end



function compute_tE_snapshot(namefile)

    # namefile = path_data*"snapshots_"*srun_data*"/time_0.0.txt"

    rt, tabr, tab_Rc, tab_Vc = get_rt_tabr(namefile)
    tabpsi, tabM = get_tabpsi_tabM(tabr)
    nbound = length(tabr)

    nbx = 500

    tabx = range(xc-2*rt, xc, length=nbx)
    tabpsix = [psi_eff(tabx[i], 0.0, tabr, tabpsi) for i=1:nbx]

    # Compute L1, L2

    ixL1 = findmax(tabpsix[1:div(nbx,2)])[2]
    xL1 = tabx[ixL1] # HU
    yL1 = 0.0

    xL4 = (1/2-mu_c) * d_host
    yL4 = sqrt(3)/2  * d_host

    xL4, yL4 = extremize_phi_eff(xL4, yL4, tabr, tabpsi)

    # Load positions and rotate frame 
    Xc = tab_Rc[1]
    Yc = tab_Rc[2]
    VXc = tab_Vc[1]
    VYc = tab_Vc[2]
    theta = atan(Yc, Xc)

    data_stars = readdlm(namefile, header=false) # x, y, z, vx, vy, vz, Uint, Uc

    datax = data_stars[:, 1] .- Xc 
    datay = data_stars[:, 2] .- Yc

    dataxc = cos(-theta) .* datax - sin(-theta) .* datay
    datayc = sin(-theta) .* datax + cos(-theta) .* datay

    Rc = sqrt(Xc^2 + Yc^2)
    datax = dataxc .+ Rc
    datay = datayc

    psiEff1 = psi_eff(xL1, yL1, tabr, tabpsi)
    psiEff4 = psi_eff(xL4, yL4, tabr, tabpsi)

    Ecrit = psiEff1 - psiEff4

    EJ_t = zeros(Float64, Threads.nthreads())

    Threads.@threads for i=1:Npart 
        tid = Threads.threadid()

        x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]
        x = datax[i]
        y = datay[i]

        v2 = (vx-VXc)^2 + (vy-VYc)^2

        Ec = 0.5 * mass * v2 + Uint

        if (Ec < 0.0)
            EJ_i = 0.5*v2 + psi_eff(x, y, tabr, tabpsi) - psiEff4
            EJ_t[tid] += EJ_i
        end
    end

    EJ = 0.0
    for tid=1:Threads.nthreads()
        EJ += EJ_t[tid]
    end

    EJ /= nbound

    tE = (EJ - Ecrit)/abs(Ecrit)

    return tE
end

function compute_tE()

    if (isfile(path_data*"snapshots_"*srun_data*"/.DS_Store"))
        rm(path_data*"snapshots_"*srun_data*"/.DS_Store")
    end

    listFile = readdir(path_data*"snapshots_"*srun_data*"/";join=true)
    nsnap = length(listFile)

    tab_time = zeros(Float64, nsnap)

    for i=1:nsnap
        interm = split(split(listFile[i],"_")[end],".")
        index = parse(Float64, interm[1]*"."*interm[2])
        tab_time[i] = index
    end

    p = sortperm(tab_time)

    tabt = zeros(Float64, nsnap)
    tabtE = zeros(Float64, nsnap)

    for i=1:1:nsnap

        println("Progress = ", i/nsnap)
        namefile = listFile[p[i]]
        data = readdlm(namefile, header=false)
        interm = split(split(namefile,"_")[end],".")
        time = parse(Float64, interm[1]*"."*interm[2]) * T_HU_in_Myr

        tabt[i] = time
        tabtE[i] = compute_tE_snapshot(namefile)

    end

    plt = plot(tabt, [tabtE], 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel=L"\widetilde{E}", 
        xlims=(0, tabt[end]),
        color=:red,
        frame=:box)

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/cluster_"*srun*"_tilde_E.pdf"
    savefig(plt, namefile_pdf)

end

compute_tE()

# compute_Lagrange_points()

