using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse
using NearestNeighbors
using Glob 
using SphericalHarmonics
using SpecialFunctions

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
    "--framerate"
    help = "Number of frames per second"
    arg_type = Int64
    default = 20
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
const Npart0 = parsed_args["N"]

const path_to_script = @__DIR__
const path_data = path_to_script * "/../../../data/"

const srun = string(run)
const nb_neigh = 6

# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_MSun_Myr = 4.49851e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_MSun_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 

const mass_avg = 1.0/Npart0

const R_HU_in_pc = R_HU_in_kpc * 1000.0

####################################################################################################
# Import data of snapshot
####################################################################################################





function sorted_namefiles()

    # https://stackoverflow.com/questions/20484581/search-for-files-in-a-folder
    listFile = glob("time_*.txt", path_data*"snapshots_"*srun*"/")

    nsnap = length(listFile)

    tab_time = zeros(Float64, nsnap)

    for i=1:nsnap
        interm = split(split(listFile[i],"_")[end],".")
        index = parse(Float64, interm[1]*"."*interm[2])
        tab_time[i] = index
    end

    p = sortperm(tab_time)

    return listFile[p], tab_time[p]
end


function get_data_snapshot(namefile)

    data_stars = readdlm(namefile) # x, y, z, vx, vy, vz, m, Uint, Uc   (in HU)
    Npart = length(data_stars[:, 1])
    
    datax = data_stars[:, 1]
    datay = data_stars[:, 2]
    dataz = data_stars[:, 3]
    datam = data_stars[:, 7]

    tab_pos = zeros(Float64, 3, Npart)
    tab_dens = zeros(Float64, Npart)

    tab_pos[1, 1:Npart] = datax[1:Npart]
    tab_pos[2, 1:Npart] = datay[1:Npart] 
    tab_pos[3, 1:Npart] = dataz[1:Npart]

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

    Xc /= rho_tot
    Yc /= rho_tot
    Zc /= rho_tot
    Vcx /= rho_tot
    Vcy /= rho_tot
    Vcz /= rho_tot

    tab_pos[1, :] = tab_pos[1, :] .- Xc
    tab_pos[2, :] = tab_pos[2, :] .- Yc
    tab_pos[3, :] = tab_pos[3, :] .- Zc

    # Parse out unbound particles
    n_bound_t = zeros(Int64, Threads.nthreads())
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

        if (Ec < 0.0)
            n_bound_t[tid] += 1
        end

    end

    n_bound = 0

    for tid=1:Threads.nthreads()
        n_bound += n_bound_t[tid]
    end

    tab_pos_bound = zeros(Float64, 3, n_bound)
    datam_bound = zeros(Float64, n_bound)

    index = 1
    for i=1:Npart

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

        if (Ec < 0.0)
            tab_pos_bound[1, index] = tab_pos[1, i]
            tab_pos_bound[2, index] = tab_pos[2, i]
            tab_pos_bound[3, index] = tab_pos[3, i] 
            datam_bound[index] = m
            index += 1
        end

    end


    return tab_pos_bound, datam_bound, n_bound, Xc, Yc
    # return tab_pos, datam, Npart
end

function compute_rh_HU(tab_pos_bound, datam_bound)

    nbound = length(tab_pos_bound[1, :])
    mtot = 0.0

    for i=1:nbound 
        mtot += datam_bound[i]
    end

    tabr2 = tab_pos_bound[1, :] .^2 + tab_pos_bound[2, :] .^2 + tab_pos_bound[3, :] .^2
    p = sortperm(tabr2)

    rh = 0.0
    menc = 0.0
    for i=1:nbound 
        menc += datam_bound[p[i]]
        if (menc >= 0.5*mtot)
            r = sqrt(tabr2[p[i]])
            rh = r 
            break 
        end
    end

    return rh

end

####################################################################################################
# Zhao basis elements
####################################################################################################

const G = 1.0
const b_default = 1.0

function _xi(r::Float64, alpha::Float64, b_length::Float64=b_default)

    return ((r/b_length)^(1/alpha)-1.0)/((r/b_length)^(1/alpha)+1.0)

end

# https://juliapackages.com/p/sphericalharmonics
function Ylm(ell::Int64, m::Int64, theta::Float64, phi::Float64)

    Y = computeYlm(theta, phi, ell, m)
    return Y[(ell,m)]
end

# Precompute Y_{l}^{0}(theta_k, 0.0), which are real numbers
# Using Y_{l}^{m}(theta, phi) = Y_{l}^{0}(theta, 0.0) * exp(1im * m * phi)
# For l=0,1,...,lmax
function precompute_Ylm(lmax, pos, Npart)

    tab_Ylm = zeros(ComplexF64, Npart, lmax+1, lmax+1)

    Threads.@threads for i=1:Npart 
        x = pos[1,i]
        y = pos[2,i]
        z = pos[3,i]
        r = sqrt(x^2 + y^2 + z^2)
        theta = acos(z/r)
        phi = atan(y,x)

        Y = computeYlm(theta, phi, lmax)
        for l=0:lmax 
            for m=0:l
                tab_Ylm[i, l+1, m+1] = Y[(l,m)]
            end
        end
    end

    return tab_Ylm

end

# https://arxiv.org/pdf/2103.10165.pdf
function _Cnl(n::Int64, alpha::Float64, xi::Float64)

    if (n==0)
        return 1
    elseif (n==1)
        return 2*alpha*xi
    else
        return 2*(n-1+alpha)/n*xi*_Cnl(n-1,alpha,xi) - (n-2+2*alpha)/n*_Cnl(n-2,alpha,xi)
    end
end

function  _Knl(n::Int64, l::Int64, alpha::Float64)

    w = (2.0*l+1)*alpha + 0.5

    return (4.0*(n+w)^2-1.0)/(16.0*pi*alpha^2)
    
end

function _invNl(n::Int64, l::Int64, alpha::Float64)

    w = (2.0*l+1)*alpha + 0.5

    pref1 = 2^(4.0*w+1)*alpha*(n+w)/(pi*(4.0*(n+w)^2 - 1.0))
    pref2 = factorial(n)*gamma(w)^2/gamma(2.0*w+n)

    return pref1*pref2 

end

function _Unl(n::Int64, l::Int64, alpha::Float64, r::Float64, b_length::Float64=b_default)

    w = (2.0*l+1)*alpha + 0.5
    xi = _xi(r,alpha,b_length)

    pref = (r/b_length)^l / (1.0 + (r/b_length)^(1/alpha))^(alpha+2.0*l*alpha)
    invNl = _invNl(n,l,alpha)

    return -sqrt(G/b_length)*sqrt(4.0*pi)*sqrt(invNl)*pref*_Cnl(n,w,xi)
end

function _Dnl(n::Int64, l::Int64, alpha::Float64, r::Float64, b_length::Float64=b_default)

    w = (2.0*l+1)*alpha + 0.5
    xi = _xi(r,alpha,b_length)

    pref = (r/b_length)^(l-2+1/alpha) / (1.0 + (r/b_length)^(1/alpha))^(2+alpha+2*l*alpha)
    invNl = _invNl(n,l,alpha)
    Knl = _Knl(n,l,alpha)

    return 1/(sqrt(G)*b_length^(5/2)) * sqrt(4.0*pi)*Knl*sqrt(invNl)*pref*_Cnl(n,w,xi)
end

####################################################################################################
# Optimization
####################################################################################################

function _a00n(n::Int64, alpha::Float64, b_length::Float64, pos, tabm, Npart)

    ap = 0.0



    for i=1:Npart
        x = pos[1,i]
        y = pos[2,i]
        z = pos[3,i]
        mass = tabm[i]
        r = sqrt(x^2 + y^2 + z^2)

        ap -= mass*_Unl(n,0,alpha,r,b_length)*1.0/sqrt(4.0*pi)

    end

    return ap

end

function _almn(l::Int64, m::Int64, n::Int64, alpha::Float64, b_length::Float64, pos, tabm, Npart, tab_Ylm)

    ap_t = zeros(ComplexF64, Threads.nthreads())



    Threads.@threads for i=1:Npart 
        x = pos[1,i]
        y = pos[2,i]
        z = pos[3,i]
        mass = tabm[i]

        tid = Threads.threadid()

        # z = r cos theta
        # R = r sin theta

        r = sqrt(x^2 + y^2 + z^2)
        theta = acos(z/r)
        phi = atan(y,x)

        
        psip = _Unl(n,l,alpha,r,b_length) * tab_Ylm[i, l+1, m+1] 

        ap_t[tid] -= mass*conj(psip)
        
    end

    ap = 0.0 + 0.0im

    for tid=1:Threads.nthreads()
        ap += ap_t[tid]
    end

    return ap

end

function grad_b_a000k(r::Float64, alpha::Float64, b::Float64)

    pref = sqrt(gamma(3.0+2.0*alpha)/(8.0*b^3*gamma(2.0+alpha)^2))
    num = -1.0 + (r/b)^(1.0/alpha)
    den = (1.0 + (r/b)^(1.0/alpha))^(1.0+alpha)

    return pref * num/den 
end

function grad_alpha_a000k(r::Float64, alpha::Float64, b::Float64)

    # term 1
    pref1 = sqrt(gamma(3.0+2.0*alpha)/(b*gamma(2.0+alpha)^2))
    den1a = (1.0 + (r/b)^(1.0/alpha))^(alpha)
    num1a = (r/b)^(1.0/alpha)*log(r/b)
    den1b = (1.0 + (r/b)^(1.0/alpha))*(alpha)
    num1b = log(1+(r/b)^(1.0/alpha))

    term1 = pref1/sqrt(2.0)/den1a * (num1a/den1b - num1b)

    # term 2

    pref2a = -2.0*gamma(3.0+2.0*alpha)*digamma(2.0+alpha)/(b*gamma(2.0+alpha)^2)
    pref2b = 2.0*gamma(3.0+2.0*alpha)*digamma(3.0+2.0*alpha)/(b*gamma(2.0+alpha)^2)
    
    term2 = 1.0/den1a/(2.0*sqrt(2.0)) * (pref2a + pref2b)/pref1 

    return term1 + term2 
end


function grad_b_a000(alpha::Float64, b::Float64, pos, tabm, Npart)



    grad = zeros(Float64, Threads.nthreads())


    Threads.@threads for i=1:Npart
        x = pos[1,i]
        y = pos[2,i]
        z = pos[3,i]
        mass = tabm[i]
        r = sqrt(x^2 + y^2 + z^2)

        tid = Threads.threadid()
        grad[tid] += mass*grad_b_a000k(r,alpha,b)

    end

    grad_tot = 0.0

    for tid=1:Threads.nthreads()
        grad_tot += grad[tid]
    end

    return grad_tot
end


function grad_alpha_a000(alpha::Float64, b::Float64, pos, tabm, Npart)



    grad = zeros(Float64, Threads.nthreads())


    Threads.@threads for i=1:Npart
        x = pos[1,i]
        y = pos[2,i]
        z = pos[3,i]
        mass = tabm[i]
        r = sqrt(x^2 + y^2 + z^2)

        tid = Threads.threadid()

        grad[tid] += mass*grad_alpha_a000k(r,alpha,b)

    end

    grad_tot = 0.0

    for tid=1:Threads.nthreads()
        grad_tot += grad[tid]
    end

    return grad_tot
end


# https://en.wikipedia.org/wiki/Gradient_descent
function optimize_a000(pos, tabm, Npart, alpha_start::Float64=0.5, b_start::Float64=0.5, nsteps::Int64=100, err::Float64=1.0*10^(-10))



    alpha = alpha_start
    b = b_start
    logb = log(b)

    # df/db = df/dlnb dlnb/db = df/dlnb 1/b
    # df/dlnb = b df/db

    # Initialize
    a000 = _a00n(0,alpha,b,pos,tabm, Npart)



    grad_alpha = grad_alpha_a000(alpha,b, pos,tabm,  Npart)
    grad_b = grad_b_a000(alpha,b,pos,tabm, Npart)
    grad_logb = b * grad_b

    gamma_n = 1.0

    delta_alpha =  gamma_n * grad_alpha
    delta_logb =  gamma_n * grad_logb 

    if (alpha + delta_alpha < 0.0)
        alpha /= 2.0
    else 
        alpha += delta_alpha
    end

    
    logb += delta_logb
    b = exp(logb)



    a000 = _a00n(0,alpha,b,pos,tabm, Npart)



    delta = sqrt(delta_alpha^2 + delta_logb^2)

    step = 1

    while (step <= nsteps && delta > err)

        grad_alpha_old = grad_alpha
        grad_logb_old = grad_logb 

        grad_alpha = grad_alpha_a000(alpha,b,pos,tabm, Npart)
        grad_b = grad_b_a000(alpha,b,pos,tabm, Npart)
        grad_logb = b * grad_b

        normSq = (grad_alpha-grad_alpha_old)^2 + (grad_logb-grad_logb_old)^2

        gamma_n = abs(delta_alpha * (grad_alpha-grad_alpha_old)+ delta_logb * (grad_logb-grad_logb_old))/normSq 

        delta_alpha =  gamma_n * grad_alpha
        delta_logb =  gamma_n * grad_logb 

        if (alpha + delta_alpha < 0.0)
            alpha /= 2.0
        else 
            alpha += delta_alpha
        end

        logb += delta_logb
        b = exp(logb)

        a000 = _a00n(0,alpha,b,pos,tabm, Npart)



        delta = sqrt(delta_alpha^2 + delta_logb^2)

        step += 1


    end



    return alpha, b, a000
end

####################################################################################################
# Fourier analysis
####################################################################################################

function fourier_snapshot(namefile::String, nmax::Int64=10, lmax::Int64=4)

    pos, tabm, Npart, Xc, Yc = get_data_snapshot(namefile)
    tab_Ylm = precompute_Ylm(lmax, pos, Npart)

    tab_E = zeros(Float64, lmax+1)

    alphaZ, bZ, a000Z =  optimize_a000(pos, tabm, Npart)

    for n=0:nmax 
        for l=0:lmax

            # m = 0
            ap = _almn(l, 0, n, alphaZ, bZ, pos, tabm, Npart, tab_Ylm)
            tab_E[l+1] += abs2(ap)

            # m>0 : |a_{l-mn}| = |a_{lmn}|
            for m=1:l
                ap = _almn(l, m, n, alphaZ, bZ, pos, tabm, Npart, tab_Ylm)


                tab_E[l+1] += 2.0 * abs2(ap)
            end
        end
    end

    theta = atan(Yc, Xc)

    dataxc = cos(-theta) .* pos[1, :] - sin(-theta) .* pos[2, :]
    datayc = sin(-theta) .* pos[1, :] + cos(-theta) .* pos[2, :]

    return tab_E, dataxc, datayc, alphaZ, bZ, pos, tabm, Npart
end

function plot_fourier_data(nmax::Int64=10, lmax::Int64=4)

    listFiles, tabt = sorted_namefiles()
    nsnap = length(listFiles)

    # nsnap = 100

    list_E = zeros(Float64, nsnap, lmax+1)
    list_alpha = zeros(Float64, nsnap)
    list_b = zeros(Float64, nsnap)
    list_avg_m = zeros(Float64, nsnap)

    rh = 0.0
    Trh = 0.0

    anim = @animate for i=1:nsnap 
        println("Time [HU] = ", tabt[i])
        namefile = listFiles[i]
        time = round(tabt[i], digits=1)

        tab_E, dataxc, datayc, alphaZ, bZ, pos, tabm, Npart = fourier_snapshot(namefile, nmax, lmax)
       
        if (i==1)
            rh = compute_rh_HU(pos, tabm)
            Trh = 0.138 * rh^(3/2) * Npart0/log(0.11*Npart0)
            println("Trh [HU] = ", Trh)
        end

        list_alpha[i] = alphaZ
        list_b[i] = bZ
        list_avg_m[i] = sum(tabm) / Npart * M_HU_in_Msun

        # Normalized Energy 
        for l=0:lmax
            list_E[i, l+1] = tab_E[l+1]/tab_E[1]
        end

        rmax = 10.0 # HU
        s = 1.0



        # https://docs.juliaplots.org/latest/generated/attributes_plot/
        # https://stackoverflow.com/questions/71992758/size-and-colour-in-julia-scatter-plot

        # plot([0, -rmax/2.0],[0.0,0.0 ],
        #         xlabel=L"x"*" [kpc]", ylabel=L"y"*" [kpc]", 
        #         framestyle=:box, labels=:false,
        #         xlims=(-rmax, rmax), ylims=(-rmax,rmax), 
        #         aspect_ratio=1, size=(800,800), 
        #         # left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
        #         background_color = :black,
        #         #markersize=s, color=:white, 
        #         title="t = "*string(time)*" Myr",
        #         arrow=true,color=:red,linewidth=2,label=:false)

        scatter(dataxc, datayc, 
                xlabel=L"x"*" [HU]", ylabel=L"y"*" [HU]", 
                framestyle=:box, 
                labels=:false,
                xlims=(-rmax, rmax), ylims=(-rmax,rmax), 
                aspect_ratio=1, size=(800,800), 
                left_margin = [2mm 0mm], right_margin = [2mm 0mm], 
                background_color = :black,
                markersize=s, color=:white, 
                title="t = "*string(time)*" HU")



    end

    mkpath(path_data*"gif/")
    namefile_gif = path_data*"gif/cluster_bound_"*srun*"_centered.gif"
    gif(anim, namefile_gif, fps = framepersec)





    # https://discourse.julialang.org/t/there-is-a-problem-with-minimum-function-can-any-one-help-me/87628/2
    minE1 = 10.0^floor(Int64, log10(minimum(x for x ∈ list_E[:, 2] if !isnan(x))))
    minE2 = 10.0^floor(Int64, log10(minimum(x for x ∈ list_E[:, 3] if !isnan(x))))
    minE3 = 10.0^floor(Int64, log10(minimum(x for x ∈ list_E[:, 4] if !isnan(x))))
    minE4 = 10.0^floor(Int64, log10(minimum(x for x ∈ list_E[:, 5] if !isnan(x))))
    minE = min(minE1, minE2, minE3, minE4)

    maxE1 = 10.0^(floor(Int64, log10(maximum(x for x ∈ list_E[:, 2] if !isnan(x)))) + 1)
    maxE2 = 10.0^(floor(Int64, log10(maximum(x for x ∈ list_E[:, 3] if !isnan(x)))) + 1)
    maxE3 = 10.0^(floor(Int64, log10(maximum(x for x ∈ list_E[:, 4] if !isnan(x)))) + 1)
    maxE4 = 10.0^(floor(Int64, log10(maximum(x for x ∈ list_E[:, 5] if !isnan(x)))) + 1)
    maxE = max(maxE1, maxE2, maxE3, maxE4)


    pabs = plot(tabt ./ Trh, [list_E[:, 2], list_E[:, 3], list_E[:, 4], list_E[:, 5]] , 
                xlabel="Time "*L" [T_{\mathrm{rh}}]",
                ylabel=L"E_{\ell}\ / \ E_{0} ",
                labels=[L"\ell=1" L"\ell=2" L"\ell=3" L"\ell=4"], 
                linestyle=:solid,
                framestyle=:box, 
                legends=:bottomright,
                xlims=(0, tabt[nsnap]/ Trh),
                # legends=:topright,
                yaxis=:log10,
                ylims=(minE, maxE),
                xticks=0:2:25,
                xminorticks=4,
                yticks=10.0 .^ (-6:1:0),
                left_margin = 5mm, right_margin = 5mm, 
                linecolor=[:cyan :red :green :magenta],
                # size=(600, 300),
                bottom_margin = 5mm)

    # Workaround for x ticks on top
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twinx(pabs),
        xlims=(0, tabt[nsnap]/ Trh),
        xticks=0:2:25,
        xminorticks=4,

        yaxis=:log10,
        ylims=(minE, maxE),
        yticks=10.0 .^ (-7:1:3),
        yminorticks=10)

    # Workaround for y ticks on the right
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twiny(pabs), 
        xlims=(0, tabt[nsnap]/ Trh),
        xticks=0:2:25,
        xminorticks=4,

        yaxis=:log10,
        ylims=(minE, maxE),
        yticks=10.0 .^ (-7:1:3),
        yminorticks=10)

    display(pabs)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/harmonics_energy_"*srun*".pdf"
    savefig(pabs, namefile_pdf)




    minalpha = 0.1*(floor(Int64, minimum(x for x ∈ list_alpha if !isnan(x))/0.1))
    maxalpha = 0.1*(1.0+floor(Int64, maximum(x for x ∈ list_alpha if !isnan(x))/0.1))

    plt = plot(tabt ./ Trh, [list_alpha] , 
                xlabel="Time "*L" [T_{\mathrm{rh}}]",
                ylabel=L"\alpha",
                labels=:false,
                linestyle=:solid,
                framestyle=:box,
                xlims=(0, tabt[nsnap]/ Trh),
                ylims=(minalpha, maxalpha),
                xticks=0:2:25,
                xminorticks=4,
                yticks=(0:0.2:5),
                yminorticks=4)
    
    # Workaround for x ticks on top
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twinx(plt),
            xlims=(0, tabt[nsnap]/ Trh),
            xticks=0:2:25,
            xminorticks=4,

            ylims=(minalpha, maxalpha),
            yticks=(0:0.2:5),
            yminorticks=4)

    # Workaround for y ticks on the right
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twiny(plt), 
            xlims=(0, tabt[nsnap]/ Trh),
            xticks=0:2:25,
            xminorticks=4,

            ylims=(minalpha, maxalpha),
            yticks=(0:0.2:5),
            yminorticks=4)


    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/evolution_alpha_"*srun*".pdf"
    savefig(plt, namefile_pdf)






    minb = 0.1*(floor(Int64, minimum(x for x ∈ list_b if !isnan(x))/0.1))
    maxb = 0.1*(1.0+floor(Int64, maximum(x for x ∈ list_b if !isnan(x))/0.1))


    plt = plot(tabt ./ Trh, [list_b] , 
                xlabel="Time "*L" [T_{\mathrm{rh}}]",
                ylabel=L"b",
                labels=:false,
                linestyle=:solid,
                framestyle=:box,
                xlims=(0, tabt[nsnap]/ Trh),
                ylims=(minb, maxb),
                xticks=0:2:25,
                xminorticks=4,
                yticks=(0:0.1:4),
                yminorticks=4)

    # Workaround for x ticks on top
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twinx(plt),
            xlims=(0, tabt[nsnap]/ Trh),
            xticks=0:2:25,
            xminorticks=4,

            ylims=(minb, maxb),
            yticks=(0:0.1:4),
            yminorticks=4)

    # Workaround for y ticks on the right
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twiny(plt), 
            xlims=(0, tabt[nsnap]/ Trh),
            xticks=0:2:25,
            xminorticks=4,

            ylims=(minb, maxb),
            yticks=(0:0.1:4),
            yminorticks=4)


    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/evolution_b_"*srun*".pdf"
    savefig(plt, namefile_pdf)




    minm = minimum(x for x ∈ list_avg_m if !isnan(x))
    maxm = maximum(x for x ∈ list_avg_m if !isnan(x))

    # https://stackoverflow.com/questions/73129565/how-to-round-to-the-next-largest-integer-in-julia
    minm = 0.5*(round(minm/0.5)-1)
    maxm = 0.5*(round(maxm/0.5)+1)

    plt = plot(tabt ./ Trh, [list_avg_m] , 
                xlabel="Time "*L" [T_{\mathrm{rh}}]",
                ylabel=L"\langle m \rangle \ [M_{\odot}]",
                labels=:false,
                linestyle=:solid,
                framestyle=:box,
                xlims=(0, tabt[nsnap]/ Trh),
                ylims=(minm, maxm),
                xticks=0:2:25,
                xminorticks=4,
                yticks=(0:0.5:20),
                yminorticks=5)

    # Workaround for x ticks on top
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twinx(plt),
            xlims=(0, tabt[nsnap]/ Trh),
            xticks=0:2:25,
            xminorticks=4,

            ylims=(minm, maxm),
            yticks=(0:0.5:20),
            yminorticks=5)

    # Workaround for y ticks on the right
    # https://discourse.julialang.org/t/plot-ticks-at-both-top-and-bottom-axis/9550/8
    plot!(twiny(plt), 
            xlims=(0, tabt[nsnap]/ Trh),
            xticks=0:2:25,
            xminorticks=4,

            ylims=(minm, maxm),
            yticks=(0:0.5:20),
            yminorticks=5)


    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/evolution_avg_m_"*srun*".pdf"
    savefig(plt, namefile_pdf)




end

plot_fourier_data(10, 4)
