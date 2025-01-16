using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse


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
    default = 2.00e-2
    "--c"
    help = "Concentration parameter of the NFW dark halo profile"
    arg_type = Float64
    default = 9.4
    "--Rs_host"
    help = "Scale radius of the NFW dark halo (in kpc)"
    arg_type = Float64
    default = 23.8
    "--M_halo_200"
    help = "Mass of the dark halo, M200, within a sphere whose mean density is 200 times the critical density of the universe (in solar masses)"
    arg_type = Float64
    default = 0.97e+12
    "--framerate"
    help = "Number of frames per second"
    arg_type = Int64
    default = 60

end
parsed_args = parse_args(tabargs)

const Mtot_Msun = parsed_args["M_cluster"]
const Rv_kpc = parsed_args["Rv_cluster"]

const c = parsed_args["c"]
const M_DH_200_Msun = parsed_args["M_halo_200"]
const Rs_kpc = parsed_args["Rs_host"]

const framepersec = parsed_args["framerate"]

const path_to_script = @__DIR__
const path_data = path_to_script * "/../../data/"


# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_Mpc_Myr = 4.49e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_Mpc_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 

# Cluster

const _G = 1.0
const _Mtot = 1.0


# Host potential
const g_1 = log(2) - 1/2
const g_c = log(1+c) - c/(1+c)

const M_DH_200 =M_DH_200_Msun/(M_HU_in_Msun)
const Rs = Rs_kpc/(R_HU_in_kpc)
const Menc = M_DH_200*g_1/g_c # Enclosed mass within a sphere of radius Rs
const rho0_host = Menc/(4*pi*Rs^3 * g_1)

function psi_NFW(r::Float64)

    pref = 4*pi*_G*rho0_host*Rs^3
    return -pref/r * log(1 + r/Rs)
end


function plot_data()

    if (isfile(path_data*"snapshots/.DS_Store"))
        rm(path_data*"snapshots/.DS_Store")
    end

    listFile = readdir(path_data*"snapshots/";join=true)
    nsnap = length(listFile)

    tab_time = zeros(Float64, nsnap)

    for i=1:nsnap
        interm = split(split(listFile[i],"_")[end],".")
        index = parse(Float64, interm[1]*"."*interm[2])
        tab_time[i] = index
    end

    p = sortperm(tab_time)

    tab_energy = zeros(Float64, nsnap)

    for i=1:nsnap

        println("Progress = ", i/nsnap)
        namefile = listFile[p[i]]
        data = readdlm(namefile, header=false)
        Npart = length(data[:,1])
        mass = _Mtot/Npart

        # Kinetic energy 

        KE_threads = zeros(Float64, Threads.nthreads())

        Threads.@threads for k=1:Npart
            x, y, z, vx, vy, vz = data[k, :]
            tid = Threads.threadid()
            KE_threads[tid] += 0.5*mass*(vx^2 + vy^2 + vz^2)
        end

        KE = 0.0
        for tid=1:Threads.nthreads()
            KE += KE_threads[tid]
        end

         # Potential energy of cluster
        PE_threads = zeros(Float64, Threads.nthreads())

        Threads.@threads for ipair=1:Npart*Npart
            
            tid = Threads.threadid()

            # ipair-1 = i1-1 + (i2-1)*Npart
            i2 = div(ipair-1, Npart) + 1
            i1 = ipair - (i2-1)*Npart 

            if (i1 < i2)
                x1, y1, z1, _ = data[i1, :]
                x2, y2, z2, _ = data[i2, :]
                r12 = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)

                PE_threads[tid] += -_G*mass*mass/r12
            end

        end

        PE = 0.0
        for tid=1:Threads.nthreads()
            PE += PE_threads[tid]
        end

        # Add host potential 
        for k=1:Npart

            x, y, z, _ = data[k, :]
            r = sqrt(x^2 + y^2 + z^2)
            PE += mass * psi_NFW(r)

        end


        tab_energy[i] = PE + KE

    end

    tab_energy_frac = zeros(Float64, nsnap-1)
    ymin = Inf

    for i=2:nsnap

        tab_energy_frac[i-1] = abs(tab_energy[i]/tab_energy[1] - 1)
        if (tab_energy_frac[i-1] < ymin)
            ymin = tab_energy_frac[i-1]
        end

    end

    log10ymin = floor(Int64, log10(ymin))
    ymin = 10^(1.0*log10ymin)


    plt = plot(tab_time[p] .* T_HU_in_Myr ./10^3, [tab_energy], 
                    xlabel="Time [Gyr]", ylabel="Energy [HU]", 
                    framestyle=:box, labels=:false, 
                    # xlims=(-rmax, rmax), ylims=(-rmax,rmax), 
                    # aspect_ratio=1, size=(800,800), 
                    left_margin = 2Plots.mm, right_margin = 2Plots.mm, 
                    title="Energy conservation")

    # display(plt)
    # readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/IOM_plummer.pdf"
    namefile_png = path_data*"plot/IOM_plummer.png"
    savefig(plt, namefile_pdf)
    savefig(plt, namefile_png)



    # Fractional energy 

    plt2 = plot(tab_time[p][2:nsnap] .* T_HU_in_Myr ./10^3, [tab_energy_frac], 
                    xlabel="Time [Gyr]", ylabel="Fractional energy", 
                    framestyle=:box, labels=:false, 
                    ylims=(ymin,1.0),
                    yaxis=:log10,
                    # xlims=(-rmax, rmax), ylims=(-rmax,rmax), 
                    # aspect_ratio=1, size=(800,800), 
                    left_margin = 2Plots.mm, right_margin = 2Plots.mm, 
                    title="Energy conservation")

    # display(plt2)
    # readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/Fractional_IOM_plummer.pdf"
    namefile_png = path_data*"plot/Fractional_IOM_plummer.png"
    savefig(plt2, namefile_pdf)
    savefig(plt2, namefile_png)

    return nothing

end

@time plot_data()