
using DelimitedFiles
using Plots 
using LaTeXStrings
using Plots.PlotMeasures
using ArgParse
using NearestNeighbors
using LaTeXStrings

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
    default = 63875410815852
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


function get_data()

    listFile = readdir(path_data*"snapshots_"*srun*"/";join=true)
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

    tab_IOM = zeros(Float64, nsnap, 8) # time, E_tot, Lz, nb_unbound

    for isnap=1:nsnap

        namefile = sortedFiles[isnap]
        time = tabtsort[isnap]

        data_stars = readdlm(namefile) # x, y, z, vx, vy, vz, Uint, Uc

        tab_vb_t = zeros(Float64, Threads.nthreads(), 3)
        K_t =  zeros(Float64, Threads.nthreads())
        Uc_t =  zeros(Float64, Threads.nthreads())
        Uh_t =  zeros(Float64, Threads.nthreads())
        L_t =  zeros(Float64, Threads.nthreads(), 3)

        Threads.@threads for i=1:Npart 

            tid = Threads.threadid()
            x, y, z, vx, vy, vz, Uint, Uc = data_stars[i, :]
            

            v2 = vx^2 + vy^2 + vz^2

            # Barycenter 
            tab_vb_t[tid, 1] += vx
            tab_vb_t[tid, 2] += vy
            tab_vb_t[tid, 3] += vz

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
       
        Vbx = 0.0
        Vby = 0.0
        Vbz = 0.0

        K = 0.0
        U = 0.0
        L = zeros(Float64, 3)


        for tid=1:Threads.nthreads()
            Vbx += tab_vb_t[tid, 1]/Npart
            Vby += tab_vb_t[tid, 2]/Npart
            Vbz += tab_vb_t[tid, 3]/Npart

            K += K_t[tid]
            U += Uh_t[tid] + Uc_t[tid]
            L[1] += L_t[tid, 1]
            L[2] += L_t[tid, 2]
            L[3] += L_t[tid, 3]
    
        end

        Etot = K + U

        # Compute number of unbound stars
        n_unbound_t = zeros(Float64, Threads.nthreads())

        Threads.@threads for i=1:Npart

            tid = Threads.threadid()
    
            # Unbound particles 
            x, y, z, vx, vy, vz, Uint = data_stars[i, :]
            vc_x = vx - Vbx
            vc_y = vy - Vby
            vc_z = vz - Vbz
    
            Ec = 0.5 * mass * (vc_x^2 + vc_y^2 + vc_z^2) + Uint
    
            if (Ec >= 0.0)
                n_unbound_t[tid] += 1.0
            end
    
        end

        n_unbound = 0.0

        for tid=1:Threads.nthreads()
            n_unbound += n_unbound_t[tid]
        end
    
        tab_IOM[isnap,1] = time
        tab_IOM[isnap,2] = Etot
        tab_IOM[isnap,3] = L[3]
        tab_IOM[isnap,4] = n_unbound

    end

    return tab_IOM

end

function plot_data!()

    tab_IOM = get_data()
    datat = tab_IOM[:, 1]  .* T_HU_in_Myr
    dataE = tab_IOM[:, 2]
    dataLz = tab_IOM[:, 3]
    dataUnbound = tab_IOM[:, 4]

    # Energy
    n = length(dataE)
    dataFracE = abs.(1.0 .- dataE[2:n] ./ dataE[1])

    plt = plot(datat, [dataFracE], 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel="Fractional energy", 
        yaxis=:log10,
        yticks=10.0 .^ (-15:1:0),
        xlims=(0, datat[n-1]),
        frame=:box)

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/IOM_frac_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

    # Unbound particles

    dataUnbound = data[:,9]
    n = length(dataUnbound)

    plt = plot(datat, [dataUnbound ./ Npart .* 100], 
        labels=:false, 
        xlabel="Time [Myr]", 
        ylabel="Fraction of unbound stars [%]", 
        xlims=(0, datat[n-1]),
        frame=:box)

    display(plt)
    readline()

    mkpath(path_data*"plot/")
    namefile_pdf = path_data*"plot/count_unbound_cluster_"*srun*".pdf"
    savefig(plt, namefile_pdf)

end 


plot_data!()