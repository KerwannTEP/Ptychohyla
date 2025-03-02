using DelimitedFiles
using CSV 
using DataFrames

############################################################################################################################################
#  Initialization
############################################################################################################################################

function initialize_stars!(tab_stars::Array{Float64})

    # Load King sphere in Henon units
    if (!HAS_MULTI_MASS) 
        # Single mass
        namefile = path_dir*"data/IC/chc_king_ics_n_"*string(Npart)*".csv"
    else 
        # Multi-mass (half of mass m1, half of mass m2=2*m1)
        namefile = path_dir*"data/IC/chc_king_ics_n_"*string(Npart)*"_multi_mass_1_2.csv"
    end
    
    df = CSV.read(namefile, DataFrame, delim=',', header=false)

    datam = df[:, 3]
    datax = df[:, 6]
    datay = df[:, 7]
    dataz = df[:, 8]
    datavx = df[:, 9]
    datavy = df[:, 10]
    datavz = df[:, 11]


    vcirc = circular_velocity(d_host)
    println("Vc [HU] = ", vcirc)
    println("Vc [km/s] = ", vcirc * V_HU_in_km_s)

    Threads.@threads for i=1:Npart 

        tab_stars[i, 1] = d_host + datax[i] # Cluster on the x>0 x-axis
        tab_stars[i, 2] = datay[i]
        tab_stars[i, 3] = dataz[i]
        tab_stars[i, 4] = datavx[i]
        tab_stars[i, 5] = vcirc + datavy[i] # Circular velocity: cluster goes in the y-direction
        tab_stars[i, 6] = datavz[i]
        tab_stars[i, 7] = datam[i]

    end


end

function initialize_stars_restart!(tab_stars::Array{Float64})

    # Load last saved snapshot 
    # Fill data
    # Return restart time

    listFiles = readdir(folder_output*"snapshots_"*srun*"/"; join=true)
    nsnap = length(listFiles)
    tabt = zeros(Float64, nsnap)

    for i=1:nsnap 
        interm = split(split(listFiles[i],"_")[end],".")
        interm = interm[1]*"."*interm[2]
        time = parse(Float64, interm)

        tabt[i] = time 
    end

    p = sortperm(tabt)

    namefile = listFiles[p[nsnap]]
    time = tabt[p[nsnap]]

    data_stars = readdlm(namefile) # x, y, z, vx, vy, vz, m, Uint, Uc
    # Array of size (Npart, 8), in Henon units

    Threads.@threads for i=1:Npart 

        tab_stars[i, 1] = data_stars[i, 1]
        tab_stars[i, 2] = data_stars[i, 2]
        tab_stars[i, 3] = data_stars[i, 3]
        tab_stars[i, 4] = data_stars[i, 4]
        tab_stars[i, 5] = data_stars[i, 5]
        tab_stars[i, 6] = data_stars[i, 6]
        tab_stars[i, 7] = data_stars[i, 7]

    end

    return time

end

############################################################################################################################################
# Accelerations and interaction energies
############################################################################################################################################



function acc_U_internal(k::Int64, tab_stars::Array{Float64})

    ax = 0.0
    ay = 0.0
    az = 0.0
    U = 0.0

    xk = tab_stars[k,1]
    yk = tab_stars[k,2]
    zk = tab_stars[k,3]
    mk = tab_stars[k,7]

    for i=1:Npart
        if (i != k)

            xi = tab_stars[i,1]
            yi = tab_stars[i,2]
            zi = tab_stars[i,3]
            mi = tab_stars[k,7]

            # Vector ri - rk
            xik = xi - xk 
            yik = yi - yk 
            zik = zi - zk 

            rik = sqrt(xik^2 + yik^2 + zik^2 + eps^2)

            intensity = _G*mi/rik^3 
            intensityU = - _G*mk*mi/rik

            ax += intensity * xik
            ay += intensity * yik
            az += intensity * zik
            U += intensityU

        end

    end



    return ax, ay, az, U

end
