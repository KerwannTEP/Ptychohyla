using DelimitedFiles

function initialize_stars!(tab_stars::Array{Float64})


    namefile = path_dir*"data/IC/plummer_sphere_n_"*string(Npart)*"_q_"*string(q)*".txt"
    data = readdlm(namefile, header=false)

    vcirc = circular_velocity(d_host)
    println("Vc = ", vcirc)
    println("Vc [km/s] = ", vcirc * V_HU_in_km_s)

    Threads.@threads for i=1:Npart 

        tab_stars[i, 1] = d_host + data[i, 1] # Cluster on the x>0 x-axis
        tab_stars[i, 2] = data[i, 2]
        tab_stars[i, 3] = data[i, 3]
        tab_stars[i, 4] = data[i, 4]
        tab_stars[i, 5] = vcirc + data[i, 5] # Circular velocity: cluster goes in the y-direction
        tab_stars[i, 6] = data[i, 6]

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

    data_stars = readdlm(namefile) # x, y, z, vx, vy, vz, Uint, Uc
    # Array of size (Npart, 8), in Henon units

    Threads.@threads for i=1:Npart 

        tab_stars[i, 1] = data_stars[i, 1]
        tab_stars[i, 2] = data_stars[i, 2]
        tab_stars[i, 3] = data_stars[i, 3]
        tab_stars[i, 4] = data_stars[i, 4]
        tab_stars[i, 5] = data_stars[i, 5]
        tab_stars[i, 6] = data_stars[i, 6]

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

    for i=1:Npart
        if (i != k)

            xi = tab_stars[i,1]
            yi = tab_stars[i,2]
            zi = tab_stars[i,3]

            # Vector ri - rk
            xik = xi - xk 
            yik = yi - yk 
            zik = zi - zk 

            rik = sqrt(xik^2 + yik^2 + zik^2 + eps^2)

            intensity = _G*mass/rik^3 
            intensityU = - _G*mass*mass/rik

            ax += intensity * xik
            ay += intensity * yik
            az += intensity * zik
            U += intensityU

        end

    end



    return ax, ay, az, U

end
