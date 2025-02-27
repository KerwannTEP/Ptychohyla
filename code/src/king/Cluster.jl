using DelimitedFiles
using CSV 
using DataFrames

############################################################################################################################################
#  Initialization
############################################################################################################################################

function initialize_stars!(tab_stars::Array{Float64})

    # Load King sphere in Henon units
    namefile = path_dir*"data/IC/chc_king_ics_n_"*string(Npart)*".csv"
    df = CSV.read(namefile, DataFrame, delim=',', header=false)


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

function acc_jerk_U_internal(k::Int64, tab_stars::Array{Float64})

    ax = 0.0
    ay = 0.0
    az = 0.0
    U = 0.0

    jx = 0.0
    jy = 0.0
    jz = 0.0

    xk = tab_stars[k,1]
    yk = tab_stars[k,2]
    zk = tab_stars[k,3]

    vxk = tab_stars[k,4]
    vyk = tab_stars[k,5]
    vzk = tab_stars[k,6]

    for i=1:Npart
        if (i != k)

            xi = tab_stars[i,1]
            yi = tab_stars[i,2]
            zi = tab_stars[i,3]

            vxk = tab_stars[i,4]
            vyk = tab_stars[i,5]
            vzk = tab_stars[i,6]

            # Vector ri - rk
            xik = xi - xk 
            yik = yi - yk 
            zik = zi - zk 

            # Vector vi - vk
            vxik = vxi - vxk 
            vyik = vyi - vyk 
            vzik = vzi - vzk 

            rik = sqrt(xik^2 + yik^2 + zik^2 + eps^2)

            intensity = _G*mass/rik^3 
            intensityU = - _G*mass*mass/rik

            ax += intensity * xik
            ay += intensity * yik
            az += intensity * zik
            U += intensityU

            vik_dot_rik = vxik*xik + vyik*yik + vzik*zik
            jx += _G*mass*(vxik/rik^3 - 3.0*(vik_dot_rik)*xk/rik^5)
            jy += _G*mass*(vyik/rik^3 - 3.0*(vik_dot_rik)*yk/rik^5)
            jz += _G*mass*(vzik/rik^3 - 3.0*(vik_dot_rik)*zk/rik^5)

        end

    end

    return ax, ay, az, jx, jy, jz, U

end

function snap_internal(k::Int64, tab_stars::Array{Float64}, tab_acc::Array{Float64}, tab_jerk::Array{Float64})

    sx = 0.0
    sy = 0.0
    sz = 0.0

    xk = tab_stars[k,1]
    yk = tab_stars[k,2]
    zk = tab_stars[k,3]

    vxk = tab_stars[k,4]
    vyk = tab_stars[k,5]
    vzk = tab_stars[k,6]

    axk = tab_acc[k,1]
    ayk = tab_acc[k,2]
    azk = tab_acc[k,3]

    jxk = tab_jerk[k,1]
    jyk = tab_jerk[k,2]
    jzk = tab_jerk[k,3]

    for i=1:Npart
        if (i != k)

            xi = tab_stars[i,1]
            yi = tab_stars[i,2]
            zi = tab_stars[i,3]

            vxk = tab_stars[i,4]
            vyk = tab_stars[i,5]
            vzk = tab_stars[i,6]

            axi = tab_acc[i,1]
            ayi = tab_acc[i,2]
            azi = tab_acc[i,3]

            jxi = tab_jerk[i,1]
            jyi = tab_jerk[i,2]
            jzi = tab_jerk[i,3]

            # Vector ri - rk
            xik = xi - xk 
            yik = yi - yk 
            zik = zi - zk 

            # Vector vi - vk
            vxik = vxi - vxk 
            vyik = vyi - vyk 
            vzik = vzi - vzk 

            # Vector ai - ak
            axik = axi - axk 
            ayik = ayi - ayk 
            azik = azi - azk 

            # Vector ji - jk
            jxik = jxi - jxk 
            jyik = jyi - jyk 
            jzik = jzi - jzk 

            rik = sqrt(xik^2 + yik^2 + zik^2 + eps^2)

            rik_dot_vik =  xik * vxik  +  yik * vyik  +  zik * vzik
            rik_dot_aik =  xik * axik  +  yik * ayik  +  zik * azik
            vik_dot_jik = vxik * jxik  + vyik * jyik  + vzik * jzik

            sx += _G*mass*(axik/rik^3 - 3.0*rik_dot_aik*xik/rik^5 - 3.0*vik_dot_jik*xik/rik^5 - 3.0*rik_dot_vik*jxik/rik^5 + 15.0*rik_dot_vik*rik_dot_aik*xik/rik^7)
            sy += _G*mass*(ayik/rik^3 - 3.0*rik_dot_aik*yik/rik^5 - 3.0*vik_dot_jik*yik/rik^5 - 3.0*rik_dot_vik*jyik/rik^5 + 15.0*rik_dot_vik*rik_dot_aik*yik/rik^7)
            sz += _G*mass*(azik/rik^3 - 3.0*rik_dot_aik*zik/rik^5 - 3.0*vik_dot_jik*zik/rik^5 - 3.0*rik_dot_vik*jzik/rik^5 + 15.0*rik_dot_vik*rik_dot_aik*zik/rik^7)

        end

    end

    return sx, sy, sz

end