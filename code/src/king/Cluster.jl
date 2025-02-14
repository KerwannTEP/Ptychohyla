using DelimitedFiles
using CSV 
using DataFrames

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


# force on star k from all other stars
# tab_stars : (k, 6) : star k : (x,y,z,vx,vy,vz) : Npart stars
function acc_internal(k::Int64, tab_stars::Array{Float64})

    ax = 0.0
    ay = 0.0
    az = 0.0

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

            ax += intensity * xik
            ay += intensity * yik
            az += intensity * zik

        end

    end

    return ax, ay, az

end


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
            intensityU = - _G*mass/rik

            ax += intensity * xik
            ay += intensity * yik
            az += intensity * zik
            U += intensityU

        end

    end



    return ax, ay, az, U

end

# force on star k from all other stars
# tab_stars : (k, 6) : star k : (x,y,z,vx,vy,vz) : Npart stars
function force_internal(k::Int64, tab_stars::Array{Float64})

    ax, ay, az = acc_internal(k, tab_stars)

    Fx = mass * ax
    Fy = mass * ay
    Fz = mass * az

    return Fx, Fy, Fz

end
