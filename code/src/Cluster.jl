using DelimitedFiles

function initialize_stars!(tab_stars::Array{Float64})

    # namefile = path_dir*"code/IC_generator/data_IC/plummer_sphere_n_10_q_0.0.txt"
    # namefile = path_dir*"code/IC_generator/data_IC/plummer_sphere_n_100_q_0.0.txt"
    namefile = path_dir*"code/IC_generator/data_IC/plummer_sphere_n_1000_q_0.0.txt"
    data = readdlm(namefile, header=false)

    vcirc = circular_velocity(d_host)
    println("Vc = ", vcirc)


    Threads.@threads for i=1:Npart 

        tab_stars[i, 1] = d_host + data[i, 1] # Cluster on the x>0 x-axis
        tab_stars[i, 2] = data[i, 2]
        tab_stars[i, 3] = data[i, 3]
        tab_stars[i, 4] = data[i, 4]
        tab_stars[i, 5] = vcirc + data[i, 5] # Circular velocity: cluster goes in the y-direction
        tab_stars[i, 6] = data[i, 6]

    end


end


# force on star k from all other stars
# tab_stars : (k, 6) : star k : (x,y,z,vx,vy,vz) : Npart stars
function force_internal(k::Int64, tab_stars::Array{Float64})

    Fx = 0.0
    Fy = 0.0
    Fz = 0.0

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

            intensity = _G*mass*mass/rik^3 

            Fx += intensity * xik
            Fy += intensity * yik
            Fz += intensity * zik

        end

    end

    return Fx, Fy, Fz

end

