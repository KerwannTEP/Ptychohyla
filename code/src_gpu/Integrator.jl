function update_tab_acc_Uint!(tab_stars::Array{Float64}, tab_acc::Array{Float64}, tab_Uint::Array{Float64}, tab_Uc::Array{Float64})

    # Compute contribution from the cluster self-gravity (GPU acceleration)
    tab_acc_Uint_int_gpu!(tab_acc, tab_Uint, tab_stars[:, 1:3])

    # Add the contribution from the host 
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz = tab_stars[i, :]
        ax_host, ay_host, az_host = acc_host(x, y, z) 

        tab_acc[i, 1] += ax_host
        tab_acc[i, 2] += ay_host
        tab_acc[i, 3] += az_host

        r = sqrt(x^2 + y^2 + z^2)
        R = sqrt(x^2 + y^2)
        psi_xyz = psi_halo(r) + psi_disk(R, z) + psi_bulge(r) 

        tab_Uc[i] = mass * psi_xyz


    end

end

function integrate_stars_leapfrog!(index::Int64, time::Float64, tab_stars::Array{Float64}, tab_acc::Array{Float64}, tab_Uint::Array{Float64}, tab_Uc::Array{Float64}, first_timestep::Bool)

    # Integrate each star
    # N-body force (GPU) + Host potential

    # Leapfrog 
    # https://en.wikipedia.org/wiki/Leapfrog_integration#Algorithm

    if (first_timestep)
        update_tab_acc_Uint!(tab_stars, tab_acc, tab_Uint, tab_Uc)
    end


    # Write snapshot data
    # If first timestep, then Uint are computed just above.
    # If not, then position, velocities and Uint were updated in the previous timestep
    if (index % N_dt == 0)
        if !(index == 0 && (RESTART)) # If not the IC of the restart (already saved in previous run)
            write_data!(time, tab_stars, tab_Uint, tab_Uc)
        end  
    end


    # v_{k} -> v_{k+1/2} = v_{k} + a_{k}*dt/2
    # a_{k} = F(x_{k})
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz = tab_stars[i, :]

        ax = tab_acc[i, 1]
        ay = tab_acc[i, 2]
        az = tab_acc[i, 3]

        dvx = dt/2 * ax
        dvy = dt/2 * ay
        dvz = dt/2 * az

        tab_stars[i, 4] = vx + dvx
        tab_stars[i, 5] = vy + dvy
        tab_stars[i, 6] = vz + dvz

    end

    # x_{k} -> x_{k+1} = x_{k} + v_{k+1/2}*dt/2
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz = tab_stars[i, :]

        dx = dt * vx
        dy = dt * vy 
        dz = dt * vz

        tab_stars[i, 1] = x + dx
        tab_stars[i, 2] = y + dy
        tab_stars[i, 3] = z + dz

    end


    update_tab_acc_Uint!(tab_stars, tab_acc, tab_Uint, tab_Uc)

    # v_{k+1/2} -> v_{k+1} = v_{k+1/2} + a_{k+1}*dt/2
    # a_{k+1} = F(x_{k+1})
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz = tab_stars[i, :]

        ax = tab_acc[i, 1]
        ay = tab_acc[i, 2]
        az = tab_acc[i, 3]

        dvx = dt/2 * ax
        dvy = dt/2 * ay
        dvz = dt/2 * az

        tab_stars[i, 4] = vx + dvx
        tab_stars[i, 5] = vy + dvy
        tab_stars[i, 6] = vz + dvz

    end

    return nothing

end
