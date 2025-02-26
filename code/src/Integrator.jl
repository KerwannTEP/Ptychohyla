
function update_tab_acc_Uint!(tab_stars::Array{Float64}, tab_acc::Array{Float64}, tab_Uint::Array{Float64}, tab_Uc::Array{Float64})


    # Add the contribution from the host 
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz, m = tab_stars[i, :]
        ax_internal, ay_internal, az_internal, U_internal = acc_U_internal(i, tab_stars)

        ax = ax_internal
        ay = ay_internal
        az = az_internal

        if (HAS_HOST)
            ax_host, ay_host, az_host = acc_host(x, y, z) 

            ax += ax_host
            ay += ay_host
            az += az_host
        end

        tab_acc[i, 1] = ax
        tab_acc[i, 2] = ay
        tab_acc[i, 3] = az
        tab_Uint[i] = U_internal

        if (HAS_HOST)
            r = sqrt(x^2 + y^2 + z^2)
            R = sqrt(x^2 + y^2)
            psi_xyz = psi_halo(r) + psi_disk(R, z) + psi_bulge(r) 

            tab_Uc[i] = m * psi_xyz
        end

        

    end

end


function integrate_stars_leapfrog!(index::Int64, time::Float64, tab_stars::Array{Float64}, tab_acc::Array{Float64}, tab_Uint::Array{Float64}, tab_Uc::Array{Float64}, first_timestep::Bool)

    # Integrate each star
    # N-body force + Host potential

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

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

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

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

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

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

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








# 4th order Yoshida integrator
# https://en.wikipedia.org/wiki/Leapfrog_integration#4th_order_Yoshida_integrator


function integrate_stars_yoshida!(index::Int64, time::Float64, tab_stars::Array{Float64}, tab_acc::Array{Float64}, tab_Uint::Array{Float64}, tab_Uc::Array{Float64}, first_timestep::Bool)

    w0 = - cbrt(2.0)/(2.0 - cbrt(2.0))
    w1 = 1.0/(2.0 - cbrt(2.0))
    c1 = 0.5*w1 
    c4 = c1 
    c2 = 0.5*(w0 + w1)
    c3 = c2
    d1 = w1 
    d3 = d1 
    d2 = w0

    # Integrate each star
    # N-body force + Host potential

    


    # Write snapshot data
    # If first timestep, then Uint are computed just above.
    # If not, then position, velocities and Uint were updated in the previous timestep
    if (index % N_dt == 0)
        if !(index == 0 && (RESTART)) # If not the IC of the restart (already saved in previous run)

            # Compute integral of motions
            update_tab_acc_Uint!(tab_stars, tab_acc, tab_Uint, tab_Uc)

            write_data!(time, tab_stars, tab_Uint, tab_Uc)
        end  
    end

    # Drift 1
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

        dx = c1 * vx * dt
        dy = c1 * vy * dt
        dz = c1 * vz * dt

        tab_stars[i, 1] = x + dx
        tab_stars[i, 2] = y + dy
        tab_stars[i, 3] = z + dz

    end

    # Kick 1
    update_tab_acc_Uint!(tab_stars, tab_acc, tab_Uint, tab_Uc)
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

        ax = tab_acc[i, 1]
        ay = tab_acc[i, 2]
        az = tab_acc[i, 3]

        dvx = d1 * ax * dt 
        dvy = d1 * ay * dt 
        dvz = d1 * az * dt 

        tab_stars[i, 4] = vx + dvx
        tab_stars[i, 5] = vy + dvy
        tab_stars[i, 6] = vz + dvz

    end

    # Drift 2
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

        dx = c2 * vx * dt
        dy = c2 * vy * dt
        dz = c2 * vz * dt

        tab_stars[i, 1] = x + dx
        tab_stars[i, 2] = y + dy
        tab_stars[i, 3] = z + dz

    end

    # Kick 2
    update_tab_acc_Uint!(tab_stars, tab_acc, tab_Uint, tab_Uc)
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

        ax = tab_acc[i, 1]
        ay = tab_acc[i, 2]
        az = tab_acc[i, 3]

        dvx = d2 * ax * dt 
        dvy = d2 * ay * dt 
        dvz = d2 * az * dt 

        tab_stars[i, 4] = vx + dvx
        tab_stars[i, 5] = vy + dvy
        tab_stars[i, 6] = vz + dvz

    end

    # Drift 3
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

        dx = c3 * vx * dt
        dy = c3 * vy * dt
        dz = c3 * vz * dt

        tab_stars[i, 1] = x + dx
        tab_stars[i, 2] = y + dy
        tab_stars[i, 3] = z + dz

    end

    # Kick 3
    update_tab_acc_Uint!(tab_stars, tab_acc, tab_Uint, tab_Uc)
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

        ax = tab_acc[i, 1]
        ay = tab_acc[i, 2]
        az = tab_acc[i, 3]

        dvx = d3 * ax * dt 
        dvy = d3 * ay * dt 
        dvz = d3 * az * dt 

        tab_stars[i, 4] = vx + dvx
        tab_stars[i, 5] = vy + dvy
        tab_stars[i, 6] = vz + dvz

    end

    # Drift 4
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz, m = tab_stars[i, :]

        dx = c4 * vx * dt
        dy = c4 * vy * dt
        dz = c4 * vz * dt

        tab_stars[i, 1] = x + dx
        tab_stars[i, 2] = y + dy
        tab_stars[i, 3] = z + dz

    end

end
