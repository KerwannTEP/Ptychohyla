
function integrate_stars_euler!(tab_stars::Array{Float64})

    # Integrate each star
    # N-body forces + Host potential

    
    # Euler integration
    # Fix tab_stars_intermediary
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz = tab_stars[i, :]
        Fx_internal, Fy_internal, Fz_internal = force_internal(i, tab_stars)
        Fx_host, Fy_host, Fz_host = force_host(x, y, z) #force_NFW(x, y, z)

        ax = (Fx_internal + Fx_host)/mass
        ay = (Fy_internal + Fy_host)/mass
        az = (Fz_internal + Fz_host)/mass

        
        dx = dt * vx
        dy = dt * vy 
        dz = dt * vz
        dvx = dt * ax
        dvy = dt * ay
        dvz = dt * az

        tab_stars[i, 1] = x + dx
        tab_stars[i, 2] = y + dy
        tab_stars[i, 3] = z + dz
        tab_stars[i, 4] = vx + dvx
        tab_stars[i, 5] = vy + dvy
        tab_stars[i, 6] = vz + dvz

    end


    return nothing

end



function integrate_stars_leapfrog!(tab_stars::Array{Float64})

    # Integrate each star
    # N-body forces + Host potential

    # Fix tab_stars_intermediary
    tab_stars_init = tab_stars
    
    # Leapfrog 
    # https://en.wikipedia.org/wiki/Leapfrog_integration#Algorithm

    # v_{k} -> v_{k+1/2} = v_{k} + a_{k}*dt/2
    # a_{k} = F(x_{k})
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz = tab_stars_init[i, :]
        Fx_internal, Fy_internal, Fz_internal = force_internal(i, tab_stars_init)
        Fx_host, Fy_host, Fz_host = force_host(x, y, z) #force_NFW(x, y, z)

        ax = (Fx_internal + Fx_host)/mass
        ay = (Fy_internal + Fy_host)/mass
        az = (Fz_internal + Fz_host)/mass


        dvx = dt/2 * ax
        dvy = dt/2 * ay
        dvz = dt/2 * az

        tab_stars[i, 4] = vx + dvx
        tab_stars[i, 5] = vy + dvy
        tab_stars[i, 6] = vz + dvz

    end

    tab_stars_init = tab_stars

    # x_{k} -> x_{k+1} = x_{k} + v_{k+1/2}*dt/2
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz = tab_stars_init[i, :]

        dx = dt * vx
        dy = dt * vy 
        dz = dt * vz

        tab_stars[i, 1] = x + dx
        tab_stars[i, 2] = y + dy
        tab_stars[i, 3] = z + dz

    end

    tab_stars_init = tab_stars


    # v_{k+1/2} -> v_{k+1} = v_{k+1/2} + a_{k+1}*dt/2
    # a_{k+1} = F(x_{k+1})
    Threads.@threads for i=1:Npart 

        x, y, z, vx, vy, vz = tab_stars_init[i, :]
        Fx_internal, Fy_internal, Fz_internal = force_internal(i, tab_stars_init)
        Fx_host, Fy_host, Fz_host = force_host(x, y, z) #force_NFW(x, y, z)

        ax = (Fx_internal + Fx_host)/mass
        ay = (Fy_internal + Fy_host)/mass
        az = (Fz_internal + Fz_host)/mass

        dvx = dt/2 * ax
        dvy = dt/2 * ay
        dvz = dt/2 * az

        tab_stars[i, 4] = vx + dvx
        tab_stars[i, 5] = vy + dvy
        tab_stars[i, 6] = vz + dvz

    end

    
    

    return nothing

end