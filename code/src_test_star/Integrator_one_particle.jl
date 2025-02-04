function integrate_one_leapfrog!(tab_stars::Array{Float64})

    # Integrate each star
    # N-body forces + Host potential

    
    # Leapfrog 
    # https://en.wikipedia.org/wiki/Leapfrog_integration#Algorithm

    # v_{k} -> v_{k+1/2} = v_{k} + a_{k}*dt/2
    # a_{k} = F(x_{k})
   

    x, y, z, vx, vy, vz = tab_stars[:]
    Fx_host, Fy_host, Fz_host = force_host(x, y, z) #force_NFW(x, y, z)

    ax = (Fx_host)
    ay = (Fy_host)
    az = (Fz_host)

    dvx = dt/2 * ax
    dvy = dt/2 * ay
    dvz = dt/2 * az

    tab_stars[4] = vx + dvx
    tab_stars[5] = vy + dvy
    tab_stars[6] = vz + dvz

    # x_{k} -> x_{k+1} = x_{k} + v_{k+1/2}*dt/2

    x, y, z, vx, vy, vz = tab_stars[:]

    dx = dt * vx
    dy = dt * vy 
    dz = dt * vz

    tab_stars[1] = x + dx
    tab_stars[2] = y + dy
    tab_stars[3] = z + dz

    # v_{k+1/2} -> v_{k+1} = v_{k+1/2} + a_{k+1}*dt/2
    # a_{k+1} = F(x_{k+1})

    x, y, z, vx, vy, vz = tab_stars[:]
    Fx_host, Fy_host, Fz_host = force_host(x, y, z) #force_NFW(x, y, z)

    ax = (Fx_host)
    ay = (Fy_host)
    az = (Fz_host)

    dvx = dt/2 * ax
    dvy = dt/2 * ay
    dvz = dt/2 * az

    tab_stars[4] = vx + dvx
    tab_stars[5] = vy + dvy
    tab_stars[6] = vz + dvz


    return nothing

end