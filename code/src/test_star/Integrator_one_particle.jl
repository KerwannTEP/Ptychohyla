function integrate_one_leapfrog!(tab_stars::Array{Float64})

    # Integrate each star
    # N-body forces + Host potential

    
    # Leapfrog 
    # https://en.wikipedia.org/wiki/Leapfrog_integration#Algorithm

    # v_{k} -> v_{k+1/2} = v_{k} + a_{k}*dt/2
    # a_{k} = F(x_{k})
   

    x, y, z, vx, vy, vz = tab_stars[:]
    ax, ay, az = acc_host(x, y, z) 

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
    ax, ay, az = acc_host(x, y, z) 


    dvx = dt/2 * ax
    dvy = dt/2 * ay
    dvz = dt/2 * az

    tab_stars[4] = vx + dvx
    tab_stars[5] = vy + dvy
    tab_stars[6] = vz + dvz


    return nothing

end



# 4th order Yoshida integrator
# https://en.wikipedia.org/wiki/Leapfrog_integration#4th_order_Yoshida_integrator

function integrate_one_yoshida!(tab_stars::Array{Float64})

    w0 = - cbrt(2.0)/(2.0 - cbrt(2.0))
    w1 = 1.0/(2.0 - cbrt(2.0))
    c1 = 0.5*w1 
    c4 = c1 
    c2 = 0.5*(w0 + w1)
    c3 = c2
    d1 = w1 
    d3 = d1 
    d2 = w0

    # Drift 1
    x, y, z, vx, vy, vz = tab_stars[:]

    dx = c1 * vx * dt
    dy = c1 * vy * dt
    dz = c1 * vz * dt

    tab_stars[1] = x + dx
    tab_stars[2] = y + dy
    tab_stars[3] = z + dz


    # Kick 1
    x, y, z, vx, vy, vz = tab_stars[:]
    ax, ay, az = acc_host(x, y, z) 

    dvx = d1 * ax * dt 
    dvy = d1 * ay * dt 
    dvz = d1 * az * dt 

    tab_stars[4] = vx + dvx
    tab_stars[5] = vy + dvy
    tab_stars[6] = vz + dvz


    # Drift 2
    x, y, z, vx, vy, vz = tab_stars[:]

    dx = c2 * vx * dt
    dy = c2 * vy * dt
    dz = c2 * vz * dt

    tab_stars[1] = x + dx
    tab_stars[2] = y + dy
    tab_stars[3] = z + dz


    # Kick 2
    x, y, z, vx, vy, vz = tab_stars[:]
    ax, ay, az = acc_host(x, y, z) 

    dvx = d2 * ax * dt 
    dvy = d2 * ay * dt 
    dvz = d2 * az * dt 

    tab_stars[4] = vx + dvx
    tab_stars[5] = vy + dvy
    tab_stars[6] = vz + dvz


    # Drift 3
    x, y, z, vx, vy, vz = tab_stars[:]

    dx = c3 * vx * dt
    dy = c3 * vy * dt
    dz = c3 * vz * dt

    tab_stars[1] = x + dx
    tab_stars[2] = y + dy
    tab_stars[3] = z + dz


    # Kick 3
    x, y, z, vx, vy, vz = tab_stars[:]
    ax, ay, az = acc_host(x, y, z) 

    dvx = d3 * ax * dt 
    dvy = d3 * ay * dt 
    dvz = d3 * az * dt 

    tab_stars[4] = vx + dvx
    tab_stars[5] = vy + dvy
    tab_stars[6] = vz + dvz


    # Drift 4
    x, y, z, vx, vy, vz = tab_stars[:]

    dx = c4 * vx * dt
    dy = c4 * vy * dt
    dz = c4 * vz * dt

    tab_stars[1] = x + dx
    tab_stars[2] = y + dy
    tab_stars[3] = z + dz

    
    return nothing

end