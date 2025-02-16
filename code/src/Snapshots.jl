using DelimitedFiles


function write_data!(time::Float64, tab_stars::Array{Float64}, tab_Uint::Array{Float64}, tab_IOM::Array{Float64}, tab_bary::Array{Float64})
    
    n_digits = floor(Int64,-log10(dt))+1
    namefile = path_dir*"data/snapshots_"*srun*"/time_"*string(round(time, digits=n_digits))*".txt"
    writedlm(namefile, [tab_stars tab_Uint])

    # https://discourse.julialang.org/t/adding-to-existing-txt-file-created-by-writedlm/6907/2
    namefile_iom = path_dir*"data/iom_snapshots_"*srun*".txt"
    io = open(namefile_iom, "a")
    writedlm(io, transpose([round(time, digits=n_digits); tab_IOM]))
    close(io)

    namefile_bary = path_dir*"data/bary_snapshots_"*srun*".txt"
    io2 = open(namefile_bary, "a")
    writedlm(io2, transpose([round(time, digits=n_digits); tab_bary]))
    close(io2)

    return nothing

end

# TODO: Density center to compute ?
function compute_bary!(tab_stars::Array{Float64}, tab_bary::Array{Float64})

    xb_t = zeros(Float64, Threads.nthreads())
    yb_t = zeros(Float64, Threads.nthreads())
    zb_t = zeros(Float64, Threads.nthreads())

    Threads.@threads for i=1:Npart

        tid = Threads.threadid()

        x, y, z, vx, vy, vz = tab_stars[i, :]

        xb_t[tid] += x
        yb_t[tid] += y
        zb_t[tid] += z

    end

    tab_bary[1] = 0.0
    tab_bary[2] = 0.0
    tab_bary[3] = 0.0

    for tid=1:Threads.nthreads()

        tab_bary[1] += xb_t[tid]/Npart 
        tab_bary[2] += yb_t[tid]/Npart 
        tab_bary[3] += zb_t[tid]/Npart 

    end

    

end


function compute_IOM!(tab_stars::Array{Float64}, tab_IOM::Array{Float64}, tab_Uint::Array{Float64}, first_timestep::Bool=false)

    # tab_IOM = zeros(Float64, 7) # K, U, Etot, Lx, Ly, Lz, L

    # Compute velocity barycenter 
    
    tab_vb_t = zeros(Float64, Threads.nthreads(), 3)

    Threads.@threads for i=1:Npart 

        tid = Threads.threadid()
        x, y, z, vx, vy, vz = tab_stars[i, :]
        tab_vb_t[tid, 1] += vx
        tab_vb_t[tid, 2] += vy
        tab_vb_t[tid, 3] += vz

    end

    Vbx = 0.0
    Vby = 0.0
    Vbz = 0.0

    for tid=1:Threads.nthreads()
        Vbx += tab_vb_t[tid, 1]/Npart
        Vby += tab_vb_t[tid, 2]/Npart
        Vbz += tab_vb_t[tid, 3]/Npart
    end

    K_t = zeros(Float64, Threads.nthreads()) # Kinetic energy
    Uh_t = zeros(Float64, Threads.nthreads()) # Potential energy from host
    Uc_t = zeros(Float64, Threads.nthreads()) # Potential energy from cluster (self-interaction)
    L_t = zeros(Float64, Threads.nthreads(), 3) # Components of the angular momentum
    n_unbound_t = zeros(Float64, Threads.nthreads())

    Threads.@threads for i=1:Npart

        tid = Threads.threadid()

        x, y, z, vx, vy, vz = tab_stars[i, :]
        r = sqrt(x^2 + y^2 + z^2)
        R = sqrt(x^2 + y^2)
        v2 = vx^2 + vy^2 + vz^2
        psi_xyz = psi_halo(r) + psi_disk(R, z) + psi_bulge(r) 

       

        # Kinetic energy 
        K_t[tid] += 0.5 * mass * v2 

        # Host potential energy
        Uh_t[tid] += mass * psi_xyz

        
        # Angular momenta 
        Lx = y*vz - z*vy 
        Ly = z*vx - x*vz 
        Lz = x*vy - y*vx 

        L_t[tid, 1] += mass * Lx
        L_t[tid, 2] += mass * Ly
        L_t[tid, 3] += mass * Lz


       
    end

    K = 0.0
    Uh = 0.0
    Uc = 0.0
    L = zeros(Float64, 3)
    n_unbound = 0.0

    if (first_timestep)

        Threads.@threads for i=1:Npart

            tid = Threads.threadid()

            Ui = 0.0 # 1/2 sum_{j != i} G mi mj/rij

            # Cluster interaction potential energy 

            xi = tab_stars[i,1]
            yi = tab_stars[i,2]
            zi = tab_stars[i,3]

            for j=1:Npart
                if (j != i)
        
                    xj = tab_stars[j,1]
                    yj = tab_stars[j,2]
                    zj = tab_stars[j,3]
        
                    # Vector ri - rj
                    xij = xi - xj 
                    yij = yi - yj 
                    zij = zi - zj 
        
                    rij = sqrt(xij^2 + yij^2 + zij^2 + eps^2)
        
                    Ui -= _G*mass*mass/rij

                end
        
            end

            tab_Uint[i] = Ui

            # Unbound particles 

            x, y, z, vx, vy, vz = tab_stars[i, :]

            vc_x = vx - Vbx
            vc_y = vy - Vby
            vc_z = vz - Vbz

            Ec = 0.5 * mass * (vc_x^2 + vc_y^2 + vc_z^2) + Ui
            if (Ec >= 0.0)
                n_unbound_t[tid] += 1.0
            end

            Uc_t[tid] += 0.5 * Ui

        end

        
 

        for tid=1:Threads.nthreads()

            Uc += Uc_t[tid]
            n_unbound += n_unbound_t[tid]

        end

    else 
        Threads.@threads for i=1:Npart

            tid = Threads.threadid()

            # Cluster potential energy 
            Uint = tab_Uint[i]
            Uc_t[tid] += 0.5*Uint

            # Unbound particles 

            x, y, z, vx, vy, vz = tab_stars[i, :]
            vc_x = vx - Vbx
            vc_y = vy - Vby
            vc_z = vz - Vbz

            Ec = 0.5 * mass * (vc_x^2 + vc_y^2 + vc_z^2) + Uint

            if (Ec >= 0.0)
                n_unbound_t[tid] += 1.0
            end
    
           

        end

        for tid=1:Threads.nthreads()

            Uc += Uc_t[tid]
            n_unbound += n_unbound_t[tid]

        end
        
 

    end

    

    for tid=1:Threads.nthreads()

        K += K_t[tid]
        Uh += Uh_t[tid]
        L[1] += L_t[tid, 1]
        L[2] += L_t[tid, 2]
        L[3] += L_t[tid, 3]

    end

    tab_IOM[1] = K 
    tab_IOM[2] = Uh + Uc 
    tab_IOM[3] = K + Uh + Uc 
    tab_IOM[4] = L[1]
    tab_IOM[5] = L[2]
    tab_IOM[6] = L[3]
    tab_IOM[7] = sqrt(L[1]^2 + L[2]^2 + L[3]^2)
    tab_IOM[8] = n_unbound

end