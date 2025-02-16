using DelimitedFiles


function write_data!(time::Float64, tab_stars::Array{Float64}, tab_Uint::Array{Float64})
    
    n_digits = floor(Int64,-log10(dt))+1
    namefile = folder_output*"snapshots_"*srun*"/time_"*string(round(time, digits=n_digits))*".txt"
    writedlm(namefile, [tab_stars tab_Uint])

    return nothing

end


# Compute along force calculation for efficiency
# GPU acceleration for Uint 
# Use IOM julia script from NbodyPostTreatment folder

function compute_IOM!(tab_stars::Array{Float64}, tab_Uint::Array{Float64}, tab_IOM::Array{Float64})

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

        # Cluster potential energy 
        Uint = tab_Uint[i]
        Uc_t[tid] += 0.5*Uint

        # Angular momenta 
        Lx = y*vz - z*vy 
        Ly = z*vx - x*vz 
        Lz = x*vy - y*vx 

        L_t[tid, 1] += mass * Lx
        L_t[tid, 2] += mass * Ly
        L_t[tid, 3] += mass * Lz

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


    K = 0.0
    Uh = 0.0
    Uc = 0.0
    L = zeros(Float64, 3)
    n_unbound = 0.0

    for tid=1:Threads.nthreads()

        K += K_t[tid]
        Uh += Uh_t[tid]
        Uc += Uc_t[tid]
        L[1] += L_t[tid, 1]
        L[2] += L_t[tid, 2]
        L[3] += L_t[tid, 3]
        n_unbound += n_unbound_t[tid]

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