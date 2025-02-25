using DelimitedFiles


function write_data!(time::Float64, tab_stars::Array{Float64}, tab_Uint::Array{Float64}, tab_Uc::Array{Float64})
    
    n_digits = floor(Int64,-log10(dt))+1
    namefile = folder_output*"snapshots_"*srun*"/time_"*string(round(time, digits=n_digits))*".txt"
    writedlm(namefile, [tab_stars tab_Uint tab_Uc])

    return nothing

end

function post_treatment!()

    # Read data

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

    sortedFiles = Array{String}(undef, nsnap)
    tabtsort = zeros(Float64, nsnap)

    for i=1:nsnap 

        tabtsort[i] = tabt[p[i]]
        sortedFiles[i] = listFiles[p[i]]

    end

    # Compute useful quantities

    tab_IOM = zeros(Float64, nsnap, 8) # time, K, U, E_tot, Lx, Ly, Lz, nb_unbound

    for isnap=1:nsnap

        namefile = sortedFiles[isnap]
        time = tabtsort[isnap]

        data_stars = readdlm(namefile) # x, y, z, vx, vy, vz, Uint

        tab_vb_t = zeros(Float64, Threads.nthreads(), 3)
        K_t =  zeros(Float64, Threads.nthreads())
        Uc_t =  zeros(Float64, Threads.nthreads())
        Uh_t =  zeros(Float64, Threads.nthreads())
        L_t =  zeros(Float64, Threads.nthreads(), 3)

        Threads.@threads for i=1:Npart 

            tid = Threads.threadid()
            x, y, z, vx, vy, vz, Uint = data_stars[i, :]
            r = sqrt(x^2 + y^2 + z^2)
            R = sqrt(x^2 + y^2)
            psi_xyz = psi_halo(r) + psi_disk(R, z) + psi_bulge(r) 

            v2 = vx^2 + vy^2 + vz^2

            # Barycenter 
            tab_vb_t[tid, 1] += vx
            tab_vb_t[tid, 2] += vy
            tab_vb_t[tid, 3] += vz

            # Kinetic energy 
            K_t[tid] += 0.5 * mass * v2 

            # Host potential energy
            Uh_t[tid] += mass * psi_xyz

            # Cluster potential energy 
            Uc_t[tid] += 0.5*Uint

                
            # Angular momenta 
            Lx = y*vz - z*vy 
            Ly = z*vx - x*vz 
            Lz = x*vy - y*vx 

            L_t[tid, 1] += mass * Lx
            L_t[tid, 2] += mass * Ly
            L_t[tid, 3] += mass * Lz

        end
       
        Vbx = 0.0
        Vby = 0.0
        Vbz = 0.0

        K = 0.0
        U = 0.0
        L = zeros(Float64, 3)


        for tid=1:Threads.nthreads()
            Vbx += tab_vb_t[tid, 1]/Npart
            Vby += tab_vb_t[tid, 2]/Npart
            Vbz += tab_vb_t[tid, 3]/Npart

            K += K_t[tid]
            U += Uh_t[tid] + Uc_t[tid]
            L[1] += L_t[tid, 1]
            L[2] += L_t[tid, 2]
            L[3] += L_t[tid, 3]
    
        end

        Etot = K + U

        # Compute number of unbound stars
        n_unbound_t = zeros(Float64, Threads.nthreads())

        Threads.@threads for i=1:Npart

            tid = Threads.threadid()
    
            # Unbound particles 
            x, y, z, vx, vy, vz, Uint = data_stars[i, :]
            vc_x = vx - Vbx
            vc_y = vy - Vby
            vc_z = vz - Vbz
    
            Ec = 0.5 * mass * (vc_x^2 + vc_y^2 + vc_z^2) + Uint
    
            if (Ec >= 0.0)
                n_unbound_t[tid] += 1.0
            end
    
        end

        n_unbound = 0.0

        for tid=1:Threads.nthreads()
            n_unbound += n_unbound_t[tid]
        end
    
        tab_IOM[isnap,1] = time
        tab_IOM[isnap,2] = K 
        tab_IOM[isnap,3] = U
        tab_IOM[isnap,4] = Etot
        tab_IOM[isnap,5] = L[1]
        tab_IOM[isnap,6] = L[2]
        tab_IOM[isnap,7] = L[3]
        tab_IOM[isnap,8] = n_unbound

        

    end

    # Save quantities
    namefile = folder_output*"/iom_snapshots_"*srun*".txt"
    writedlm(namefile, tab_IOM)

end
