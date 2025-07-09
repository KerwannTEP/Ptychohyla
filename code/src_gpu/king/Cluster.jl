using DelimitedFiles
using CSV 
using DataFrames

println("Compiling CUDA...")
@time using CUDA 

############################################################################################################################################
#  Initialization
############################################################################################################################################

function initialize_stars!(tab_stars::Array{Float64})

    # # Load King sphere in Henon units
    # if (!HAS_MULTI_MASS) 
    #     # Single mass
    #     namefile = path_dir*"data/IC/chc_king_ics_n_"*string(Npart)*".csv"
    # else 
    #     # Multi-mass (half of mass m1, half of mass m2=2*m1)
    #     namefile = path_dir*"data/IC/chc_king_ics_n_"*string(Npart)*"_multi_mass_1_2.csv"
    # end

    namefile = file_IC

    df = CSV.read(namefile, DataFrame, delim=',', header=false)

    datam = df[:, 3]
    datax = df[:, 6]
    datay = df[:, 7]
    dataz = df[:, 8]
    datavx = df[:, 9]
    datavy = df[:, 10]
    datavz = df[:, 11]

    if (CUSTOM_IC)

        custom_IC = readdlm(path_dir * "job_slurm/custom_IC.txt"; header=false)

        println("x0  [kpc]  = ", custom_IC[1,1])
        println("y0  [kpc]  = ", custom_IC[1,2])
        println("z0  [kpc]  = ", custom_IC[1,3])
        println("vx0 [km/s] = ", custom_IC[2,1])
        println("vy0 [km/s] = ", custom_IC[2,2])
        println("vz0 [km/s] = ", custom_IC[2,3])

        x0 = custom_IC[1,1] / R_HU_in_kpc  # HU
        y0 = custom_IC[1,2] / R_HU_in_kpc  # HU
        z0 = custom_IC[1,3] / R_HU_in_kpc  # HU

        vx0 = custom_IC[2,1] / V_HU_in_km_s  # HU
        vy0 = custom_IC[2,2] / V_HU_in_km_s  # HU
        vz0 = custom_IC[2,3] / V_HU_in_km_s  # HU

        Threads.@threads for i=1:Npart 

            tab_stars[i, 1] = x0 + datax[i]
            tab_stars[i, 2] = y0 + datay[i]
            tab_stars[i, 3] = z0 + dataz[i]
            tab_stars[i, 4] = vx0 + datavx[i]
            tab_stars[i, 5] = vy0 + datavy[i]
            tab_stars[i, 6] = vz0 + datavz[i]
            tab_stars[i, 7] = datam[i]

        end

    else

        vcirc = host.circular_velocity(d_host)
        println("Vc [HU] = ", vcirc)
        println("Vc [km/s] = ", vcirc * V_HU_in_km_s)

        Threads.@threads for i=1:Npart 

            tab_stars[i, 1] = d_host + datax[i] # Cluster on the x>0 x-axis
            tab_stars[i, 2] = datay[i]
            tab_stars[i, 3] = dataz[i]
            tab_stars[i, 4] = datavx[i]
            tab_stars[i, 5] = vcirc + datavy[i] # Circular velocity: cluster goes in the y-direction
            tab_stars[i, 6] = datavz[i]
            tab_stars[i, 7] = datam[i]

        end

    end


end

function initialize_stars_restart!(tab_stars::Array{Float64})

    # Load last saved snapshot 
    # Fill data
    # Return restart time

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

    namefile = listFiles[p[nsnap]]
    time = tabt[p[nsnap]]

    data_stars = readdlm(namefile) # x, y, z, vx, vy, vz, m, Uint, Uc
    # Array of size (Npart, 8), in Henon units

    Threads.@threads for i=1:Npart 

        tab_stars[i, 1] = data_stars[i, 1]
        tab_stars[i, 2] = data_stars[i, 2]
        tab_stars[i, 3] = data_stars[i, 3]
        tab_stars[i, 4] = data_stars[i, 4]
        tab_stars[i, 5] = data_stars[i, 5]
        tab_stars[i, 6] = data_stars[i, 6]
        tab_stars[i, 7] = data_stars[i, 7]

    end

    return time

end


############################################################################################################################################
# Accelerations and interaction energies
############################################################################################################################################

# https://developer.nvidia.com/gpugems/gpugems3/part-v-physics-simulation/chapter-31-fast-n-body-simulation-cuda
function compute_acc_Uint_int_gpu!(tab_pos_x::CuDeviceArray{T}, tab_pos_y::CuDeviceArray{T}, tab_pos_z::CuDeviceArray{T}, tab_m::CuDeviceArray{T}, 
                            tab_acc_x::CuDeviceArray{T}, tab_acc_y::CuDeviceArray{T}, tab_acc_z::CuDeviceArray{T}, tab_Uint::CuDeviceArray{T}) where T


    # We are within a thread, within a block 
    # Each thread correspond for one star, for which we compute the acceleration 
    # The acceleration is induced by the star's interaction with all N stars (contribution with itself vanishes for a softened potential)
    # Cut each thread block into square tiles of size blockSize
    # Tile decomposition avoids loading huge chuncks of data (of order Npart^2)

    gtid = threadIdx().x  + (blockIdx().x - 1) * blockDim().x # Index of the particle for which we compute the acceleration

    # Do not go out-of-bounds for the particles for which we compute the acceleration
    if (gtid <= Npart)

        x = tab_pos_x[gtid]
        y = tab_pos_y[gtid]
        z = tab_pos_z[gtid]
        m = tab_m[gtid]

    end

    ax = 0.0
    ay = 0.0
    az = 0.0
    Uint = 0.0

    tile = 1
    i = 1

    # https://cuda.juliagpu.org/stable/development/kernel/#Shared-memory
    # https://cuda.juliagpu.org/stable/api/kernel/#Shared-memory
    shPosition = CuDynamicSharedArray(T, (blockDim().x, 4))

    # Interaction of the star with the N stars of the clusters (include itself)
    while (i <= Npart)

        idx = threadIdx().x  + (tile - 1) * blockDim().x

        # Filling shared memory of tile
        # Be careful not to go out-of-bounds for the last tile (i.e. when looking at particles accelerating this thread's particle)
        
        if (idx <= Npart)

            shPosition[threadIdx().x, 1] = tab_pos_x[idx]
            shPosition[threadIdx().x, 2] = tab_pos_y[idx]
            shPosition[threadIdx().x, 3] = tab_pos_z[idx]
            shPosition[threadIdx().x, 4] = tab_m[idx]

        end

        length_tile = min(blockDim().x, Npart - i + 1)

        sync_threads()

        # Calculation of acceleration within tile sequentially
        # Each thread is performed in parallel
        # Use shared memory for faster access

        # Do not go out-of-bounds for the particles for which we compute the acceleration
        if (gtid <= Npart)

            # Be careful not to go out-of-bounds for the last tile (i.e. when looking at particles accelerating this thread's particle)
            for j=1:length_tile

                id_part = j + (tile - 1) * blockDim().x # Index of particle which interacts with gtid

                if (gtid != id_part) # no self-interaction
                    dx = shPosition[j,1] - x
                    dy = shPosition[j,2] - y
                    dz = shPosition[j,3] - z
                    mj = shPosition[j,4]

                    dr = sqrt(dx*dx + dy*dy + dz*dz + eps*eps)

                    ar = _G*mj/dr^3 
                    Uij = - _G*m*mj/dr # Interaction potential between gtid and 
                    
                    ax += ar * dx
                    ay += ar * dy
                    az += ar * dz
                    Uint += Uij
                end


            end

        end
        
        sync_threads()

        # Go to next tile
        i += blockDim().x
        tile += 1
    end

    # Save acceleration into global array for the star in the current thread
    # Do not go out-of-bounds for the particles for which we compute the acceleration
    if (gtid <= Npart)
        tab_acc_x[gtid] = ax
        tab_acc_y[gtid] = ay
        tab_acc_z[gtid] = az
        tab_Uint[gtid] = Uint
    end

    return nothing

end 


function tab_acc_Uint_int_gpu!(tab_acc::Array{Float64}, tab_Uint::Array{Float64}, tab_pos::Array{Float64}, tab_m::Array{Float64})

    tab_pos_x = tab_pos[:,1]
    tab_pos_y = tab_pos[:,2]
    tab_pos_z = tab_pos[:,3]

    numblocks = ceil(Int64, Npart/nbThreadsPerBlocks)

    dev_tab_pos_x = CuArray(tab_pos_x)
    dev_tab_pos_y = CuArray(tab_pos_y)
    dev_tab_pos_z = CuArray(tab_pos_z)
    dev_tab_m = CuArray(tab_m)

    dev_tab_acc_x = CuArray(zeros(Float64, Npart))
    dev_tab_acc_y = CuArray(zeros(Float64, Npart))
    dev_tab_acc_z = CuArray(zeros(Float64, Npart))
    dev_tab_Uint = CuArray(zeros(Float64, Npart))


    @cuda threads=nbThreadsPerBlocks blocks=numblocks shmem=4*nbThreadsPerBlocks*sizeof(Float64) compute_acc_Uint_int_gpu!(dev_tab_pos_x, dev_tab_pos_y, dev_tab_pos_z, dev_tab_m,
                                                                                                                dev_tab_acc_x, dev_tab_acc_y, dev_tab_acc_z, dev_tab_Uint)

    tab_acc_x = Array(dev_tab_acc_x)
    tab_acc_y = Array(dev_tab_acc_y)
    tab_acc_z = Array(dev_tab_acc_z)
    tab_Uint_load = Array(dev_tab_Uint)

    Threads.@threads for i=1:Npart 
        tab_acc[i, 1] = tab_acc_x[i]
        tab_acc[i, 2] = tab_acc_y[i]
        tab_acc[i, 3] = tab_acc_z[i]
        tab_Uint[i] = tab_Uint_load[i]

    end 

    return nothing

end