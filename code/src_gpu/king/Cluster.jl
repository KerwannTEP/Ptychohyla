using DelimitedFiles
using CSV 
using DataFrames

println("Compiling CUDA...")
@time using CUDA 

function initialize_stars!(tab_stars::Array{Float64})

    # Load King sphere in Henon units
    namefile = path_dir*"data/IC/chc_king_ics_n_"*string(Npart)*".csv"
    df = CSV.read(namefile, DataFrame, delim=',', header=false)


    datax = df[:, 6]
    datay = df[:, 7]
    dataz = df[:, 8]
    datavx = df[:, 9]
    datavy = df[:, 10]
    datavz = df[:, 11]


    vcirc = circular_velocity(d_host)
    println("Vc [HU] = ", vcirc)
    println("Vc [km/s] = ", vcirc * V_HU_in_km_s)

    Threads.@threads for i=1:Npart 

        tab_stars[i, 1] = d_host + datax[i] # Cluster on the x>0 x-axis
        tab_stars[i, 2] = datay[i]
        tab_stars[i, 3] = dataz[i]
        tab_stars[i, 4] = datavx[i]
        tab_stars[i, 5] = vcirc + datavy[i] # Circular velocity: cluster goes in the y-direction
        tab_stars[i, 6] = datavz[i]

    end


end

# ADAPT FOR GPU
# Compute IOM separately
# https://developer.nvidia.com/gpugems/gpugems3/part-v-physics-simulation/chapter-31-fast-n-body-simulation-cuda
function compute_acc_int_gpu!(tab_pos_x::CuDeviceArray{T}, tab_pos_y::CuDeviceArray{T}, tab_pos_z::CuDeviceArray{T}, 
                            tab_acc_x::CuDeviceArray{T}, tab_acc_y::CuDeviceArray{T}, tab_acc_z::CuDeviceArray{T}, 
                            tab_Uint::CuDeviceArray{T}, 
                            tab_potential_energy_blocks::CuDeviceArray{T}, nbThreadsPerBlocks::Int64, N::Int64, nbPairs::Int64) where T


    # We are within a thread, within a block 
    # Each thread correspond for one star, for which we compute the acceleration 
    # The acceleration is induced by the star's interaction with all N stars (contribution with itself vanishes for a softened potential)
    # Cut each thread block into square tiles of size blockSize
    # Tile decomposition avoids loading huge chuncks of data (of order Npart^2)

    gtid = threadIdx().x  + (blockIdx().x - 1) * blockDim().x
    x = tab_pos_x[gtid]
    y = tab_pos_y[gtid]
    z = tab_pos_z[gtid]

    ax = 0.0
    ay = 0.0
    az = 0.0

    tile = 1
    i = 1

    shPositionX = CuDynamicSharedArray(T, nbThreadsPerBlocks)
    shPositionY = CuDynamicSharedArray(T, nbThreadsPerBlocks)
    shPositionZ = CuDynamicSharedArray(T, nbThreadsPerBlocks)

    # Interaction of the star with the N stars of the clusters (include itself)
    while (i <= Npart)

        idx = threadIdx().x  + (tile - 1) * blockDim().x

        # Filling shared memory of tile

        shPositionX[threadIdx().x] = tab_pos_x[idx]
        shPositionY[threadIdx().x] = tab_pos_y[idx]
        shPositionZ[threadIdx().x] = tab_pos_z[idx]

        sync_threads()

        # Calculation of acceleration within tile sequentially
        # Each thread is performed in parallel
        # Use shared memory for faster access

        for j=1:blockDim().x

            dx = x - shPositionX[j]
            dy = y - shPositionY[j]
            dz = z - shPositionZ[j]

            dr = sqrt(dx*dx + dy*dy + dz*dz + eps*eps)

            ar = _G*mass/dr^3 
            
            ax += ar * dx
            ay += ar * dy
            az += ar * dz


        end
        
        sync_threads()

        # Go to next tile
        i += blockDim().x
        tile += 1
    end

    # Save acceleration into global array for the star in the current thread
    tab_acc_x[gtid] = ax
    tab_acc_y[gtid] = ay
    tab_acc_z[gtid] = az

    return nothing

end 


function tab_acc_int_gpu!(tab_acc::Array{Float64}, tab_pos::Array{Float64})

    tab_pos_x = tab_pos[:,1]
    tab_pos_y = tab_pos[:,2]
    tab_pos_z = tab_pos[:,3]

    # # Number of pairs (i,j) with repetitions (i.e. we count {ij} and {ji} as two distinct pairs, and append a 0.5 in front of the sum)
    # # For convenient we include {i,i} twice, since a_{ii} = 0 for a softened potential anyway
    # nbPairs = Npart*Npart

    # numblocks = min(40, ceil(Int64, nbPairs/nbThreadsPerBlocks))
    numblocks = min(40, ceil(Int64, Npart/nbThreadsPerBlocks))

    dev_tab_pos_x = CuArray(tab_pos_x)
    dev_tab_pos_y = CuArray(tab_pos_y)
    dev_tab_pos_z = CuArray(tab_pos_z)

    dev_tab_acc_x = CuArray(zeros(Float64, Npart))
    dev_tab_acc_y = CuArray(zeros(Float64, Npart))
    dev_tab_acc_z = CuArray(zeros(Float64, Npart))


    # dev_tab_potential_energy_block = CuArray(zeros(Float64, numblocks))

    @cuda threads=nbThreadsPerBlocks blocks=numblocks shmem=3*nbThreadsPerBlocks*sizeof(Float64) compute_acc_int_gpu!(dev_tab_pos_x, dev_tab_pos_y, dev_tab_pos_z,
                                                                                                                dev_tab_acc_x, dev_tab_acc_y, dev_tab_acc_z,
                                                                                                                dev_tab_Uint,
                                                                                                                dev_tab_potential_energy_block, nbThreadsPerBlocks, Npart, nbPairs)



    tab_acc_x = Array(dev_tab_acc_x)
    tab_acc_y = Array(dev_tab_acc_y)
    tab_acc_z = Array(dev_tab_acc_z)

    Threads.@threads for i=1:Npart 
        tab_acc[i, 1] = tab_acc_x[i]
        tab_acc[i, 2] = tab_acc_y[i]
        tab_acc[i, 3] = tab_acc_z[i]

        tid = Threads.threadid()

    end 

    return nothing

end
