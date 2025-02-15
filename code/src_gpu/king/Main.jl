using Pkg; Pkg.activate("../../../packages/julia_1-9-3/JuEnvGPU")

println("Nb threads = ", Threads.nthreads())

include("Args.jl")
include("Constants.jl")
include("../Host.jl")
include("Cluster.jl")
include("../Snapshots.jl")
include("../Integrator.jl")

println("Run ID :", srun)


function main()

    timing_start = now()

    tab_stars = zeros(Float64, Npart,6) # (x, y, z, vx, vy, vz)
    tab_IOM = zeros(Float64, 7) # K, U, Etot, Lx, Ly, Lz, L
    time = 0.0
    initialize_stars!(tab_stars)

    # Compute barycenter

    tab_bary = zeros(Float64, 3)

    index = 0

    mkpath(folder_output*"snapshots_"*srun*"/")

    tab_acc = zeros(Float64, Npart, 3)  

    first_timestep = true


    while (time < time_end)
        

        # compute_IOM!(tab_stars, tab_IOM)
        # compute_bary!(tab_stars, tab_bary)


        # Snapshots 
        if (index % N_dt == 0)
            write_data!(time, tab_stars, tab_IOM, tab_bary)
            # println("Progress = ", round(time/time_end, digits=4), " | Energy = ", tab_IOM[3])
      
        end

        integrate_stars_leapfrog!(tab_stars, tab_acc, first_timestep)
        time += dt
        index += 1

        if (first_timestep)
            first_timestep = false
        end

    end

    # Last snapshot
    # compute_IOM!(tab_stars, tab_IOM)
    # compute_bary!(tab_stars, tab_bary)
    write_data!(time, tab_stars, tab_IOM, tab_bary)

    timing_end = now()

    # https://stackoverflow.com/questions/41293747/round-julias-millisecond-type-to-nearest-second-or-minute
    dtim = timing_end - timing_start

    println("-----------------------")

    # https://discourse.julialang.org/t/how-to-convert-period-in-milisecond-to-minutes-seconds-hour-etc/2423/6
    dt_v = Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(dtim)))
    println("Simulation took : ", dt_v)

end


main()