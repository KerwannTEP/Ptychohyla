using Pkg; Pkg.activate("../../../packages/julia_1-9-3/JuEnvGPU")

println("Nb threads = ", Threads.nthreads())

include("Args.jl")
include("Constants.jl")
include("../Host.jl")
include("Cluster.jl")
include("../Snapshots.jl")
include("../Integrator.jl")

println("Run ID : ", srun)


function main()

    mkpath(folder_output*"snapshots_"*srun*"/")

    timing_start = now()

    tab_stars = zeros(Float64, Npart, 6) # (x, y, z, vx, vy, vz)
    tab_acc = zeros(Float64, Npart, 3)  
    tab_Uint = zeros(Float64, Npart) # (Ui, Ei_wrt_cluster)
    
    index = 0
    time = 0.0
    first_timestep = true

    initialize_stars!(tab_stars)

    while (time < time_end)
        

        integrate_stars_leapfrog!(index, time, tab_stars, tab_acc, tab_Uint, first_timestep)
        time += dt
        index += 1

        if (first_timestep)
            first_timestep = false
        end

    end

    # Last snapshot
    write_data!(time, tab_stars, tab_Uint)  

    timing_end = now()

    # https://stackoverflow.com/questions/41293747/round-julias-millisecond-type-to-nearest-second-or-minute
    dtim = timing_end - timing_start

    println("-----------------------")

    # https://discourse.julialang.org/t/how-to-convert-period-in-milisecond-to-minutes-seconds-hour-etc/2423/6
    dt_v = Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(dtim)))
    println("Simulation took : ", dt_v)

    # println("-----------------------")
    # println("Post-treatment (energy, etc) ...")

    # WRITE IOM IN POST-PROCESSING

end


main()