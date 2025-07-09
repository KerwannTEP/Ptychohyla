include("Args.jl")
include("Constants.jl")
include("../Host.jl")
include("Snapshots.jl")
include("Integrator_one_particle.jl")

println("Run ID : ", srun)

function main()

    mkpath(folder_output*"snapshots_"*srun*"/")

    timing_start = now()

    tab_stars = zeros(Float64, 6) # (x, y, z, vx, vy, vz)
    tab_IOM = zeros(Float64, 7) # K, U, Etot, Lx, Ly, Lz, L
    time = 0.0

    # Initialize 

    tab_stars[1] = d_host
    tab_stars[5] = host.circular_velocity(d_host)

    index = 0

    while (time < time_end)
        

        compute_IOM!(tab_stars, tab_IOM)

        # Snapshots 
        if (index % N_dt == 0)
            write_data!(time, tab_stars, tab_IOM)
            println("Progress = ", round(time/time_end, digits=4), " | Energy = ", tab_IOM[3])
        end

        integrate_one_leapfrog!(tab_stars, host)
        # integrate_one_yoshida!(tab_stars, host)
        time += dt
        index += 1

    end

    timing_end = now()

    # https://stackoverflow.com/questions/41293747/round-julias-millisecond-type-to-nearest-second-or-minute
    dtim = timing_end - timing_start

    println("-----------------------")

    # https://discourse.julialang.org/t/how-to-convert-period-in-milisecond-to-minutes-seconds-hour-etc/2423/6
    dt_v = Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(dtim)))
    println("Simulation took : ", dt_v)

end


@time main()