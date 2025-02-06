include("Args.jl")
include("Constants.jl")
include("Host.jl")
include("Snapshots.jl")
include("Integrator_one_particle.jl")


function main()

    tab_stars = zeros(Float64, 6) # (x, y, z, vx, vy, vz)
    tab_IOM = zeros(Float64, 7) # K, U, Etot, Lx, Ly, Lz, L
    time = 0.0

    # Initialize 

    tab_stars[1] = d_host
    tab_stars[5] = circular_velocity(d_host)


    index = 0

    # mkpath(path_dir*"data/snapshots_"*srun*"/")

    while (time < time_end)
        

        compute_IOM!(tab_stars, tab_IOM)

        # Snapshots 
        if (index % N_dt == 0)
            write_data!(time, tab_stars, tab_IOM)
            println("Progress = ", round(time/time_end, digits=4), " | Energy = ", tab_IOM[3])
        end

        integrate_one_leapfrog!(tab_stars)
        time += dt
        index += 1

    end

    # Last snapshot
    compute_IOM!(tab_stars, tab_IOM)
    write_data!(time, tab_stars, tab_IOM)

end


@time main()