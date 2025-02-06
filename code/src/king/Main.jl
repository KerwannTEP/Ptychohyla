include("Args.jl")
include("Constants.jl")
include("Host.jl")
include("Cluster.jl")
include("Snapshots.jl")
include("Integrator.jl")


function main()

    tab_stars = zeros(Float64, Npart,6) # (x, y, z, vx, vy, vz)
    tab_IOM = zeros(Float64, 7) # K, U, Etot, Lx, Ly, Lz, L
    time = 0.0
    initialize_stars!(tab_stars)

    # Compute barycenter

    tab_bary = zeros(Float64, 3)

    index = 0

    mkpath(path_dir*"data/snapshots_"*srun*"/")

    while (time < time_end)
        

        compute_IOM!(tab_stars, tab_IOM)
        compute_bary!(tab_stars, tab_bary)


        # Snapshots 
        if (index % N_dt == 0)
            write_data!(time, tab_stars, tab_IOM, tab_bary)
            println("Progress = ", round(time/time_end, digits=4), " | Energy = ", tab_IOM[3])
      
        end

        integrate_stars_leapfrog!(tab_stars)
        # integrate_stars_yoshida!(tab_stars)
        time += dt
        index += 1

    end

    # Last snapshot
    compute_IOM!(tab_stars, tab_IOM)
    compute_bary!(tab_stars, tab_bary)
    write_data!(time, tab_stars, tab_IOM, tab_bary)

end


@time main()