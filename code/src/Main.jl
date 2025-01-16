include("Args.jl")
include("Constants.jl")
include("Host.jl")
include("Cluster.jl")
include("Snapshots.jl")
include("Integrator.jl")


function main()

    tab_stars = zeros(Float64, Npart, 6)
    time = 0.0
    initialize_stars!(tab_stars)

    index = 0

    mkpath(path_dir*"data/snapshots_"*sdate*"/")

    while (time < time_end)

        println("Progress = ", time/time_end)

        # Snapshots 
        if (index % N_dt == 0)
            write_data!(time, tab_stars)
        end

        integrate_stars_leapfrog!(tab_stars)
        time += dt
        index += 1

    end


end


@time main()