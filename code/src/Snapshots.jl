using DelimitedFiles

function write_data!(time::Float64, tab_stars::Array{Float64})
    
    n_digits = floor(Int64,-log10(dt))+1
    namefile = path_dir*"data/snapshots_"*srun*"/time_"*string(round(time, digits=n_digits))*".txt"
    writedlm(namefile, tab_stars)

    return nothing

end