using DelimitedFiles

function write_data!(time::Float64, tab_stars::Array{Float64})

    # println(sdate)
    
    n_digits = floor(Int64,-log10(dt))+1
    # println((time, n_digits))
    namefile = path_dir*"data/snapshots_"*sdate*"/time_"*string(round(time, digits=n_digits))*".txt"
    writedlm(namefile, tab_stars)

    return nothing

end