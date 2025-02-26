using DelimitedFiles


function write_data!(time::Float64, tab_stars::Array{Float64}, tab_Uint::Array{Float64}, tab_Uc::Array{Float64})
    
    n_digits = floor(Int64,-log10(dt))+1
    namefile = folder_output*"snapshots_"*srun*"/time_"*string(round(time, digits=n_digits))*".txt"
    writedlm(namefile, [tab_stars tab_Uint tab_Uc])

    return nothing

end