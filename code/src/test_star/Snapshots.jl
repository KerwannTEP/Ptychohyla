using DelimitedFiles


function write_data!(time::Float64, tab_stars::Array{Float64}, tab_IOM::Array{Float64})
    
    n_digits = floor(Int64,-log10(dt))+1
    namefile = folder_output*"snapshots_"*srun*"/pos_snapshots_"*srun*".txt"

    io = open(namefile, "a")
    writedlm(io, transpose([round(time, digits=n_digits); tab_stars]))
    close(io)

    # https://discourse.julialang.org/t/adding-to-existing-txt-file-created-by-writedlm/6907/2
    namefile_iom = folder_output*"snapshots_"*srun*"/iom_snapshots_"*srun*".txt"
    io2 = open(namefile_iom, "a")
    writedlm(io2, transpose([round(time, digits=n_digits); tab_IOM]))
    close(io2)

    return nothing

end


function compute_IOM!(tab_stars::Array{Float64}, tab_IOM::Array{Float64})

    # tab_IOM = zeros(Float64, 7) # K, U, Etot, Lx, Ly, Lz, L

    x, y, z, vx, vy, vz = tab_stars[:]

    r = sqrt(x^2 + y^2 + z^2)
    R = sqrt(x^2 + y^2)
    v2 = vx^2 + vy^2 + vz^2
    psi_xyz = psi_halo(r) + psi_disk(R, z) + psi_bulge(r) ##

    # Kinetic energy 
    K = 0.5 * v2 

    # Host potential energy
    Uh = psi_xyz

    # Angular momenta 
    Lx = y*vz - z*vy 
    Ly = z*vx - x*vz 
    Lz = x*vy - y*vx 

    tab_IOM[1] = K 
    tab_IOM[2] = Uh 
    tab_IOM[3] = K + Uh
    tab_IOM[4] = Lx
    tab_IOM[5] = Ly
    tab_IOM[6] = Lz
    tab_IOM[7] = sqrt(Lx^2 + Ly^2 + Lz^2)

end