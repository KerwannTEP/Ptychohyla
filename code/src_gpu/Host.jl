include("host/MW2014.jl")
include("host/MW2022.jl")

# https://discourse.julialang.org/t/functions-as-a-struct-field/107013
struct Host{T1<:Function, T2<:Function, T3<:Function}

    psi_host::T1
    circular_velocity::T2
    acc_host::T3

end

if (HAS_HOST)

    if (HOST_TYPE == "MW2014")
        const host = Host(psi_host_2014, circular_velocity_2014, acc_host_2014)
    elseif (HOST_TYPE == "MW2022")
        const host = Host(psi_host_2022, circular_velocity_2022, acc_host_2022)
    end

end