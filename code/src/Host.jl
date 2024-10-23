
# Host potential (Henon units)
# https://en.wikipedia.org/wiki/Navarro%E2%80%93Frenk%E2%80%93White_profile




function psi_NFW(r::Float64)

    pref = 4*pi*_G*rho0_host*Rs^3
    return -pref/r * log(1 + r/Rs)
end

function gradr_psi_NFW(r::Float64)

    pref = 4*pi*_G*rho0_host*Rs^3
    return pref * (-1.0/(Rs*r*(1+r/Rs)) + log(1+r/Rs)/r^2)
end

# F_{NFW,k} = - m_k gradr_psi_NFW(r_k) \vec{r_k}
function force_NFW(x::Float64, y::Float64, z::Float64)

    r = sqrt(x*x + y*y + z*z)
    gradr_psi = gradr_psi_NFW(r)

    Fx = -mass * gradr_psi * x/r
    Fy = -mass * gradr_psi * y/r
    Fz = -mass * gradr_psi * z/r

    return Fx, Fy, Fz 
end


function circular_velocity(r::Float64)

    return sqrt(r*gradr_psi_NFW(r))

end