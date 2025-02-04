using SpecialFunctions

# Host potential (Henon units)

####################################################################################
# Bulge
# CHC II appendix C1
####################################################################################

# CAREFUL: gamma(s, x) = int_x^infty ... is the upper incomplete gamma function

# Enclosed mass of the bulge

# Something wrong here with potential

function gamma_low_sx(xr)

    return gamma_s - gamma(s_bulge, xr)

end

function Menc_bulge(r::Float64)

    xr = (r/rc_bulge)^2
    return M_bulge * gamma_low_sx(xr)/gamma_s
end

function force_bulge(x::Float64, y::Float64, z::Float64)

    r = sqrt(x*x + y*y + z*z)
    xr = (r/rc_bulge)^2
    gamma_s_x = gamma_low_sx(xr)



    Menc = M_bulge * gamma_s_x/gamma_s
    fr = - _G * Menc/r^2

    # println(fr)

    Fx = fr * x/r
    Fy = fr * y/r
    Fz = fr * z/r

    return Fx, Fy, Fz 

end

function psi_bulge(r::Float64)

    xr = (r/rc_bulge)^2
    gamma_upper = gamma(s_bulge, xr)
    gamma_lower = gamma_s - gamma_upper

    Menc = M_bulge * gamma_lower/gamma_s
    rho0 = M_bulge/(2.0 * gamma_s)


    # Inner 
    psi_inner = -_G * Menc/r 

    # outer
    psi_outer = -2.0*pi*_G * rho0*rc_bulge^2 * gamma_upper

    return psi_inner + psi_outer

end

function dpsidr_bulge(r::Float64)

    xr = (r/rc_bulge)^2
    gamma_s_x = gamma_low_sx(xr)
    return _G * M_bulge/r^2 * gamma_s_x/gamma_s
end


####################################################################################
# Disk
# Miyamoto-Nagai disk
# CHC II appendix C1
####################################################################################

function force_disk(x::Float64, y::Float64, z::Float64)

    R = sqrt(x*x + y*y)

    fR = - _G * M_disk * R/(R^2 + (sqrt(z^2 + b_disk^2)+a_disk)^2)^(3/2)
    fz = - _G * M_disk * z * (sqrt(z^2 + b_disk^2) + a_disk)/(sqrt(z^2+b_disk^2)*(R^2 + (sqrt(z^2 + b_disk^2)+a_disk)^2)^(3/2))

    Fx = fR * x/R 
    Fy = fR * y/R 
    Fz = fz

    return Fx, Fy, Fz

end

function psi_disk(R::Float64, z::Float64)

    return _G * M_disk/sqrt(R^2 + (sqrt(z^2 + b_disk^2) + a_disk)^2)

end

function dpsidR_disk(R::Float64, z::Float64)

    return _G * M_disk * R/(R^2 + (sqrt(z^2 + b_disk^2)+a_disk)^2)^(3/2)
end



####################################################################################
# Dark halo (NFW)
# https://en.wikipedia.org/wiki/Navarro%E2%80%93Frenk%E2%80%93White_profile
####################################################################################

function _g(x::Float64)

    return log(1.0+x) - x/(1+x)
end

function Menc_halo(r::Float64)

    return Mvir * _g(r/Rs)/g_c
    
end

function force_halo(x::Float64, y::Float64, z::Float64)

    r = sqrt(x*x + y*y + z*z)

    Menc = Menc_halo(r)
    fr = - _G * Menc/r^2

    Fx = fr * x/r
    Fy = fr * y/r
    Fz = fr * z/r

    return Fx, Fy, Fz 

end

function psi_halo(r::Float64)

    x = r/Rs
    cst = _G*Mvir/Rs * 1.0/g_c

    return -cst/x * log(1+x)

end

function dpsidr_halo(r::Float64)

    Menc = Menc_halo(r)
    return _G * Menc/r^2
end


####################################################################################
# Total force
####################################################################################

function force_host(x::Float64, y::Float64, z::Float64)

    Fx_b, Fy_b, Fz_b = force_bulge(x, y, z)
    Fx_d, Fy_d, Fz_d = force_disk(x, y, z)
    Fx_h, Fy_h, Fz_h = force_halo(x, y, z)


    Fx = Fx_h + Fx_d + Fx_b
    Fy = Fy_h + Fy_d + Fy_b
    Fz = Fz_h + Fz_d + Fz_b

    return Fx, Fy, Fz 

end

####################################################################################
# Circular velocity
####################################################################################

# z=0
function circular_velocity(R::Float64)

    dpsidR_b = dpsidr_bulge(R)
    dpsidR_d = dpsidR_disk(R, 0.0)
    dpsidR_h = dpsidr_halo(R)

    dpsidR = dpsidR_h + dpsidR_d + dpsidR_b

    return sqrt(R*dpsidR)

end