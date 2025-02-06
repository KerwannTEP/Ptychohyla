using SpecialFunctions

# Host potential (Henon units)

####################################################################################
# Bulge
# CHC II appendix C1
####################################################################################

# CAREFUL: gamma(s, x) = int_x^infty ... is the upper incomplete gamma function

# Enclosed mass of the bulge

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

    Fx = mass * fr * x/r
    Fy = mass * fr * y/r
    Fz = mass * fr * z/r

    return Fx, Fy, Fz 

end

# This is wrong ?
function psi_bulge(r::Float64)

    xr = (r/rc_bulge)^2

    # Inner potential
    gamma_lower = gamma_s - gamma(s_bulge, xr)
    psi_inner = -_G*M_bulge/r * gamma_lower/gamma_s

    # Outer potential 
    gamma_upper = gamma(1.0-alpha_bulge/2.0, xr)
    psi_outer = -_G*M_bulge/rc_bulge * gamma_upper/gamma_s
     
    return psi_inner + psi_outer

end

function dpsidr_bulge(r::Float64)

    xr = (r/rc_bulge)^2
    gamma_s_x = gamma_low_sx(xr)
    return _G * M_bulge/r^2 * gamma_s_x/gamma_s
end

function rho_bulge(r::Float64)

    rho0 = M_bulge/(2.0*pi*rc_bulge^3*gamma_s)
    return rho0 * (rc_bulge/r)^(alpha_bulge) * exp(-(r/rc_bulge)^2)

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

    Fx = mass * fR * x/R 
    Fy = mass * fR * y/R 
    Fz = mass * fz

    return Fx, Fy, Fz

end

function psi_disk(R::Float64, z::Float64)

    return -_G * M_disk/sqrt(R^2 + (sqrt(z^2 + b_disk^2) + a_disk)^2)

end

function dpsidR_disk(R::Float64, z::Float64)

    return _G * M_disk * R/(R^2 + (sqrt(z^2 + b_disk^2)+a_disk)^2)^(3/2)
end

# https://galaxiesbook.org/chapters/II-01.-Flattened-Mass-Distributions.html
# eq (8.20)
function rho_disk(R::Float64, z::Float64)

    num = a_disk*R^2 + (3.0*sqrt(z^2+b_disk^2)+a_disk) * (sqrt(z^2+b_disk^2)+a_disk)^2
    den = (R^2 + (sqrt(z^2+b_disk^2)+a_disk)^2)^(5/2) * (z^2+b_disk^2)^(3/2)

    return b_disk^2*M_disk/(4*pi) * num/den

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

    Fx = mass * fr * x/r
    Fy = mass * fr * y/r
    Fz = mass * fr * z/r

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

# https://galaxiesbook.org/chapters/I-01.-Potential-Theory-and-Spherical-Mass-Distributions.html
# NFW 
# alpha=1, beta=3
function rho_halo(r::Float64)

    #4*pi*rho0*Rs^3 = Mvir/g_c

    rho0 = Mvir / (g_c * 4*pi*Rs^3)
    return rho0 * Rs/(r * (1.0 + (r/Rs)^2))

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