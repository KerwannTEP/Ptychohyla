using SpecialFunctions
using DelimitedFiles

include("../king/Constants.jl")

const params_MW2014 = readdlm("parameters/MW2014.txt", header=false)


# Dark halo
const Mvir_2014 = params_MW2014[1,3]/(M_HU_in_Msun)
const Rs_2014 = params_MW2014[1,2]/(R_HU_in_kpc)
const g_c_2014 = log(1+params_MW2014[1,1]) - c/(1.0+params_MW2014[1,1])

# Bulge
const M_bulge_2014 = params_MW2014[2,3]/(M_HU_in_Msun)
const alpha_bulge_2014 = params_MW2014[2,1]
const rc_bulge_2014 = params_MW2014[2,2]/(R_HU_in_kpc)
const s_bulge_2014 = 0.5*(3.0-alpha_bulge_2014)
const gamma_s_2014 = gamma(s_bulge_2014)

# Disk
const M_disk_2014 = params_MW2014[3,3]/(M_HU_in_Msun)
const a_disk_2014 = params_MW2014[3,1]/(R_HU_in_kpc)
const b_disk_2014 = params_MW2014[3,2]/(R_HU_in_kpc)



# Host potential (Henon units)

####################################################################################
# Bulge
# CHC II appendix C1
####################################################################################

# CAREFUL: gamma(s, x) = int_x^infty ... is the upper incomplete gamma function
# Enclosed mass of the bulge

function gamma_low_sx(xr)

    return gamma_s_2014 - gamma(s_bulge, xr)

end

function Menc_bulge_2014(r::Float64)

    xr = (r/rc_bulge_2014)^2
    return M_bulge_2014 * gamma_low_sx(xr)/gamma_s_2014
end

function acc_bulge_2014(x::Float64, y::Float64, z::Float64)

    r = sqrt(x*x + y*y + z*z)
    xr = (r/rc_bulge_2014)^2
    gamma_s_x = gamma_low_sx(xr)

    Menc = M_bulge_2014 * gamma_s_x/gamma_s_2014
    ar = - _G * Menc/r^2 

    ax = ar * x/r
    ay = ar * y/r
    az = ar * z/r

    return ax, ay, az 
end

function psi_bulge_2014(r::Float64)

    xr = (r/rc_bulge_2014)^2

    # Inner potential
    gamma_lower = gamma_s_2014 - gamma(s_bulge, xr)
    psi_inner = -_G*M_bulge_2014/r * gamma_lower/gamma_s_2014

    # Outer potential 
    gamma_upper = gamma(1.0-alpha_bulg_2014/2.0, xr)
    psi_outer = -_G*M_bulge/rc_bulge * gamma_upper/gamma_s
     
    return psi_inner + psi_outer

end

function dpsidr_bulge_2014(r::Float64)

    xr = (r/rc_bulge_2014)^2
    gamma_s_x = gamma_low_sx(xr)
    return _G * M_bulge_2014/r^2 * gamma_s_x/gamma_s_2014
end

function rho_bulge_2014(r::Float64)

    rho0 = M_bulge_2014/(2.0*pi*rc_bulge_2014^3*gamma_s_2014)
    return rho0 * (rc_bulge_2014/r)^(alpha_bulge_2014) * exp(-(r/rc_bulge_2014)^2)

end

####################################################################################
# Disk
# Miyamoto-Nagai disk
# CHC II appendix C1
####################################################################################

function acc_disk_2014(x::Float64, y::Float64, z::Float64)

    R = sqrt(x*x + y*y)

    aR = - _G * M_disk_2014 * R/(R^2 + (sqrt(z^2 + b_disk_2014^2)+a_disk_2014)^2)^(3/2)
    az = - _G * M_disk_2014 * z * (sqrt(z^2 + b_disk_2014^2) + a_disk_2014)/(sqrt(z^2+b_disk_2014^2)*(R^2 + (sqrt(z^2 + b_disk_2014^2)+a_disk_2014)^2)^(3/2))

    ax = aR * x/R 
    ay = aR * y/R 

    return ax, ay, az

end

function psi_disk_2014(R::Float64, z::Float64)

    return -_G * M_disk_2014/sqrt(R^2 + (sqrt(z^2 + b_disk_2014^2) + a_disk_2014)^2)

end

function dpsidR_disk_2014(R::Float64, z::Float64)

    return _G * M_disk_2014 * R/(R^2 + (sqrt(z^2 + b_disk_2014^2)+a_disk_2014)^2)^(3/2)
end

# https://galaxiesbook.org/chapters/II-01.-Flattened-Mass-Distributions.html
# eq (8.20)
function rho_disk_2014(R::Float64, z::Float64)

    num = a_disk_2014*R^2 + (3.0*sqrt(z^2+b_disk_2014^2)+a_disk_2014) * (sqrt(z^2+b_disk_2014^2)+a_disk_2014)^2
    den = (R^2 + (sqrt(z^2+b_disk_2014^2)+a_disk_2014)^2)^(5/2) * (z^2+b_disk_2014^2)^(3/2)

    return b_disk_2014^2*M_disk_2014/(4*pi) * num/den

end

####################################################################################
# Dark halo (NFW)
# https://en.wikipedia.org/wiki/Navarro%E2%80%93Frenk%E2%80%93White_profile
####################################################################################

function _g(x::Float64)

    return log(1.0+x) - x/(1+x)
end

function Menc_halo_2014(r::Float64)

    return Mvir * _g(r/Rs_2014)/g_c_2014
    
end

function acc_halo_2014(x::Float64, y::Float64, z::Float64)

    r = sqrt(x*x + y*y + z*z)

    Menc = Menc_halo_2014(r)
    ar = - _G * Menc/r^2

    ax = ar * x/r
    ay = ar * y/r
    az = ar * z/r

    return ax, ay, az 

end

function psi_halo_2014(r::Float64)

    x = r/Rs_2014
    cst = _G*Mvir_2014/Rs_2014 * 1.0/g_c_2014

    return -cst/x * log(1+x)

end

function dpsidr_halo_2014(r::Float64)

    Menc = Menc_halo_2014(r)
    return _G * Menc/r^2
end

# https://galaxiesbook.org/chapters/I-01.-Potential-Theory-and-Spherical-Mass-Distributions.html
# NFW 
# alpha=1, beta=3
function rho_halo_2014(r::Float64)

    #4*pi*rho0*Rs^3 = Mvir/g_c

    rho0 = Mvir_2014 / (g_c_2014 * 4*pi*Rs_2014^3)
    return rho0 * Rs_2014/(r * (1.0 + (r/Rs_2014)^2))

end

####################################################################################
# Total density
####################################################################################

function rho_host_2014(R::Float64, z::Float64)

    rho_b = rho_bulge_2014(sqrt(R^2 + z^2))
    rho_d = rho_disk_2014(R, z)
    rho_h = rho_halo_2014(sqrt(R^2 + z^2))
    
    rho = rho_h + rho_d + rho_b

    return rho

end

####################################################################################
# Total potential
####################################################################################

function psi_host_2014(R::Float64, z::Float64)

    psi_b = psi_bulge_2014(sqrt(R^2 + z^2))
    psi_d = psi_disk_2014(R, z)
    psi_h = psi_halo_2014(sqrt(R^2 + z^2))
    
    psi = psi_h + psi_d + psi_b

    return psi

end

####################################################################################
# Total acceleration
####################################################################################

function acc_host_2014(x::Float64, y::Float64, z::Float64)

    ax_b, ay_b, az_b = acc_bulge_2014(x, y, z)
    ax_d, ay_d, az_d = acc_disk_2014(x, y, z)
    ax_h, ay_h, az_h = acc_halo_2014(x, y, z)

    ax = ax_h + ax_d + ax_b
    ay = ay_h + ay_d + ay_b
    az = az_h + az_d + az_b

    return ax, ay, az 

end

####################################################################################
# Circular velocity
####################################################################################

# z=0
function circular_velocity_2014(R::Float64)

    dpsidR_b = dpsidr_bulge(R)
    dpsidR_d = dpsidR_disk(R, 0.0)
    dpsidR_h = dpsidr_halo(R)

    dpsidR = dpsidR_h + dpsidR_d + dpsidR_b

    return sqrt(R*dpsidR)

end