using DelimitedFiles

include("../king/Constants.jl")

const params_MW2022 = readdlm("parameters/MW2022.txt", header=false)

# Nucleus
const a_nucleus_2022 = params_MW2022[1,1]/(R_HU_in_kpc)
const M_nucleus_2022 = params_MW2022[1,2]/(M_HU_in_Msun)

# Bulge 

const a_bulge_2022 = params_MW2022[2,1]/(R_HU_in_kpc)
const M_bulge_2022 = params_MW2022[2,2]/(M_HU_in_Msun)

# Disk
const a_disk1_2022 = params_MW2022[3,1]/(R_HU_in_kpc)
const b_disk1_2022 = params_MW2022[3,2]/(R_HU_in_kpc)
const M_disk1_2022 = params_MW2022[3,3]/(M_HU_in_Msun)

const a_disk2_2022 = params_MW2022[4,1]/(R_HU_in_kpc)
const b_disk2_2022 = params_MW2022[4,2]/(R_HU_in_kpc)
const M_disk2_2022 = params_MW2022[4,3]/(M_HU_in_Msun)

const a_disk3_2022 = params_MW2022[5,1]/(R_HU_in_kpc)
const b_disk3_2022 = params_MW2022[5,2]/(R_HU_in_kpc)
const M_disk3_2022 = params_MW2022[5,3]/(M_HU_in_Msun)

# Dark halo 
const rs_2022 = params_MW2022[6,1]/(R_HU_in_kpc)
const M_halo_2022 = params_MW2022[6,2]/(M_HU_in_Msun)



# Host potential (Henon units)

####################################################################################
# Nucleus and bulge 
# Hernquist potentials
####################################################################################


function psi_Hernquist(r::Float64, MH::Float64, aH::Float64)

    return -_G * MH/aH * 1/(1 + r/aH)

end

function dpsidr_Hernquist(r::Float64, MH::Float64, aH::Float64)

    return _G * MH/aH^2 * 1/(1 + r/aH)^2

end


# Nucleus

function psi_nucleus_2022(r::Float64)

    return psi_Hernquist(r, M_nucleus_2022, a_nucleus_2022)

end

function dpsidr_nucleus_2022(r::Float64)

    return dpsidr_Hernquist(r, M_nucleus_2022, a_nucleus_2022)

end

function acc_nucleus_2022(x::Float64, y::Float64, z::Float64)

    ar = -dpsidr_nucleus_2022(r)

    ax = ar * x/r
    ay = ar * y/r
    az = ar * z/r

    return ax, ay, az 
end


# Bulge 

function psi_bulge_2022(r::Float64)

    return psi_Hernquist(r, M_bulge_2022, a_bulge_2022)

end

function dpsidr_bulge_2022(r::Float64)

    return dpsidr_Hernquist(r, M_bulge_2022, a_bulge_2022)

end

function acc_bulge_2022(x::Float64, y::Float64, z::Float64)

    ar = -dpsidr_bulge_2022(r)

    ax = ar * x/r
    ay = ar * y/r
    az = ar * z/r

    return ax, ay, az 
end

####################################################################################
# Disk
# Miyamoto-Nagai potential components
####################################################################################


function psi_MN(R::Float64, z::Float64, MMN::Float64, aMN::Float64, bMN::Float64)

    return - _G * MMN/sqrt(R^2 + (sqrt(z^2 + bMN^2) + aMN)^2)

end

function dpsidR_MN(R::Float64, z::Float64, MMN::Float64, aMN::Float64, bMN::Float64)

    return _G * MMN * R/(R^2 + (sqrt(z^2 + bMN^2) + aMN)^2)^(3/2)

end

function dpsidz_MN(R::Float64, z::Float64, MMN::Float64, aMN::Float64, bMN::Float64)

    return _G * MMN * z * (sqrt(z^2 + bMN^2) + aMN)/(sqrt(z^2 + bMN)*(R^2 + (sqrt(z^2 + bMN^2) + aMN)^2)^(3/2))

end



function psi_disk_2022(R::Float64, z::Float64)

    psi1 = psi_MN(R, z, M_disk1_2022, a_disk1_2022, b_disk1_2022)
    psi2 = psi_MN(R, z, M_disk2_2022, a_disk2_2022, b_disk2_2022)
    psi3 = psi_MN(R, z, M_disk3_2022, a_disk3_2022, b_disk3_2022)

    return psi1 + psi2 + psi3 

end

function dpsidR_disk_2022(R::Float64, z::Float64)

    dpsidR1 = dpsidR_MN(R, z, M_disk1_2022, a_disk1_2022, b_disk1_2022)
    dpsidR2 = dpsidR_MN(R, z, M_disk2_2022, a_disk2_2022, b_disk2_2022)
    dpsidR3 = dpsidR_MN(R, z, M_disk3_2022, a_disk3_2022, b_disk3_2022)

    return dpsidR1 + dpsidR2 + dpsidR3 

end

function dpsidz_disk_2022(R::Float64, z::Float64)

    dpsidz1 = dpsidz_MN(R, z, M_disk1_2022, a_disk1_2022, b_disk1_2022)
    dpsidz2 = dpsidz_MN(R, z, M_disk2_2022, a_disk2_2022, b_disk2_2022)
    dpsidz3 = dpsidz_MN(R, z, M_disk3_2022, a_disk3_2022, b_disk3_2022)

    return dpsidz1 + dpsidz2 + dpsidz3 

end

function acc_disk_2022(x::Float64, y::Float64, z::Float64)

    R = sqrt(x*x + y*y)

    aR = -dpsidR_disk_2022(R, z)
    az = -dpsidz_disk_2022(R, z)

    ax = aR * x/R 
    ay = aR * y/R 

    return ax, ay, az

end

####################################################################################
# Halo
# NFW potential
####################################################################################

function psi_halo_2022(r::Float64)

    return -_G * M_halo_2022/r * log(1 + r/rs_2022)

end

function dpsidr_halo_2022(r::Float64)

    g_x = log(1 + r/rs_2022) - (r/rs_2022)/(1 + r/rs_2022)
    Menc = M_halo_2022 * g_x 

    return _G * Menc/r^2

end

function acc_halo_2022(r::Float64)

    ar = -dpsidr_halo_2022(r)

    ax = ar * x/r
    ay = ar * y/r
    az = ar * z/r

    return ax, ay, az 
end


####################################################################################
# Total potential
####################################################################################

function psi_host_2022(R::Float64, z::Float64)
    
    r = sqrt(R^2 + z^2)

    psin = psi_nucleus_2022(r)
    psib = psi_bulge_2022(r)
    psid = psi_disk_2022(R, z)
    psih = psi_halo_2022(r)

    psi =  psin + psib + psid + psih 

    return psi

end


####################################################################################
# Total acceleration
####################################################################################

function acc_host_2022(x::Float64, y::Float64, z::Float64)

    axn, ayn, azn = acc_nucleus_2022(x, y, z)
    axb, ayb, azb = acc_bulge_2022(x, y, z)
    axd, ayd, azd = acc_disk_2022(x, y, z)
    axh, ayh, azh = acc_halo_2022(x, y, z)

    ax = axn + axb + axd + axh
    ay = ayn + ayb + ayd + ayh
    az = azn + azb + azd + azh

    return ax, ay, az 

end


####################################################################################
# Circular velocity
####################################################################################

# z=0
function circular_velocity_2022(R::Float64)

    dpsidRn = dpsidr_nucleus_2022(R)
    dpsidRb = dpsidr_bulge_2022(R)
    dpsidRd = dpsidR_disk_2022(R, 0.0)
    dpsidRh = dpsidr_halo_2022(R)

    dpsidR = dpsidRn + dpsidRb + dpsidRd + dpsidRh

    return sqrt(R*dpsidR)

end
