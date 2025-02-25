using Dates
using SpecialFunctions

# Path to src/

const path_to_src = @__DIR__
const path_dir = path_to_src * "/../../../"

# Current time for file saving 
const date = now()
const sdate = Dates.format(date, "yyyy-mm-dd_HH-MM-SS")
const run = Dates.value(date)
const srun = string(run)

# Henon units
const _Mtot = 1.0
const R_vir = 1.0
const _G = 1.0

# Conversion HU to astrophysical units
# >>> l=(u.m).to(u.kpc)
# >>> l
# 3.2407792894443654e-20
# >>> m=(u.kg).to(u.M_sun)
# >>> m
# 5.029144215870041e-31
# >>> t=(u.s).to(u.Myr)
# >>> t
# 3.168808781402895e-14
# >>> l**3/(m*t**2)
# 0.0674003588611473
# >>> G * l**3/(m*t**2)
# <Quantity 4.49850215e-12

# Conversion HU to astrophysical units
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_MSun_Myr = 4.49851e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_MSun_Myr*M_HU_in_Msun)) # Myr

println("1 time HU = ", T_HU_in_Myr, " Myr")

# Conversion kpc/Myr to km/s (using astropy)
# x kpc/Myr = x kpc/km s/Myr km/s = y km/s ; y = x kpc/km s/Myr
# >>> l=1*u.kpc
# >>> l.to(u.km)
# <Quantity 3.08567758e+16 km>
# >>> t=u.Myr
# >>> t.to(u.s)
# <Quantity 3.15576e+13 s>
# lkm=l.to(u.km)
# ts=t.to(u.s)
# >>> lkm/ts
# <Quantity 977.79222168 km / s>

const V_HU_in_kpc_Myr = sqrt((G_in_kpc_MSun_Myr*M_HU_in_Msun)/R_HU_in_kpc)
const V_HU_in_km_s = V_HU_in_kpc_Myr * 977.79222168
# x kpc/Myr = x kpc/km s/Myr km/s = y km/s ; y = x kpc/km s/Myr

# Cluster potential
const d_host = d_kpc/R_HU_in_kpc # Distance to host's centre


# Dark halo
const Mvir = Mvir_Msun/(M_HU_in_Msun)
const Rs = Rs_kpc/(R_HU_in_kpc)
const g_c = log(1+c) - c/(1.0+c)

# Bulge
const M_bulge = M_bulge_Msun/(M_HU_in_Msun)
const s_bulge = 0.5*(3.0-alpha_bulge)
const rc_bulge = rc_bulge_kpc/(R_HU_in_kpc)
const gamma_s = gamma(s_bulge)

# Disk
const M_disk = M_disk_Msun/(M_HU_in_Msun)
const a_disk = a_disk_kpc/(R_HU_in_kpc)
const b_disk = b_disk_kpc/(R_HU_in_kpc)