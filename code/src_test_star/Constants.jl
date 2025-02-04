using Dates
using SpecialFunctions

# Path to src/

const path_to_src = @__DIR__
const path_dir = path_to_src * "/../../"

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
const M_HU_in_Msun = Mtot_Msun # Value of 1 HU mass in solar masses
const R_HU_in_kpc = Rv_kpc # Value of 1 HU length in kpc
const G_in_kpc_Mpc_Myr = 4.49e-12
const T_HU_in_Myr = sqrt(R_HU_in_kpc^3/(G_in_kpc_Mpc_Myr*M_HU_in_Msun)) # Myr # T = sqrt(Rv^3/(G*M)) = 4.22 

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