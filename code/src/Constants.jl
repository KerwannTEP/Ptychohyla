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
const mass = _Mtot/Npart


# Host potential
const g_1 = log(2) - 1/2
const g_c = log(1+c) - c/(1+c)

const M_DH_200 =M_DH_200_Msun/(M_HU_in_Msun)
const Rs = Rs_kpc/(R_HU_in_kpc)
const Menc = M_DH_200*g_1/g_c # Enclosed mass within a sphere of radius Rs
const rho0_host = Menc/(4*pi*Rs^3 * g_1)

