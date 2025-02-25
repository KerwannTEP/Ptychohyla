using Dates
using SpecialFunctions
using CSV 
using DataFrames

function get_Rv_in_kpc()

    # Load King sphere in Henon units
    namefile = path_dir*"data/IC/chc_king_ics_n_"*string(Npart)*".csv"
    df = CSV.read(namefile, DataFrame, delim=',', header=false)


    datax = df[:, 6]
    datay = df[:, 7]
    dataz = df[:, 8]

    datar = sqrt.(datax .^2 + datay .^2 + dataz .^2)
    datar = sort(datar)

    Rh_HU = datar[div(Npart, 2)]
    Rt_HU = datar[Npart]

    # Rv_in_kpc/Rh_in_kpc = Rv_in_HU/Rh_in_HU
    # => Rv_in_kpc = Rh_in_kpc/Rh_in_HU

    Rv_in_kpc = Rh_kpc/Rh_HU
    Rt_in_kpc = Rt_HU * Rv_in_kpc

    println("Rv [kpc] = ", Rv_in_kpc)
    println("Rt [kpc] = ", Rt_in_kpc)

    return Rv_in_kpc

end

# Path to src/

const path_to_src = @__DIR__
const path_dir = path_to_src * "/../../../"

# Current time for file saving 
const date = now()
const sdate = Dates.format(date, "yyyy-mm-dd_HH-MM-SS")

# Ternary operator 
# https://stackoverflow.com/questions/39790031/does-julia-have-a-ternary-conditional-operator
# https://en.wikibooks.org/wiki/Introducing_Julia/Controlling_the_flow#Ternary_expressions
const run = (id_default >= 0) ? id_default : Dates.value(date)

const srun = string(run)

# Henon units
const _Mtot = 1.0
const R_vir = 1.0
const _G = 1.0

# Conversion HU to astrophysical units (using astropy)
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
const R_HU_in_kpc = get_Rv_in_kpc() # Value of 1 HU length in kpc
const G_in_kpc_MSun_Myr = 4.49850215e-12
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


# Cluster potential
const d_host = d_kpc/R_HU_in_kpc # Distance to host's centre
const mass = _Mtot/Npart

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