using Dates
using SpecialFunctions
using CSV 
using DataFrames

function get_Rv_in_kpc()

    # # Load King sphere in Henon units
    # if (!HAS_MULTI_MASS) 
    #     # Single mass
    #     namefile = path_dir*"data/IC/chc_king_ics_n_"*string(Npart)*".csv"
    # else 
    #     # Multi-mass (half of mass m1, half of mass m2=2*m1)
    #     namefile = path_dir*"data/IC/chc_king_ics_n_"*string(Npart)*"_multi_mass_1_2.csv"
    # end

    namefile = file_IC
    
    df = CSV.read(namefile, DataFrame, delim=',', header=false)

    datam = df[:, 3]
    datax = df[:, 6]
    datay = df[:, 7]
    dataz = df[:, 8]

    datar = sqrt.(datax .^2 + datay .^2 + dataz .^2)
    p = sortperm(datar)
    datar = datar[p]
    datam = datam[p]

    # Mtot = 1 HU 
    m_enc = 0.0
    Rh_HU = 0.0
    for i=1:Npart 
        r = datar[i]
        m_enc += datam[i]
        if (m_enc >= 0.5)
            Rh_HU = r 
            break
        end
    end

    # Rh_HU = datar[div(Npart, 2)]
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
const mass_avg = _Mtot/Npart
