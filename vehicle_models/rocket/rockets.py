import math
from pathlib import Path

import pandas as pd

current_dir = Path(__file__).parent

OPENROCKET_TIME_FIELD_NAME = "# Time (s)"
OPENROCKET_THRUST_FIELD_NAME = "Thrust (N)"

def RocketDev():

    mmtom = lambda x: x * 0.001
    gtokg = lambda x: x * 0.001

    # Name of sphere 
    vehicle_name = "RocketDev"
    
    # Short name for file names
    short_name = "RocketDev"

    m_sphere_kg = gtokg(518)
    
    # Density of a single blueberry
    # rho_blueberry_kgpm3 = 786
    
    # Radius of sphere
    r_sphere_m  = mmtom(62.4/2)
    
    # Reference wing span
    b_m = r_sphere_m
    
    # Reference wing chord
    c_m = r_sphere_m

    # Approximate drag coefficient for subsonic regimes
    CD_approx = 0.0
    
    # Roll damping from roll rate
    Clp = 0.0
    
    # Roll damping from yaw rate
    Clr = 0.0
    
    # Pitch damping from pitch rate
    Cmq = 0.0
    
    # Yaw damping from roll rate
    Cnp = 0.0
    
    # Yaw damping from yaw rate
    Cnr = 0.0

    # Compute sphere properties
    # vol_sphere_m3, m_sphere_kg, J_sphere_kgm2, A_ref_m2 = CalcSphereProps(r_sphere_m, rho_blueberry_kgpm3)
    J_sphere_kgm2 = 1.

    A_ref_m2 = math.pi * r_sphere_m ** 2
    Vterm_mps = 0.0

    # Terminal velocity approximation
    # Vterm_mps = math.sqrt((2*m_sphere_kg*9.81)/(1.2*CD_approx*A_ref_m2))o

    # OpenRocket simulation data
    df = pd.read_csv(current_dir / "openrocket_sim_data_3_motors.csv")

    vmod = {
        "V_name"     : vehicle_name,
        "m_kg"       : m_sphere_kg,
        "Jxz_b_kgm2" : 0,
        "Jxx_b_kgm2" : J_sphere_kgm2,
        "Jyy_b_kgm2" : J_sphere_kgm2,
        "Jzz_b_kgm2" : J_sphere_kgm2,
        "r_sphere_m" : r_sphere_m,
        "m_sphere_kg": m_sphere_kg,
        "CD_approx"  : CD_approx,
        "Clp"        : Clp,
        "Clr"        : Clr,
        "Cmq"        : Cmq,
        "Cnp"        : Cnp,
        "Cnr"        : Cnr,
        "A_ref_m2"   : A_ref_m2,
        "Vterm_mps"  : Vterm_mps,
        "b_m"        : b_m,
        "c_m"        : c_m,
        "short_name" : short_name,
        "motor_thrust_time_axis": df[OPENROCKET_TIME_FIELD_NAME].values,
        "motor_thrust_N": df[OPENROCKET_THRUST_FIELD_NAME].values
    } 

    return vmod

# def CalcSphereProps(r_sphere_m, rho_sphere_kgpm3):
#     """ FUNCTION CalcSphereProps calculates the mass 
#     of a sphere for flight simulation"""

#     vol_sphere_m3 = 4/3*math.pi*r_sphere_m**3
#     m_sphere_kg = rho_sphere_kgpm3*vol_sphere_m3
#     J_sphere_kgm2   = 0.4*m_sphere_kg*r_sphere_m**2
#     A_ref_m2         = math.pi*r_sphere_m**2

#     return vol_sphere_m3, m_sphere_kg, J_sphere_kgm2, A_ref_m2
    