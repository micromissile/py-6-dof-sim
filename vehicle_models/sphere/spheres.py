import math

def Musketball50cal():

    # Inches to meters conversion
    in2m = 0.0254

    # Name of sphere 
    vehicle_name = "50 Cal Lead Ball"
    
    # Short name for file names
    short_name = "50cal"
    
    # Density of sphere
    rho_lead_kgpm3      = 11300
    
    # Radius of sphere
    r_sphere_in  = 0.495
    r_sphere_m  = r_sphere_in*in2m
    
    # Reference wing span
    b_m = r_sphere_m
    
    # Reference wing chord
    c_m = r_sphere_m

    # Approximate drag coefficient for subsonic regimes
    CD_approx = 0.5
    
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
    vol_sphere_m3, m_sphere_kg, J_sphere_kgm2, A_ref_m2 = CalcSphereProps(r_sphere_m, rho_lead_kgpm3)

    # Terminal velocity approximation
    Vterm_mps = math.sqrt((2*m_sphere_kg*9.81)/(1.2*CD_approx*A_ref_m2))

    vmod = {"V_name"     : vehicle_name, \
            "m_kg"       : m_sphere_kg, \
            "Jxz_b_kgm2" : 0, \
            "Jxx_b_kgm2" : J_sphere_kgm2, \
            "Jyy_b_kgm2" : J_sphere_kgm2, \
            "Jzz_b_kgm2" : J_sphere_kgm2, \
            "r_sphere_m" : r_sphere_m, \
            "m_sphere_kg": m_sphere_kg, \
            "CD_approx"  : CD_approx, \
            "Clp"        : Clp, \
            "Clr"        : Clr, \
            "Cmq"        : Cmq, \
            "Cnp"        : Cnp, \
            "Cnr"        : Cnr, \
            "A_ref_m2"   : A_ref_m2, \
            "Vterm_mps"  : Vterm_mps, \
            "b_m" : b_m, \
            "c_m" : c_m, \
            "short_name" : short_name} 
    
    return vmod

def Carronade12lb():

    # Inches to meters conversion
    in2m = 0.0254

    # Name of sphere 
    vehicle_name = "Carronade 12 lb (5.4 kg) Cannonball"
    
    # Short name for file names
    short_name = "12lb"
    
    # Density of sphere
    rho_castiron_kgpm3 = 7000
    
    # Radius of sphere
    r_sphere_in  = 4.40
    r_sphere_m  = r_sphere_in*in2m
    
    # Reference wing span
    b_m = r_sphere_m
    
    # Reference wing chord
    c_m = r_sphere_m

    # Approximate drag coefficient for subsonic regimes
    CD_approx = 0.5
    
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
    vol_sphere_m3, m_sphere_kg, J_sphere_kgm2, A_ref_m2 = CalcSphereProps(r_sphere_m, rho_castiron_kgpm3)

    # Terminal velocity approximation
    Vterm_mps = math.sqrt((2*m_sphere_kg*9.81)/(1.2*CD_approx*A_ref_m2))

    vmod = {"V_name"     : vehicle_name, \
            "m_kg"       : m_sphere_kg, \
            "Jxz_b_kgm2" : 0, \
            "Jxx_b_kgm2" : J_sphere_kgm2, \
            "Jyy_b_kgm2" : J_sphere_kgm2, \
            "Jzz_b_kgm2" : J_sphere_kgm2, \
            "r_sphere_m" : r_sphere_m, \
            "m_sphere_kg": m_sphere_kg, \
            "CD_approx"  : CD_approx, \
            "Clp"        : Clp, \
            "Clr"        : Clr, \
            "Cmq"        : Cmq, \
            "Cnp"        : Cnp, \
            "Cnr"        : Cnr, \
            "A_ref_m2"   : A_ref_m2, \
            "Vterm_mps"  : Vterm_mps, \
            "b_m" : b_m, \
            "c_m" : c_m, \
            "short_name" : short_name} 

    return vmod

def Blueberry():

    # Inches to meters conversion
    in2m = 0.0254

    # Name of sphere 
    vehicle_name = "A Blueberry"
    
    # Short name for file names
    short_name = "Blueberry"
    
    # Density of a single blueberry
    rho_blueberry_kgpm3 = 786
    
    # Radius of sphere
    r_sphere_in  = 0.3
    r_sphere_m  = r_sphere_in*in2m
    
    # Reference wing span
    b_m = r_sphere_m
    
    # Reference wing chord
    c_m = r_sphere_m

    # Approximate drag coefficient for subsonic regimes
    CD_approx = 0.5
    
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
    vol_sphere_m3, m_sphere_kg, J_sphere_kgm2, A_ref_m2 = CalcSphereProps(r_sphere_m, rho_blueberry_kgpm3)

    # Terminal velocity approximation
    Vterm_mps = math.sqrt((2*m_sphere_kg*9.81)/(1.2*CD_approx*A_ref_m2))

    vmod = {"V_name"     : vehicle_name, \
            "m_kg"       : m_sphere_kg, \
            "Jxz_b_kgm2" : 0, \
            "Jxx_b_kgm2" : J_sphere_kgm2, \
            "Jyy_b_kgm2" : J_sphere_kgm2, \
            "Jzz_b_kgm2" : J_sphere_kgm2, \
            "r_sphere_m" : r_sphere_m, \
            "m_sphere_kg": m_sphere_kg, \
            "CD_approx"  : CD_approx, \
            "Clp"        : Clp, \
            "Clr"        : Clr, \
            "Cmq"        : Cmq, \
            "Cnp"        : Cnp, \
            "Cnr"        : Cnr, \
            "A_ref_m2"   : A_ref_m2, \
            "Vterm_mps"  : Vterm_mps, \
            "b_m" : b_m, \
            "c_m" : c_m, \
            "short_name" : short_name} 

    return vmod

def BowlingBall():

    # Inches to meters conversion
    in2m = 0.0254

    # Name of sphere 
    vehicle_name = "Bowling Ball"
    
    # Short name for file names
    short_name = "Bowling"
    
    # Density of sphere
    rho_bowlingball_kgpm3 = 1500
    
    # Radius of sphere
    r_sphere_in  = 4.40
    r_sphere_m  = r_sphere_in*in2m
    
    # Reference wing span
    b_m = r_sphere_m
    
    # Reference wing chord
    c_m = r_sphere_m

    # Approximate drag coefficient for subsonic regimes
    CD_approx = 0.5
    
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
    vol_sphere_m3, m_sphere_kg, J_sphere_kgm2, A_ref_m2 = CalcSphereProps(r_sphere_m, rho_bowlingball_kgpm3)

    # Terminal velocity approximation
    Vterm_mps = math.sqrt((2*m_sphere_kg*9.81)/(1.2*CD_approx*A_ref_m2))

    vmod = {"V_name"     : vehicle_name, \
            "m_kg"       : m_sphere_kg, \
            "Jxz_b_kgm2" : 0, \
            "Jxx_b_kgm2" : J_sphere_kgm2, \
            "Jyy_b_kgm2" : J_sphere_kgm2, \
            "Jzz_b_kgm2" : J_sphere_kgm2, \
            "r_sphere_m" : r_sphere_m, \
            "m_sphere_kg": m_sphere_kg, \
            "CD_approx"  : CD_approx, \
            "Clp"        : Clp, \
            "Clr"        : Clr, \
            "Cmq"        : Cmq, \
            "Cnp"        : Cnp, \
            "Cnr"        : Cnr, \
            "A_ref_m2"   : A_ref_m2, \
            "Vterm_mps"  : Vterm_mps, \
            "b_m" : b_m, \
            "c_m" : c_m, \
            "short_name" : short_name} 

    return vmod

def TsarCannonball():

    # Inches to meters conversion
    in2m = 0.0254

    # Name of sphere 
    vehicle_name = "Tsar Cannonball"
    
    # Short name for file names
    short_name = "Tsar"
    
    # Density of sphere
    rho_castiron_kgpm3 = 7000
    
    # Radius of sphere
    r_sphere_in  = 35
    r_sphere_m  = r_sphere_in*in2m
    
    # Reference wing span
    b_m = r_sphere_m
    
    # Reference wing chord
    c_m = r_sphere_m

    # Approximate drag coefficient for subsonic regimes
    CD_approx = 0.5
    
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
    vol_sphere_m3, m_sphere_kg, J_sphere_kgm2, A_ref_m2 = CalcSphereProps(r_sphere_m, rho_castiron_kgpm3)

    # Terminal velocity approximation
    Vterm_mps = math.sqrt((2*m_sphere_kg*9.81)/(1.2*CD_approx*A_ref_m2))

    vmod = {"V_name"     : vehicle_name, \
            "m_kg"       : m_sphere_kg, \
            "Jxz_b_kgm2" : 0, \
            "Jxx_b_kgm2" : J_sphere_kgm2, \
            "Jyy_b_kgm2" : J_sphere_kgm2, \
            "Jzz_b_kgm2" : J_sphere_kgm2, \
            "r_sphere_m" : r_sphere_m, \
            "m_sphere_kg": m_sphere_kg, \
            "CD_approx"  : CD_approx, \
            "Clp"        : Clp, \
            "Clr"        : Clr, \
            "Cmq"        : Cmq, \
            "Cnp"        : Cnp, \
            "Cnr"        : Cnr, \
            "A_ref_m2"   : A_ref_m2, \
            "Vterm_mps"  : Vterm_mps, \
            "b_m" : b_m, \
            "c_m" : c_m, \
            "short_name" : short_name}  

    return vmod

def NASA_Atmos01_Sphere():

    # Inches to meters conversion
    in2m = 0.0254

    # Name of sphere 
    vehicle_name = "NASA Atmos01 1-Slug Cannonball"
    
    # Short name for file names
    short_name = "Atmos01"
    
    # Density of sphere
    rho_NASA_Sphere_kgpm3 = 7868.36
    
    # Radius of sphere
    r_sphere_in  = 3
    r_sphere_m  = r_sphere_in*in2m
    
    # Reference wing span
    b_m = r_sphere_m
    
    # Reference wing chord
    c_m = r_sphere_m

    # Approximate drag coefficient for subsonic regimes
    CD_approx = 0.5
    
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
    vol_sphere_m3, m_sphere_kg, J_sphere_kgm2, A_ref_m2 = CalcSphereProps(r_sphere_m, rho_NASA_Sphere_kgpm3)

    # Terminal velocity approximation
    Vterm_mps = math.sqrt((2*m_sphere_kg*9.81)/(1.2*CD_approx*A_ref_m2))

    vmod = {"V_name"     : vehicle_name, \
            "m_kg"       : m_sphere_kg, \
            "Jxz_b_kgm2" : 0, \
            "Jxx_b_kgm2" : J_sphere_kgm2, \
            "Jyy_b_kgm2" : J_sphere_kgm2, \
            "Jzz_b_kgm2" : J_sphere_kgm2, \
            "r_sphere_m" : r_sphere_m, \
            "m_sphere_kg": m_sphere_kg, \
            "CD_approx"  : CD_approx, \
            "Clp"        : Clp, \
            "Clr"        : Clr, \
            "Cmq"        : Cmq, \
            "Cnp"        : Cnp, \
            "Cnr"        : Cnr, \
            "A_ref_m2"   : A_ref_m2, \
            "Vterm_mps"  : Vterm_mps, \
            "b_m" : b_m, \
            "c_m" : c_m, \
            "short_name" : short_name} 

    return vmod

def CalcSphereProps(r_sphere_m, rho_sphere_kgpm3):
    """ FUNCTION CalcSphereProps calculates the mass 
    of a sphere for flight simulation"""

    vol_sphere_m3 = 4/3*math.pi*r_sphere_m**3
    m_sphere_kg = rho_sphere_kgpm3*vol_sphere_m3
    J_sphere_kgm2   = 0.4*m_sphere_kg*r_sphere_m**2
    A_ref_m2         = math.pi*r_sphere_m**2

    return vol_sphere_m3, m_sphere_kg, J_sphere_kgm2, A_ref_m2
    