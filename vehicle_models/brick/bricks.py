import math

'''A simple, made-up constant mass matrix for a US face brick (8" x 4" x 2.25") for use in the
2012 NASA NESC EOM verification data project.'''

def NASA_Atmos02_Brick():

        # Inches to meters conversion
        in2m = 0.0254

        # Slug to kg conversion
        slug2kg = 14.5939
        kg2slug = 1/slug2kg

        # feet to meters
        ft2m = 0.304878

        # Name of sphere 
        vehicle_name = "Tumbling Brick (No Damping or Drag)"
        
        # Short name
        short_name = "Atmos02"

        # Mass of brick
        m_brick_slug = 0.1554048
        m_brick_kg = kg2slug*m_brick_slug

        # Moments and products of inertia
        Jxx_slugft2 = 0.00189422
        Jxx_kgm2 = slug2kg*(ft2m**2)*Jxx_slugft2
        
        Jyy_slugft2 = 0.00621102
        Jyy_kgm2 = slug2kg*(ft2m**2)*Jyy_slugft2
        
        Jzz_slugft2 = 0.00719467
        Jzz_kgm2 = slug2kg*(ft2m**2)*Jzz_slugft2
        
        Jzx_slugft2 = 0
        Jzx_kgm2 = slug2kg*(ft2m**2)*Jzx_slugft2

        # Body position of center of mass with respect to the geometric center of the brick
        dxcg_ft = 0
        dycg_ft = 0
        dzcg_ft = 0
        
        # Dimensions of brick
        length_brick_m = 8*in2m
        width_brick_m  = 4*in2m
        
        # Reference area
        A_ref_m2 = length_brick_m*width_brick_m
        
        # Reference wing span
        b_m = 0.33333*ft2m
        
        # Reference wing chord
        c_m = 0.66667*ft2m
        
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

        vmod = {"V_name"     : vehicle_name, \
                "m_kg"       : m_brick_kg, \
                "Jxz_b_kgm2" : Jzx_kgm2, \
                "Jxx_b_kgm2" : Jxx_kgm2, \
                "Jyy_b_kgm2" : Jyy_kgm2, \
                "Jzz_b_kgm2" : Jzz_kgm2, \
                "A_ref_m2" : A_ref_m2, \
                "Clp" : Clp, \
                "Clr" : Clr, \
                "Cmq" : Cmq, \
                "Cnp" : Cnp, \
                "Cnr" : Cnr, \
                "b_m" : b_m, \
                "c_m" : c_m, \
                "short_name" : short_name} 

        return vmod

def NASA_Atmos03_Brick():

        # Inches to meters conversion
        in2m = 0.0254

        # Slug to kg conversion
        slug2kg = 14.5939
        kg2slug = 1/slug2kg

        # feet to meters
        ft2m = 0.304878

        # Name of sphere 
        vehicle_name = "Tumbling Brick (with Aerodynamic Damping)"
        
        # Short name
        short_name = "Atmos03"

        # Mass of brick
        m_brick_slug = 0.1554048
        m_brick_kg = kg2slug*m_brick_slug

        # Moments and products of inertia
        Jxx_slugft2 = 0.00189422
        Jxx_kgm2 = slug2kg*(ft2m**2)*Jxx_slugft2
        
        Jyy_slugft2 = 0.00621102
        Jyy_kgm2 = slug2kg*(ft2m**2)*Jyy_slugft2
        
        Jzz_slugft2 = 0.00719467
        Jzz_kgm2 = slug2kg*(ft2m**2)*Jzz_slugft2
        
        Jzx_slugft2 = 0
        Jzx_kgm2 = slug2kg*(ft2m**2)*Jzx_slugft2

        # Body position of center of mass with respect to the geometric center of the brick
        dxcg_ft = 0
        dycg_ft = 0
        dzcg_ft = 0
        
        # Dimensions of brick
        length_brick_m = 8*in2m
        width_brick_m  = 4*in2m
        
        # Reference area
        A_ref_m2 = length_brick_m*width_brick_m
        
        # Reference wing span
        b_m = 0.33333*ft2m
        
        # Reference wing chord
        c_m = 0.66667*ft2m
        
        # Roll damping from roll rate
        Clp = -1.0
        
        # Roll damping from yaw rate
        Clr = 0.0
        
        # Pitch damping from pitch rate
        Cmq = -1.0
        
        # Yaw damping from roll rate
        Cnp = 0.0
        
        # Yaw damping from yaw rate
        Cnr = -1.0

        vmod = {"V_name"     : vehicle_name, \
                "m_kg"       : m_brick_kg, \
                "Jxz_b_kgm2" : Jzx_kgm2, \
                "Jxx_b_kgm2" : Jxx_kgm2, \
                "Jyy_b_kgm2" : Jyy_kgm2, \
                "Jzz_b_kgm2" : Jzz_kgm2, \
                "A_ref_m2" : A_ref_m2, \
                "Clp" : Clp, \
                "Clr" : Clr, \
                "Cmq" : Cmq, \
                "Cnp" : Cnp, \
                "Cnr" : Cnr, \
                "b_m" : b_m, \
                "c_m" : c_m, \
                "short_name" : short_name}  

        return vmod

