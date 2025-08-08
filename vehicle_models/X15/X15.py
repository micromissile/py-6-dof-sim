import math
import numpy as np

def X15_at_Mach_3():
        
        # Rocket motor burned out
        # Speed brake closed
        # Data taken from WT line on plot
        # Lower rudder on

        # Conversions
        ft2m = 0.3048

        # Name of vehicle (used in plot titles)
        vehicle_name = "X15 Space Plane at Mach 3.4"

        # Short name (used in file names)
        short_name = "X15_at_Mach_3"
        
        # Mass
        m_kg = 11765
        m_kg = 14000

        # Reference wing span
        b_m = 22.36*ft2m

        # Reference mean aerodynamic chord
        c_m = 10.27*ft2m
        
        # Reference area
        A_ref_m2 = 18.6
        
        # Aerodynamic coefficients will be set to Mach 3 for an initial study
        #--------------------------------------------------------------------------------
        
        # Lift coefficients
        CLwb         = 0.00
        CLalpha_pdeg = 0.03
        CLdele_pdeg  = 0.005
        
        # Drag coefficients
        CDwb         = 0.062
        CDdele_pdeg  = 0
        CDdelsb_pdeg = 0
        
        # Side force coefficients
        CYbeta_prad    = -0.02           # Weathercock effect
        CYp_prps       =  0.017
        CYr_prps       = -0.92
        CYbetadot_pdps =  0
        CYdela_pdeg    = -0.0003
        CYdelr_pdeg    =  0.0057
        
        # Roll moment coefficients
        Clbeta_prad    =  0.00          # Dihedral effect, no natural tendency to roll back to level flight
        Clp_prps       = -0.185         # Roll damping, decent amount
        Clr_prps       = -0.080         # Reverse yaw-roll coupling 
        Clbetadot_pdps =  0             # Assumed zero
        Cldela_pdeg    =  0.0007        # Dominant in bank angle maneuverability
        Cldelr_pdeg    =  0.00025       # Cross-control derivative, rudder creates roll
        
        # Pitch moment coefficients
        Cmwb            =  0.00          
        Cmalpha_pdeg    = -0.01         # Stability derivative, stable in pitch
        Cmq_prps        = -4.4          # Pitch damping, strong effect
        Cmalphadot_prps = -0.1
        Cmdele_pdeg     = -0.008        # Control derivative, good pitch effectiveness
        Cmdelsb_pdeg    =  0
        
        # Yaw moment coefficients
        Cnbeta_pdeg     =  0.0046       # Stability derivative, stable (weathercock stability)
        Cnp_prps        = -0.012        # Yaw due to roll rate, negative value contributes to adverse yaw
        Cnr_prps        = -1.19         # Yaw damping, strong value
        Cnbetadot_pdps  =  0            # Assumed zero
        Cndela_pdeg     =  0.00045      # Cross-control derivative, positive value contributes to adverse yaw
        Cndelr_pdeg     = -0.0039       # Control derivative, good yaw control effectiveness
        
        # Moments and products of inertia fully burned out (Yancey64)
        Jxx_b_slugft2 = 3600
        Jxx_b_kgm2 = Jxx_b_slugft2*1.355
        Jxz_b_slugft2 = -700
        Jxz_b_kgm2 = Jxz_b_slugft2*1.355
        Jyy_b_slugft2 = 86000
        Jyy_b_kgm2 = Jyy_b_slugft2*1.355
        Jzz_b_slugft2 = 88500
        Jzz_b_kgm2 = Jzz_b_slugft2*1.355

        vmod = {"V_name"          : vehicle_name, 
                "m_kg"            : m_kg, 
                "Jxz_b_kgm2"      : Jxz_b_kgm2, 
                "Jxx_b_kgm2"      : Jxx_b_kgm2, 
                "Jyy_b_kgm2"      : Jyy_b_kgm2, 
                "Jzz_b_kgm2"      : Jzz_b_kgm2, 
                "CDwb"            : CDwb, 
                "CDdele_pdeg"     : CDdele_pdeg, 
                "CDdelsb_pdeg"    : CDdelsb_pdeg, 
                "CLwb"            : CLwb, 
                "CLalpha_pdeg"    : CLalpha_pdeg,
                "CLdele_pdeg"     : CLdele_pdeg, 
                "CYbeta_prad"     : CYbeta_prad, 
                "CYp_prps"        : CYp_prps, 
                "CYr_prps"        : CYr_prps, 
                "CYbetadot_pdps"  : CYbetadot_pdps, 
                "CYdela_pdeg"     : CYdela_pdeg, 
                "CYdelr_pdeg"     : CYdelr_pdeg, 
                "Clbeta_prad"     : Clbeta_prad, 
                "Clp_prps"        : Clp_prps, 
                "Clr_prps"        : Clr_prps, 
                "Clbetadot_pdps"  : Clbetadot_pdps, 
                "Cldela_pdeg"     : Cldela_pdeg, 
                "Cldelr_pdeg"     : Cldelr_pdeg, 
                "Cmwb"            : Cmwb,
                "Cmalpha_pdeg"    : Cmalpha_pdeg, 
                "Cmq_prps"        : Cmq_prps, 
                "Cmalphadot_prps" : Cmalphadot_prps, 
                "Cmdele_pdeg"     : Cmdele_pdeg, 
                "Cmdelsb_pdeg"    : Cmdelsb_pdeg, 
                "Cnbeta_pdeg"     : Cnbeta_pdeg, 
                "Cnp_prps"        : Cnp_prps,
                "Cnr_prps"        : Cnr_prps, 
                "Cnbetadot_pdps"  : Cnbetadot_pdps, 
                "Cndela_pdeg"     : Cndela_pdeg,   
                "Cndelr_pdeg"     : Cndelr_pdeg, 
                "A_ref_m2"        : A_ref_m2, 
                "b_m"             : b_m, 
                "c_m"             : c_m, 
                "short_name"      : short_name} 

        return vmod

def X15_Rev2():
        
        # Rocket motor burned out
        # Speed brake closed
        # Data taken from WT line on plot
        # Lower rudder on
        
        # Load aerodynamic database
        X15_aerodynamic_database = np.load('vehicle_models/X15/aerodynamic_model/X15_aerodynamic_database.npz')

        # Conversions
        ft2m = 0.3048

        # Name of vehicle (used in plot titles)
        vehicle_name = "X15 Space Plane"

        # Short name (used in file names)
        short_name = "X15"
        
        # Mass
        m_kg = 11765

        # Reference wing span
        b_m = 22.36*ft2m

        # Reference mean aerodynamic chord
        c_m = 10.27*ft2m
        
        # Reference area
        A_ref_m2 = 18.6
        
        # Aerodynamic coefficients 
        #--------------------------------------------------------------------------------
        
        # Independent variable breakpoints
        alpha_bps_deg = X15_aerodynamic_database["alpha_bps_deg"]
        Mach_bps      = X15_aerodynamic_database["Mach_bps"]
        
        # Lift coefficients
        CL_table_alpha_deg_Mach = X15_aerodynamic_database["CL_table_alpha_deg_Mach"]
        CLalpha_pdeg = 0.03
        CLdele_pdeg  = 0.005
        
        # Drag coefficients
        CDwb         = 0.062
        CDdele_pdeg  = 0
        CDdelsb_pdeg = 0
        
        # Side force coefficients
        CYbeta_prad    = -0.02           # Weathercock effect
        CYp_prps       =  0.017
        CYr_prps       = -0.92
        CYbetadot_pdps =  0
        CYdela_pdeg    = -0.0003
        CYdelr_pdeg    =  0.0057
        
        # Roll moment coefficients
        Clbeta_pdeg_table_Mach_alpha_deg = X15_aerodynamic_database["Clbeta_pdeg_table_Mach_alpha_deg"] # Dihedral effect, no natural tendency to roll back to level flight
        Clp_prps       = -0.185         # Roll damping, decent amount
        Clr_prps       = -0.080         # Reverse yaw-roll coupling 
        Clbetadot_pdps =  0             # Assumed zero
        Cldela_pdeg    =  0.0007        # Dominant in bank angle maneuverability
        Cldelr_pdeg    =  0.00025       # Cross-control derivative, rudder creates roll
        
        # Pitch moment coefficients
        Cmwb            =  0.00          
        Cmalpha_pdeg    = -0.01         # Stability derivative, stable in pitch
        Cmq_prps        = -4.4          # Pitch damping, strong effect
        Cmalphadot_prps = -0.1
        Cmdele_pdeg     = -0.008        # Control derivative, good pitch effectiveness
        Cmdelsb_pdeg    =  0
        
        # Yaw moment coefficients
        Cnbeta_pdeg     =  0.0046       # Stability derivative, stable (weathercock stability)
        Cnp_prps        = -0.012        # Yaw due to roll rate, negative value contributes to adverse yaw
        Cnr_prps        = -1.19         # Yaw damping, strong value
        Cnbetadot_pdps  =  0            # Assumed zero
        Cndela_pdeg     =  0.00045      # Cross-control derivative, positive value contributes to adverse yaw
        Cndelr_pdeg     = -0.0039       # Control derivative, good yaw control effectiveness
        
        # Moments and products of inertia fully burned out (Yancey64)
        #----------------------------------------------------------------------------------------------
        
        Jxx_b_slugft2 = 3600
        Jxx_b_kgm2 = Jxx_b_slugft2*1.355
        Jxz_b_slugft2 = -700
        Jxz_b_kgm2 = Jxz_b_slugft2*1.355
        Jyy_b_slugft2 = 86000
        Jyy_b_kgm2 = Jyy_b_slugft2*1.355
        Jzz_b_slugft2 = 88500
        Jzz_b_kgm2 = Jzz_b_slugft2*1.355

        vmod = {"V_name"                        : vehicle_name, 
                "m_kg"                          : m_kg, 
                "Jxz_b_kgm2"                    : Jxz_b_kgm2, 
                "Jxx_b_kgm2"                    : Jxx_b_kgm2, 
                "Jyy_b_kgm2"                    : Jyy_b_kgm2, 
                "Jzz_b_kgm2"                    : Jzz_b_kgm2, 
                "alpha_bps_deg"                 : alpha_bps_deg,
                "Mach_bps"                      : Mach_bps,
                "CDwb"                          : CDwb, 
                "CDdele_pdeg"                   : CDdele_pdeg, 
                "CDdelsb_pdeg"                  : CDdelsb_pdeg, 
                "CL_table_alpha_deg_Mach"       : CL_table_alpha_deg_Mach, 
                "CLalpha_pdeg"                  : CLalpha_pdeg,
                "CLdele_pdeg"                   : CLdele_pdeg, 
                "CYbeta_prad"                   : CYbeta_prad, 
                "CYp_prps"                      : CYp_prps, 
                "CYr_prps"                      : CYr_prps, 
                "CYbetadot_pdps"                : CYbetadot_pdps, 
                "CYdela_pdeg"                   : CYdela_pdeg, 
                "CYdelr_pdeg"                   : CYdelr_pdeg, 
                "Clbeta_pdeg_table_Mach_alpha_deg" : Clbeta_pdeg_table_Mach_alpha_deg, 
                "Clp_prps"                      : Clp_prps, 
                "Clr_prps"                      : Clr_prps, 
                "Clbetadot_pdps"                : Clbetadot_pdps, 
                "Cldela_pdeg"                   : Cldela_pdeg, 
                "Cldelr_pdeg"                   : Cldelr_pdeg, 
                "Cmwb"                          : Cmwb,
                "Cmalpha_pdeg"                  : Cmalpha_pdeg, 
                "Cmq_prps"                      : Cmq_prps, 
                "Cmalphadot_prps"               : Cmalphadot_prps, 
                "Cmdele_pdeg"                   : Cmdele_pdeg, 
                "Cmdelsb_pdeg"                  : Cmdelsb_pdeg, 
                "Cnbeta_pdeg"                   : Cnbeta_pdeg, 
                "Cnp_prps"                      : Cnp_prps,
                "Cnr_prps"                      : Cnr_prps, 
                "Cnbetadot_pdps"                : Cnbetadot_pdps, 
                "Cndela_pdeg"                   : Cndela_pdeg,   
                "Cndelr_pdeg"                   : Cndelr_pdeg, 
                "A_ref_m2"                      : A_ref_m2, 
                "b_m"                           : b_m, 
                "c_m"                           : c_m, 
                "short_name"                    : short_name} 

        return vmod

def X15_Rev3():
        
        # Rocket motor burned out
        # Speed brake closed
        # Data taken from WT line on plot
        # Lower rudder on
        
        # Load aerodynamic database
        X15_aerodynamic_database = np.load('vehicle_models/X15/aerodynamic_model/X15_aerodynamic_database.npz')

        # Conversions
        ft2m = 0.3048

        # Name of vehicle (used in plot titles)
        vehicle_name = "X15 Space Plane"

        # Short name (used in file names)
        short_name = "X15"
        
        # Mass
        m_kg = 11765

        # Reference wing span
        b_m = 22.36*ft2m

        # Reference mean aerodynamic chord
        c_m = 10.27*ft2m
        
        # Reference area
        A_ref_m2 = 18.6
        
        # Aerodynamic coefficients 
        #--------------------------------------------------------------------------------
        
        # Independent variable breakpoints
        alpha_bps_deg                       = X15_aerodynamic_database["alpha_bps_deg"]
        alpha_p_dele_bps_deg                = X15_aerodynamic_database["alpha_p_dele_bps_deg"]
        Mach_bps                            = X15_aerodynamic_database["Mach_bps"]
        
        # Lift coefficients
        CL_table_alpha_deg_Mach             = X15_aerodynamic_database["CL_table_alpha_deg_Mach"]
        CLalpha_pdeg                        = 0                                                             # Not used since CL_table above captures lift due to change in AoA
        CLdele_table_pdeg_v_Mach_AoA_p_dele = X15_aerodynamic_database["CLdele_table_pdeg_v_Mach_AoA_p_dele"]
        
        # Drag coefficients
        CDwb_table_v_alpha_Mach             = X15_aerodynamic_database["CD_table_v_alpha_Mach"]
        CDdele_pdeg                         = 0
        CDdelsb_pdeg                        = 0
        
        # Side force coefficients
        CYbeta_table_prad_v_alpha_deg_Mach = X15_aerodynamic_database["CYbeta_table_prad_v_alpha_deg_Mach"]
        CYp_table_prps_v_alpha_deg_Mach    = X15_aerodynamic_database["CYp_table_prps_v_alpha_deg_Mach"]
        CYr_table_prps_v_alpha_deg_Mach    = X15_aerodynamic_database["CYr_table_prps_v_alpha_deg_Mach"]
        CYbetadot_pdps                     =  0
        CYdela_table_pdeg_v_alpha_deg_Mach = X15_aerodynamic_database["CYdela_table_pdeg_v_alpha_deg_Mach"]
        CYdelr_table_pdeg_v_alpha_deg_Mach = X15_aerodynamic_database["CYdelr_table_pdeg_v_alpha_deg_Mach"]
        
        # Roll moment coefficients
        Clbeta_pdeg_table_Mach_alpha_deg   = X15_aerodynamic_database["Clbeta_pdeg_table_Mach_alpha_deg"] 
        Clp_prps_table_Mach_alpha_deg      = X15_aerodynamic_database["Clp_prps_table_Mach_alpha_deg"]
        Clr_prps_table_Mach_alpha_deg      = X15_aerodynamic_database["Clr_prps_table_Mach_alpha_deg"]
        Clbetadot_pdps                     =  0         
        Cldela_table_pdeg_v_alpha_deg_Mach = X15_aerodynamic_database["Cldela_table_pdeg_v_alpha_deg_Mach"]
        Cldelr_table_pdeg_v_alpha_deg_Mach = X15_aerodynamic_database["Cldelr_table_pdeg_v_alpha_deg_Mach"]
        
        # Pitch moment coefficients
        Cm_table_alpha_deg_Mach             = X15_aerodynamic_database["Cm_table_alpha_deg_Mach"]
        Cmalpha_pdeg_table_alpha_deg_Mach   = X15_aerodynamic_database["Cmalpha_pdeg_table_alpha_deg_Mach"] 
        Cmq_pdps_table_alpha_deg_Mach       = X15_aerodynamic_database["Cmq_pdps_table_alpha_deg_Mach"] 
        Cmalphadot_prps                     = 0
        Cmdele_pdeg_table_AoApdele_deg_Mach = X15_aerodynamic_database["Cmdele_pdeg_table_AoApdele_deg_Mach"]
        Cmdelsb_pdeg                        = 0
        
        # Yaw moment coefficients
        Cnbeta_table_prad_v_alpha_deg_Mach = X15_aerodynamic_database["Cnbeta_table_prad_v_alpha_deg_Mach"]
        Cnp_prps_table_Mach_alpha_deg      = X15_aerodynamic_database["Cnp_prps_table_Mach_alpha_deg"]
        Cnr_prps_table_Mach_alpha_deg      = X15_aerodynamic_database["Cnr_prps_table_Mach_alpha_deg"]
        Cnbetadot_pdps                     = 0         
        Cndela_table_pdeg_v_alpha_deg_Mach = X15_aerodynamic_database["Cndela_table_pdeg_v_alpha_deg_Mach"]
        Cndelr_table_pdeg_v_alpha_deg_Mach = X15_aerodynamic_database["Cndelr_table_pdeg_v_alpha_deg_Mach"]
        
        # Moments and products of inertia fully burned out (Yancey64)
        #----------------------------------------------------------------------------------------------
        Jxx_b_slugft2 = 3600
        Jxx_b_kgm2 = Jxx_b_slugft2*1.355
        Jxz_b_slugft2 = -700
        Jxz_b_kgm2 = Jxz_b_slugft2*1.355
        Jyy_b_slugft2 = 86000
        Jyy_b_kgm2 = Jyy_b_slugft2*1.355
        Jzz_b_slugft2 = 88500
        Jzz_b_kgm2 = Jzz_b_slugft2*1.355

        vmod = {"V_name"                                : vehicle_name, 
                "m_kg"                                  : m_kg, 
                "Jxz_b_kgm2"                            : Jxz_b_kgm2, 
                "Jxx_b_kgm2"                            : Jxx_b_kgm2, 
                "Jyy_b_kgm2"                            : Jyy_b_kgm2, 
                "Jzz_b_kgm2"                            : Jzz_b_kgm2, 
                "alpha_bps_deg"                         : alpha_bps_deg,
                "alpha_p_dele_bps_deg"                  : alpha_p_dele_bps_deg,
                "Mach_bps"                              : Mach_bps,
                "CDwb_table_v_alpha_Mach"               : CDwb_table_v_alpha_Mach,
                "CDdele_pdeg"                           : CDdele_pdeg, 
                "CDdelsb_pdeg"                          : CDdelsb_pdeg, 
                "CL_table_alpha_deg_Mach"               : CL_table_alpha_deg_Mach, 
                "CLalpha_pdeg"                          : CLalpha_pdeg,
                "CLdele_table_pdeg_v_Mach_AoA_p_dele"   : CLdele_table_pdeg_v_Mach_AoA_p_dele,  
                "CYbeta_table_prad_v_alpha_deg_Mach"    : CYbeta_table_prad_v_alpha_deg_Mach,
                "CYp_table_prps_v_alpha_deg_Mach"       : CYp_table_prps_v_alpha_deg_Mach,
                "CYr_table_prps_v_alpha_deg_Mach"       : CYr_table_prps_v_alpha_deg_Mach,
                "CYbetadot_pdps"                        : CYbetadot_pdps, 
                "CYdela_table_pdeg_v_alpha_deg_Mach"    : CYdela_table_pdeg_v_alpha_deg_Mach,
                "CYdelr_table_pdeg_v_alpha_deg_Mach"    : CYdelr_table_pdeg_v_alpha_deg_Mach,
                "Clbeta_pdeg_table_Mach_alpha_deg"      : Clbeta_pdeg_table_Mach_alpha_deg, 
                "Clp_prps_table_Mach_alpha_deg"         : Clp_prps_table_Mach_alpha_deg,
                "Clr_prps_table_Mach_alpha_deg"         : Clr_prps_table_Mach_alpha_deg,
                "Clbetadot_pdps"                        : Clbetadot_pdps, 
                "Cldela_table_pdeg_v_alpha_deg_Mach"    : Cldela_table_pdeg_v_alpha_deg_Mach,
                "Cldelr_table_pdeg_v_alpha_deg_Mach"    : Cldelr_table_pdeg_v_alpha_deg_Mach,
                "Cm_table_alpha_deg_Mach"               : Cm_table_alpha_deg_Mach,
                "Cmalpha_pdeg_table_alpha_deg_Mach"     : Cmalpha_pdeg_table_alpha_deg_Mach,
                "Cmq_pdps_table_alpha_deg_Mach"         : Cmq_pdps_table_alpha_deg_Mach, 
                "Cmalphadot_prps"                       : Cmalphadot_prps, 
                "Cmdele_pdeg_table_AoApdele_deg_Mach"   : Cmdele_pdeg_table_AoApdele_deg_Mach, 
                "Cmdelsb_pdeg"                          : Cmdelsb_pdeg, 
                "Cnbeta_table_prad_v_alpha_deg_Mach"    : Cnbeta_table_prad_v_alpha_deg_Mach, 
                "Cnp_prps_table_Mach_alpha_deg"         : Cnp_prps_table_Mach_alpha_deg,
                "Cnr_prps_table_Mach_alpha_deg"         : Cnr_prps_table_Mach_alpha_deg,
                "Cnbetadot_pdps"                        : Cnbetadot_pdps, 
                "Cndela_table_pdeg_v_alpha_deg_Mach"    : Cndela_table_pdeg_v_alpha_deg_Mach,
                "Cndelr_table_pdeg_v_alpha_deg_Mach"    : Cndelr_table_pdeg_v_alpha_deg_Mach,
                "A_ref_m2"                              : A_ref_m2, 
                "b_m"                                   : b_m, 
                "c_m"                                   : c_m, 
                "short_name"                            : short_name} 

        return vmod

def X15_Dutch_Roll_Manipulated():
        
         # Rocket motor burned out
        # Speed brake closed
        # Data taken from WT line on plot
        # Lower rudder on
        
        # Load aerodynamic database
        X15_aerodynamic_database = np.load('vehicle_models/X15/aerodynamic_model/X15_aerodynamic_database.npz')

        # Conversions
        ft2m = 0.3048

        # Name of vehicle (used in plot titles)
        vehicle_name = "X15 Space Plane"

        # Short name (used in file names)
        short_name = "X15"
        
        # Mass
        m_kg = 11765

        # Reference wing span
        b_m = 22.36*ft2m

        # Reference mean aerodynamic chord
        c_m = 10.27*ft2m
        
        # Reference area
        A_ref_m2 = 18.6
        
        # Aerodynamic coefficients 
        #--------------------------------------------------------------------------------
        
        # Independent variable breakpoints
        alpha_bps_deg        = X15_aerodynamic_database["alpha_bps_deg"]
        alpha_p_dele_bps_deg = X15_aerodynamic_database["alpha_p_dele_bps_deg"]
        Mach_bps             = X15_aerodynamic_database["Mach_bps"]
        
        # Lift coefficients
        CL_table_alpha_deg_Mach = X15_aerodynamic_database["CL_table_alpha_deg_Mach"]
        CLalpha_pdeg = 0.03
        CLdele_table_pdeg_v_Mach_AoA_p_dele = X15_aerodynamic_database["CLdele_table_pdeg_v_Mach_AoA_p_dele"]
        
        # Drag coefficients
        CDwb         = 0.062
        CDdele_pdeg  = 0
        CDdelsb_pdeg = 0
        
        # Side force coefficients
        CYbeta_prad    = -0.60          # Weathercock effect
        CYp_prps       =  0.017
        CYr_prps       = -0.92
        CYbetadot_pdps =  0
        CYdela_pdeg    = -0.0003
        CYdelr_pdeg    =  0.0057
        
        # Roll moment coefficients
        Clbeta_pdeg_table_Mach_alpha_deg = X15_aerodynamic_database["Clbeta_pdeg_table_Mach_alpha_deg"] # Dihedral effect, no natural tendency to roll back to level flight
        Clp_prps       = -0.185*2       # Roll damping, decent amount
        Clr_prps       =  0.080         # Yaw-roll coupling 
        Clbetadot_pdps =  0             # Assumed zero
        Cldela_pdeg    =  0.0007        # Dominant in bank angle maneuverability
        Cldelr_pdeg    =  0.00025       # Cross-control derivative, rudder creates roll
        
        # Pitch moment coefficients
        Cmwb            =  0.00          
        Cmalpha_pdeg    = -0.01         # Stability derivative, stable in pitch
        Cmq_prps        = -4.4          # Pitch damping, strong effect
        Cmalphadot_prps = -0.1
        Cmdele_pdeg     = -0.008        # Control derivative, good pitch effectiveness
        Cmdelsb_pdeg    =  0
        
        # Yaw moment coefficients
        Cnbeta_pdeg     =  0.19         # Stability derivative, stable (weathercock stability)
        Cnp_prps        = -0.012        # Yaw due to roll rate, negative value contributes to adverse yaw
        Cnr_prps        = -1.190        # Yaw damping, strong value
        Cnbetadot_pdps  =  0            # Assumed zero
        Cndela_pdeg     =  0.00045      # Cross-control derivative, positive value contributes to adverse yaw
        Cndelr_pdeg     = -0.0039       # Control derivative, good yaw control effectiveness
        
        # Moments and products of inertia fully burned out (Yancey64)
        #----------------------------------------------------------------------------------------------
        
        Jxx_b_slugft2 = 3600
        Jxx_b_kgm2 = Jxx_b_slugft2*1.355
        Jxz_b_slugft2 = -700
        Jxz_b_kgm2 = Jxz_b_slugft2*1.355
        Jyy_b_slugft2 = 86000
        Jyy_b_kgm2 = Jyy_b_slugft2*1.355
        Jzz_b_slugft2 = 88500
        Jzz_b_kgm2 = Jzz_b_slugft2*1.355

        vmod = {"V_name"                        : vehicle_name, 
                "m_kg"                          : m_kg, 
                "Jxz_b_kgm2"                    : Jxz_b_kgm2, 
                "Jxx_b_kgm2"                    : Jxx_b_kgm2, 
                "Jyy_b_kgm2"                    : Jyy_b_kgm2, 
                "Jzz_b_kgm2"                    : Jzz_b_kgm2, 
                "alpha_bps_deg"                 : alpha_bps_deg,
                "alpha_p_dele_bps_deg"          : alpha_p_dele_bps_deg,
                "Mach_bps"                      : Mach_bps,
                "CDwb"                          : CDwb, 
                "CDdele_pdeg"                   : CDdele_pdeg, 
                "CDdelsb_pdeg"                  : CDdelsb_pdeg, 
                "CL_table_alpha_deg_Mach"       : CL_table_alpha_deg_Mach, 
                "CLalpha_pdeg"                  : CLalpha_pdeg,
                "CLdele_table_pdeg_v_Mach_AoA_p_dele" : CLdele_table_pdeg_v_Mach_AoA_p_dele, 
                "CYbeta_prad"                   : CYbeta_prad, 
                "CYp_prps"                      : CYp_prps, 
                "CYr_prps"                      : CYr_prps, 
                "CYbetadot_pdps"                : CYbetadot_pdps, 
                "CYdela_pdeg"                   : CYdela_pdeg, 
                "CYdelr_pdeg"                   : CYdelr_pdeg, 
                "Clbeta_pdeg_table_Mach_alpha_deg" : Clbeta_pdeg_table_Mach_alpha_deg, 
                "Clp_prps"                      : Clp_prps, 
                "Clr_prps"                      : Clr_prps, 
                "Clbetadot_pdps"                : Clbetadot_pdps, 
                "Cldela_pdeg"                   : Cldela_pdeg, 
                "Cldelr_pdeg"                   : Cldelr_pdeg, 
                "Cmwb"                          : Cmwb,
                "Cmalpha_pdeg"                  : Cmalpha_pdeg, 
                "Cmq_prps"                      : Cmq_prps, 
                "Cmalphadot_prps"               : Cmalphadot_prps, 
                "Cmdele_pdeg"                   : Cmdele_pdeg, 
                "Cmdelsb_pdeg"                  : Cmdelsb_pdeg, 
                "Cnbeta_pdeg"                   : Cnbeta_pdeg, 
                "Cnp_prps"                      : Cnp_prps,
                "Cnr_prps"                      : Cnr_prps, 
                "Cnbetadot_pdps"                : Cnbetadot_pdps, 
                "Cndela_pdeg"                   : Cndela_pdeg,   
                "Cndelr_pdeg"                   : Cndelr_pdeg, 
                "A_ref_m2"                      : A_ref_m2, 
                "b_m"                           : b_m, 
                "c_m"                           : c_m, 
                "short_name"                    : short_name} 

        return vmod