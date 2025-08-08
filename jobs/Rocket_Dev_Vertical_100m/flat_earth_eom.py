import math
import numpy as np 
from tools.Interpolators import fastInterp1
from vehicle_models.brick.aerodynamics.roll_coef_brick  import Cl_brick
from vehicle_models.brick.aerodynamics.pitch_coef_brick import Cm_brick
from vehicle_models.brick.aerodynamics.yaw_coef_brick   import Cn_brick

# Example usage:
def flat_earth_eom(t, x, vmod, amod, cmod):
  """FUNCTION flat_earth_eom.py contains the essential elements of a six-degree of freedom 
  simulation. The purpose of this function is to allow the numerical approximation of 
  solutions of the governing equations for an aircraft.
  
  The naming convention is <variable name>_<coordinate system if applicable>_<units>. For 
  example, the pitch rate, q, resolved in the body fixed frame, bf, with units of radians
  per second is named, q_b_rps.

  Arguments: 
    t - time [s], scalar
    x - state vector at time t [various units], numpy array
      x[0]  = u_b_mps, axial velocity of CM wrt inertial CS resolved in aircraft body fixed CS
      x[1]  = v_b_mps, lateral velocity of CM wrt inertial CS resolved in aircraft body fixed CS
      x[2]  = w_b_mps, vertical velocity of CM wrt inertial CS resolved in aircraft body fixed CS
      x[3]  = p_b_rps, roll angular velocity of body fixed CS with respect to inertial CS
      x[4]  = q_b_rps, pitch angular velocity of body fixed CS with respect to inertial CS
      x[5]  = r_b_rps, yaw angular velocity of body fixed CS with respect to inertial CS
      x[6]  = phi_rad, roll angle
      x[7]  = theta_rad, pitch angle
      x[8]  = psi_rad, yaw angle
      x[9]  = p1_n_m, x-axis position of aircraft resolved in NED CS
      x[10] = p2_n_m, y-axis position of aircraft resolved in NED CS
      x[11] = p3_n_m, z-axis position of aircraft resolved in NED CS
    vmod - vehicle model data stored as a dictionary containing various parameters
    amod - atmosphere and gravity model data stored as a dictionary 

  Returns:
    dx - the time derivative of each state in x (RHS of governing equations)

  History:
    Written by Ben Dickinson
      - March 2024 Six degree of freedom equations written 
      - April 2024 Euler kinematic equations added 
      - May 2024   Added USSA 1976 tabular data for atmosphere and gravity
  """
  
  # Preallocate left hand side of state equations
  dx = np.empty( (12,), dtype=float ) 
  auxillary_data = np.empty( (1,), dtype=float )

  # Assign current state values to variable names
  u_b_mps   = x[0]
  v_b_mps   = x[1]
  w_b_mps   = x[2]
  p_b_rps   = x[3]
  q_b_rps   = x[4]
  r_b_rps   = x[5]
  phi_rad   = x[6]
  theta_rad = x[7]
  psi_rad   = x[8]
  p1_n_m    = x[9]
  p2_n_m    = x[10]
  p3_n_m    = x[11]
  
  # Compute trig operations on Euler angles 
  c_phi   = math.cos(phi_rad)
  c_theta = math.cos(theta_rad)
  c_psi   = math.cos(psi_rad)
  s_phi   = math.sin(phi_rad)
  s_theta = math.sin(theta_rad)
  s_psi   = math.sin(psi_rad)
  t_theta = math.tan(theta_rad)
  
  # Compute DCM
  C_b2n_11 =  c_theta*c_psi
  C_b2n_12 = -c_phi*s_psi + s_phi*s_theta*c_psi
  C_b2n_13 =  s_phi*s_psi + c_phi*s_theta*c_psi
  C_b2n_21 =  c_theta*s_psi
  C_b2n_22 =  c_phi*c_psi + s_phi*s_theta*s_psi
  C_b2n_23 = -s_phi*c_psi + c_phi*s_theta*s_psi
  C_b2n_31 = -s_theta
  C_b2n_32 =  s_phi*c_theta
  C_b2n_33 =  c_phi*c_theta
  C_n2b_13 = C_b2n_31
  C_n2b_23 = C_b2n_32
  C_n2b_33 = C_b2n_33

  # Get mass and moments of inertia
  m_kg       = vmod["m_kg"]
  Jxz_b_kgm2 = vmod["Jxz_b_kgm2"]
  Jxx_b_kgm2 = vmod["Jxx_b_kgm2"]
  Jyy_b_kgm2 = vmod["Jyy_b_kgm2"]
  Jzz_b_kgm2 = vmod["Jzz_b_kgm2"]
  
  # Get reference dimensions
  A_ref_m2 = vmod["A_ref_m2"]
  b_m     = vmod["b_m"]
  c_m     = vmod["c_m"]
  
  # Get aerodynamic coefficients / tables
  Clp = vmod["Clp"]
  Clr = vmod["Clr"]
  Cmq = vmod["Cmq"]
  Cnp = vmod["Cnp"]
  Cnr = vmod["Cnr"]
                
  # Get current altitude
  h_m = -p3_n_m
  
  # US Standard Atmosphere 1976 data
  rho_interp_kgpm3 = fastInterp1(amod["alt_m"], amod["rho_kgpm3"], h_m)
  c_interp_mp2     = fastInterp1(amod["alt_m"], amod["c_mps"], h_m)
  
  # Air data calculation (Mach, AoA, AoS)
  true_airspeed_mps = math.sqrt(u_b_mps**2 + v_b_mps**2 + w_b_mps**2)
  qbar_kgpms2       = 0.5*rho_interp_kgpm3*true_airspeed_mps**2
  Mach              = true_airspeed_mps/c_interp_mp2
  alpha_rad         = math.atan2(w_b_mps, u_b_mps)
  beta_rad          = math.asin(v_b_mps/true_airspeed_mps)
  
  # Trig calcs for body to wind/stability axes DCM
  s_alpha   = math.sin(alpha_rad)
  c_alpha   = math.cos(alpha_rad)
  s_beta    = math.sin(beta_rad)
  c_beta    = math.cos(beta_rad)
  C_w2b_11 =  c_alpha*c_beta
  C_w2b_12 = -c_alpha*s_beta
  C_w2b_13 = -s_alpha
  C_w2b_21 =  s_beta
  C_w2b_22 =  c_beta
  C_w2b_23 =  0
  C_w2b_31 =  s_alpha*c_beta
  C_w2b_32 = -s_alpha*s_beta
  C_w2b_33 =  c_alpha

  # Gravity acts normal to earth tangent CS 
  gz_interp_n_mps2 = fastInterp1(amod["alt_m"], amod['g_mps2'], h_m)

  # Resolve gravity in body coordinate system
  gx_b_mps2 =  C_n2b_13*gz_interp_n_mps2
  gy_b_mps2 =  C_n2b_23*gz_interp_n_mps2
  gz_b_mps2 =  C_n2b_33*gz_interp_n_mps2

  # Propulsive forces 
  t_interp_N = fastInterp1(vmod["motor_thrust_time_axis"], vmod["motor_thrust_N"], t)

  FPx_b_N = t_interp_N
  FPy_b_N = 0
  FPz_b_N = 0

  # Aerodynamic forces
  drag_kgmps2 = 0
  side_kgmps2 = 0
  lift_kgmps2 = 0

  # External forces 
  Fx_b_kgmps2 = -(C_w2b_11*drag_kgmps2 + C_w2b_12*side_kgmps2 + C_w2b_13*lift_kgmps2) + FPx_b_N
  Fy_b_kgmps2 = -(C_w2b_21*drag_kgmps2 + C_w2b_22*side_kgmps2 + C_w2b_23*lift_kgmps2) + FPy_b_N
  Fz_b_kgmps2 = -(C_w2b_31*drag_kgmps2 + C_w2b_32*side_kgmps2 + C_w2b_33*lift_kgmps2) + FPz_b_N
  
  # External moments
  l_b_kgm2ps2 = Cl_brick(Clp, Clr, p_b_rps, r_b_rps, b_m, true_airspeed_mps)*qbar_kgpms2*A_ref_m2*b_m
  m_b_kgm2ps2 = Cm_brick(Cmq, q_b_rps, c_m, true_airspeed_mps)*qbar_kgpms2*A_ref_m2*c_m
  n_b_kgm2ps2 = Cn_brick(Cnp, Cnr, p_b_rps, r_b_rps, b_m, true_airspeed_mps)*qbar_kgpms2*A_ref_m2*b_m
  
  # Denominator in roll and yaw rate equations
  Gamma_inv = 1/(Jxx_b_kgm2*Jzz_b_kgm2 - Jxz_b_kgm2**2)
  
  # x-axis (roll-axis) velocity equation
  #  State: u_b_mps
  dx[0] = 1/m_kg*Fx_b_kgmps2 + gx_b_mps2 - w_b_mps*q_b_rps + v_b_mps*r_b_rps

  # y-axis (pitch-axis) velocity equation
  #  State: v_b_mps
  dx[1] = 1/m_kg*Fy_b_kgmps2 + gy_b_mps2 - u_b_mps*r_b_rps + w_b_mps*p_b_rps
  
  # z-axis (yaw-axis) velocity equation
  #  State: w_b_mps
  dx[2] = 1/m_kg*Fz_b_kgmps2 + gz_b_mps2 - v_b_mps*p_b_rps + u_b_mps*q_b_rps
  
  # Roll equation
  #  State: p_b_rps
  dx[3] = (Jxz_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2  + Jzz_b_kgm2)    * p_b_rps * q_b_rps - \
          (Jzz_b_kgm2 * (Jzz_b_kgm2 - Jyy_b_kgm2) + Jxz_b_kgm2**2) * q_b_rps * r_b_rps +  \
           Jzz_b_kgm2 * l_b_kgm2ps2 + \
           Jxz_b_kgm2 * n_b_kgm2ps2)*Gamma_inv
            
  # Pitch equation
  #  State: q_b_rps
  dx[4] = ((Jzz_b_kgm2 - Jxx_b_kgm2) * p_b_rps * r_b_rps - \
           Jxz_b_kgm2 * (p_b_rps**2 - r_b_rps**2) + m_b_kgm2ps2)/Jyy_b_kgm2
  
  # Yaw equation
  #  State: r_b_rps
  dx[5] = ((Jxx_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2) + Jxz_b_kgm2**2) * p_b_rps * q_b_rps - \
            Jxz_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2 + Jzz_b_kgm2)     * q_b_rps * r_b_rps + \
            Jxz_b_kgm2 * l_b_kgm2ps2 + \
            Jxx_b_kgm2 * n_b_kgm2ps2)*Gamma_inv
  
  # Kinematic equations
  dx[6] = p_b_rps + s_phi*t_theta*q_b_rps + c_phi*t_theta*r_b_rps
  dx[7] =                   c_phi*q_b_rps -         s_phi*r_b_rps
  dx[8] =           s_phi/c_theta*q_b_rps + c_phi/c_theta*r_b_rps

  # Position (Navigation) equations 
  dx[9]   =  C_b2n_11*u_b_mps + C_b2n_12*v_b_mps + C_b2n_13*w_b_mps
  dx[10]  =  C_b2n_21*u_b_mps + C_b2n_22*v_b_mps + C_b2n_23*w_b_mps
  dx[11]  =  C_b2n_31*u_b_mps + C_b2n_32*v_b_mps + C_b2n_33*w_b_mps

  # No auxillary data in this one
  auxillary_data[0] = 0

  return dx, auxillary_data