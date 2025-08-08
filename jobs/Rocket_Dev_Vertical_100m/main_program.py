import math

import ussa1976
import numpy as np
import matplotlib.pyplot as plt

from . import sim_config
from . import flat_earth_eom
from tools.Interpolators import fastInterp1
from numerical_integrators import numerical_integration_methods


def run_job():
    
    # Conversions
    r2d = 180/math.pi
    d2r = 1/r2d

    # Set time conditions
    t0_s = 0.0
    tf_s = 5
    h_s  = 0.01

    # Set initial conditions
    u0_bf_mps  =   1e-10                    # Velocity x body axis
    v0_bf_mps  =   0                        # Velocity y body axis
    w0_bf_mps  =   0                        # Velocity z body axis
    p0_bf_rps  =   0*d2r                    # Rotation x body axis
    q0_bf_rps  =   0*d2r                    # Rotation y body axis
    r0_bf_rps  =   0*d2r                    # Rotation z body axis
    phi0_rad   =   0*d2r                    # Roll angle
    theta0_rad =   0*d2r                   # Pitch angle
    psi0_rad   =   0*d2r                    # Yaw angle
    p10_n_m    =   0
    p20_n_m    =   0
    p30_n_m    =   -100

    #==============================================================================
    # Part 1: Initialization of simulation
    #==============================================================================

    # A. Atmospheric data
    atmosphere = ussa1976.compute()

    # Get essential gravity and atmospheric data
    alt_m     = atmosphere["z"].values
    rho_kgpm3 = atmosphere["rho"].values
    c_mps     = atmosphere["cs"].values
    g_mps2    = ussa1976.core.compute_gravity(alt_m)

    amod = {
        "alt_m"    : alt_m,
        "rho_kgpm3" : rho_kgpm3,
        "c_mps"     : c_mps,
        "g_mps2"    : g_mps2
    } # Atmosphere and gravity data

    x0 = np.array([
        u0_bf_mps,  # x-axis body-fixed velocity (m/s)
        v0_bf_mps,  # y-axis body-fixed velocity (m/s)
        w0_bf_mps,  # z-axis body-fixed velocity (m/s)
        p0_bf_rps,  # roll rate (rad/s)
        q0_bf_rps,  # pitch rate (rad/s)
        r0_bf_rps,  # yaw rate (rad/s)
        phi0_rad,   # roll angle (rad)
        theta0_rad, # pitch angle (rad)
        psi0_rad,   # yaw angle (rad)
        p10_n_m,    # x-axis position (N*m)
        p20_n_m,    # y-axis position (N*m)
        p30_n_m,    # z-axis position (N*m)
    ])

    # Get number of elements in x0
    nx0 = x0.size

    #==============================================================================
    # Part 2: Numerically approximate solutions to the governing equations
    #==============================================================================

    # Preallocate the solution array
    t_s = np.arange( t0_s, tf_s + h_s, h_s ); nt_s = t_s.size
    x   = np.zeros((nx0, nt_s))

    # Preallocate auxillary data (the number of auxillary variables may change and it's
    # important to make sure that its dimension is consistent with the auxillary_data
    # that is returned by flat_earth_eom)
    No_of_Auxillary_Variables = 1
    Auxillary_Data_Accumulated = np.zeros((No_of_Auxillary_Variables, nt_s))

    # Assign the initial condition, x0, to solution array, x
    x[:, 0] = x0 

    # A. Perform forward Euler integration
    t_s, x, Auxillary_Data_Accumulated = numerical_integration_methods.RK4(flat_earth_eom.flat_earth_eom, \
        t_s, x, h_s, sim_config.vmod, amod, sim_config.cmod, Auxillary_Data_Accumulated)

    #==============================================================================
    # Part 3: Post-process simulation data
    #==============================================================================
        
    # Preallocate variables
    Altitude_m  = np.zeros((nt_s,1))
    Cs_mps      = np.zeros((nt_s,1))
    Rho_kgpm3   = np.zeros((nt_s,1))
    C_phi       = np.zeros((nt_s,1))
    C_theta     = np.zeros((nt_s,1))
    C_psi       = np.zeros((nt_s,1))
    S_phi       = np.zeros((nt_s,1))
    S_theta     = np.zeros((nt_s,1))
    S_psi       = np.zeros((nt_s,1))
    T_theta     = np.zeros((nt_s,1))
    C_b2n_11    = np.zeros((nt_s,1))
    C_b2n_12    = np.zeros((nt_s,1))
    C_b2n_13    = np.zeros((nt_s,1))
    C_b2n_21    = np.zeros((nt_s,1))
    C_b2n_22    = np.zeros((nt_s,1))
    C_b2n_23    = np.zeros((nt_s,1))
    C_b2n_31    = np.zeros((nt_s,1))
    C_b2n_32    = np.zeros((nt_s,1))
    C_b2n_33    = np.zeros((nt_s,1))
    u_n_mps     = np.zeros((nt_s,1)) 
    v_n_mps     = np.zeros((nt_s,1))
    w_n_mps     = np.zeros((nt_s,1)) 
    phi_2_rad   = np.zeros((nt_s,1)) 
    theta_2_rad = np.zeros((nt_s,1))
    psi_2_rad   = np.zeros((nt_s,1)) 

    for i, element in enumerate(t_s):
        Altitude_m[i,0] = -x[11,i]
        Cs_mps[i,0]      = fastInterp1(amod["alt_m"], amod["c_mps"],    Altitude_m[i,0])
        Rho_kgpm3[i,0]   = fastInterp1(amod["alt_m"], amod["rho_kgpm3"], Altitude_m[i,0])
        C_phi[i,0]       = math.cos(x[6,i])
        C_theta[i,0]     = math.cos(x[7,i])
        C_psi[i,0]       = math.cos(x[8,i])
        S_phi[i,0]       = math.sin(x[6,i])
        S_theta[i,0]     = math.sin(x[7,i])
        S_psi[i,0]       = math.sin(x[8,i])
        T_theta[i,0]     = math.tan(x[7,i])
        C_b2n_11[i,0]    =  C_theta[i,0]*C_psi[i,0]
        C_b2n_12[i,0]    = -C_phi[i,0]*S_psi[i,0] + S_phi[i,0]*S_theta[i,0]*C_psi[i,0]
        C_b2n_13[i,0]    =  S_phi[i,0]*S_psi[i,0] + C_phi[i,0]*S_theta[i,0]*C_psi[i,0]
        C_b2n_21[i,0]    =  C_theta[i,0]*S_psi[i,0]
        C_b2n_22[i,0]    =  C_phi[i,0]*C_psi[i,0] + S_phi[i,0]*S_theta[i,0]*S_psi[i,0]
        C_b2n_23[i,0]    = -S_phi[i,0]*C_psi[i,0] + C_phi[i,0]*S_theta[i,0]*S_psi[i,0]
        C_b2n_31[i,0]    = -S_theta[i,0]
        C_b2n_32[i,0]    =  S_phi[i,0]*C_theta[i,0]
        C_b2n_33[i,0]    =  C_phi[i,0]*C_theta[i,0]
        u_n_mps[i,0]     =  C_b2n_11[i,0]*x[0,i] + C_b2n_12[i,0]*x[1,i] + C_b2n_13[i,0]*x[2,i]
        v_n_mps[i,0]     =  C_b2n_21[i,0]*x[0,i] + C_b2n_22[i,0]*x[1,i] + C_b2n_23[i,0]*x[2,i]
        w_n_mps[i,0]     =  C_b2n_31[i,0]*x[0,i] + C_b2n_32[i,0]*x[1,i] + C_b2n_33[i,0]*x[2,i]
        phi_2_rad[i,0]   =  math.atan2(C_b2n_32[i,0],C_b2n_33[i,0]) 
        theta_2_rad[i,0] = -math.asin(C_b2n_31[i,0])
        psi_2_rad[i,0]   =  math.atan2(C_b2n_21[i,0],C_b2n_11[i,0])
        
        
    # Airspeed
    True_Airspeed_mps  = np.zeros((nt_s,1))
    for i, element in enumerate(t_s):
        True_Airspeed_mps[i,0] = math.sqrt(x[0,i]**2 + x[1,i]**2 + x[2,i]**2)
        
    # Angle of attack
    Alpha_rad = np.zeros((nt_s,1))
    for i, element in enumerate(t_s):  
        if x[0,i] == 0:
            w_over_u = 0
        else:
            w_over_u = x[2,i]/x[0,i]    
        Alpha_rad[i,0] = math.atan(w_over_u)
        
    # Angle of side slip
    Beta_rad = np.zeros((nt_s,1))
    for i, element in enumerate(t_s):  
        if True_Airspeed_mps[i,0] == 0:
            v_over_VT = 0
        else:
            v_over_VT = x[1,i]/True_Airspeed_mps[i,0]
        Beta_rad[i,0] = math.asin(v_over_VT)
        
    # Mach Number
    Mach = np.zeros((nt_s,1))
    for i, element in enumerate(t_s):
        Mach[i,0] = True_Airspeed_mps[i,0]/Cs_mps[i,0]

    # Ground hit distance
    Ground_hit_idx = np.argwhere(Altitude_m == Altitude_m[Altitude_m >= 0][-1])[0][0]
    Ground_hit_distance_north_m =  float(x[9,:][Ground_hit_idx])
    Ground_hit_s = float(t_s[Ground_hit_idx])

    #==============================================================================
    # Part 5: Plot data in figures
    #==============================================================================

    print(f"{Ground_hit_distance_north_m = }")
    print(f"{Ground_hit_s = }")

    # A. Create figure of translational and rotation states (6 DOF States in the body coordinate system)
    fig, axes = plt.subplots(2, 3, figsize=(10, 6))
    fig.set_facecolor('black')  
    title_suffix = f"{sim_config.vmod['V_name']}"
    title_prefix = "6 Degree of Freedom States\n"
    fig.suptitle(title_prefix + title_suffix, fontsize=14, fontweight='normal', color='yellow')

    # Axial velocity u^b_CM/n
    axes[0, 0].plot(t_s, x[0,:], color='yellow')
    axes[0, 0].set_xlabel('Time [s]', color='white')
    axes[0, 0].set_ylabel('u Velocity Body [m/s]', color='white')
    axes[0, 0].grid(True)
    axes[0, 0].set_facecolor('black')
    axes[0, 0].tick_params(colors = 'white')

    # y-axis velocity v^b_CM/n
    axes[0, 1].plot(t_s, x[1,:], color='yellow')
    axes[0, 1].set_xlabel('Time [s]', color='white')
    axes[0, 1].set_ylabel('v Velocity Body [m/s]', color='white')
    axes[0, 1].grid(True)
    axes[0, 1].set_facecolor('black')
    axes[0, 1].tick_params(colors = 'white')

    # z-axis velocity w^b_CM/n
    axes[0, 2].plot(t_s, x[2,:], color='yellow')
    axes[0, 2].set_xlabel('Time [s]', color='white')
    axes[0, 2].set_ylabel('w Velocity Body [m/s]', color='white')
    if np.linalg.norm(x[2,:]) < 1e-5:
        axes[0,2].set_ylim(-0.05,0.05)
    axes[0, 2].grid(True)
    axes[0, 2].set_facecolor('black')
    axes[0, 2].tick_params(colors = 'white')

    # Roll rate p^b_b/n
    axes[1, 0].plot(t_s, r2d*x[3,:], color='yellow')
    axes[1, 0].set_xlabel('Time [s]', color='white')
    axes[1, 0].set_ylabel('Roll Rate Body [deg/s]', color='white')
    axes[1, 0].grid(True)
    axes[1, 0].set_facecolor('black')
    axes[1, 0].tick_params(colors = 'white')

    # Pitch rate q^b_b/n
    axes[1, 1].plot(t_s, r2d*x[4,:], color='yellow')
    axes[1, 1].set_xlabel('Time [s]', color='white')
    axes[1, 1].set_ylabel('Pitch Rate Body [deg/s]', color='white')
    axes[1, 1].grid(True)
    axes[1, 1].set_facecolor('black')
    axes[1, 1].tick_params(colors = 'white')

    # Yaw rate r^b_b/n
    axes[1, 2].plot(t_s, r2d*x[5,:], color='yellow')
    axes[1, 2].set_xlabel('Time [s]', color='white')
    axes[1, 2].set_ylabel('Yaw Rate Body [deg/s]', color='white')
    axes[1, 2].grid(True)
    axes[1, 2].set_facecolor('black')
    axes[1, 2].tick_params(colors = 'white')

    plt.tight_layout()
    if sim_config.save_6dof_plot == "on":  
        plt.savefig(sim_config.save_6dof_plot_dir)
    plt.show(block=False)

    # B. Create figure of Euler angles
    fig, axes = plt.subplots(1, 3, figsize=(10, 3))
    fig.set_facecolor('black')  
    title_suffix = f"{sim_config.vmod['V_name']}"
    title_prefix = "Euler Angles (Attitude)\n"
    fig.suptitle(title_prefix + title_suffix, fontsize=14, fontweight='normal', color='yellow')

    # Roll angle, phi
    axes[0].plot(t_s, r2d*phi_2_rad, color='yellow')
    axes[0].set_xlabel('Time [s]', color='white')
    axes[0].set_ylabel('Roll Angle [deg]', color='white')
    axes[0].grid(True)
    axes[0].set_facecolor('black')
    axes[0].tick_params(colors = 'white')

    # Pitch angle, theta
    axes[1].plot(t_s, r2d*theta_2_rad, color='yellow')
    axes[1].set_xlabel('Time [s]', color='white')
    axes[1].set_ylabel('Pitch Angle [deg]', color='white')
    axes[1].grid(True)
    axes[1].set_facecolor('black')
    axes[1].tick_params(colors = 'white')

    # Yaw angle, theta
    axes[2].plot(t_s, r2d*psi_2_rad, color='yellow')
    axes[2].set_xlabel('Time [s]', color='white')
    axes[2].set_ylabel('Yaw Angle [deg]', color='white')
    axes[2].grid(True)
    axes[2].set_facecolor('black')
    axes[2].tick_params(colors = 'white')

    plt.tight_layout()
    if sim_config.save_euler_angle_plot == "on":  
        plt.savefig(sim_config.save_euler_angle_plot_dir)
    plt.show(block=False)

    # C. Create figure of aircraft positions and velocities resolved in the NED coordinate system
    fig, axes = plt.subplots(2, 3, figsize=(10, 6))
    fig.set_facecolor('black') 
    title_suffix = f"{sim_config.vmod['V_name']}"
    title_prefix = "Inertial Position and Velocity\n"
    fig.suptitle(title_prefix + title_suffix, fontsize=14, fontweight='normal', color='cyan')

    # North position p1^n_CM/T
    axes[0,0].plot(t_s, x[9,:], color='cyan')
    axes[0,0].set_xlabel('Time [s]', color='white')
    axes[0,0].set_ylabel('North Position [m]', color='white')
    if np.linalg.norm(x[9,:]) < 1e-5:
        axes[0,0].set_ylim(-0.05,0.05)
    axes[0,0].grid(True)
    axes[0,0].set_facecolor('black')
    axes[0,0].tick_params(colors = 'white')

    # East position p2^n_CM/T
    axes[0,1].plot(t_s, x[10,:], color='cyan')
    axes[0,1].set_xlabel('Time [s]', color='white')
    axes[0,1].set_ylabel('East Position [m]', color='white')
    if np.linalg.norm(x[10,:]) < 1e-5:
        axes[0,1].set_ylim(-0.05,0.05)
    axes[0,1].grid(True)
    axes[0,1].set_facecolor('black')
    axes[0,1].tick_params(colors = 'white')

    # Altitude
    axes[0,2].plot(t_s, -x[11,:], color='cyan')
    axes[0,2].set_xlabel('Time [s]', color='white')
    axes[0,2].set_ylabel('Altitude [m]', color='white')
    if np.linalg.norm(x[11,:]) < 1e-5:
        axes[0,2].set_ylim(-0.05,0.05)
    axes[0,2].grid(True)
    axes[0,2].set_facecolor('black')
    axes[0,2].tick_params(colors = 'white')

    # u_n_mps
    axes[1,0].plot(t_s, u_n_mps, color='cyan')
    axes[1,0].set_xlabel('Time [s]', color='white')
    axes[1,0].set_ylabel('u Velocity NED [ft/s]', color='white')
    if np.linalg.norm(u_n_mps) < 1e-5:
        axes[1,0].set_ylim(-0.05,0.05)
    axes[1,0].grid(True)
    axes[1,0].set_facecolor('black')
    axes[1,0].tick_params(colors = 'white')

    # v_n_mps
    axes[1,1].plot(t_s, v_n_mps, color='cyan')
    axes[1,1].set_xlabel('Time [s]', color='white')
    axes[1,1].set_ylabel('v Velocity NED [ft/s]', color='white')
    if np.linalg.norm(v_n_mps) < 1e-5:
        axes[1,1].set_ylim(-0.05,0.05)
    axes[1,1].grid(True)
    axes[1,1].set_facecolor('black')
    axes[1,1].tick_params(colors = 'white')

    # w_n_mps
    axes[1,2].plot(t_s, w_n_mps, color='cyan')
    axes[1,2].set_xlabel('Time [s]', color='white')
    axes[1,2].set_ylabel('w Velocity NED [ft/s]', color='white')
    if np.linalg.norm(w_n_mps) < 1e-5:
        axes[1,2].set_ylim(-0.05,0.05)
    axes[1,2].grid(True)
    axes[1,2].set_facecolor('black')
    axes[1,2].tick_params(colors = 'white')

    plt.tight_layout()
    if sim_config.save_position_plot == "on":  
        plt.savefig(sim_config.save_position_plot_dir)
    plt.show(block=False)

    # D. Create figure of air data
    fig, axes = plt.subplots(1, 4, figsize=(12, 4))
    fig.set_facecolor('black')  
    title_suffix = f"{sim_config.vmod['V_name']}"
    title_prefix = "Air Data\n"
    fig.suptitle(title_prefix + title_suffix, fontsize=14, fontweight='normal', color='magenta')

    # Angle of attack
    axes[0].plot(t_s, Alpha_rad*180/3.14, color='magenta')
    axes[0].set_xlabel('Time [s]', color='white')
    axes[0].set_ylabel('Angle of Attack [deg]', color='white')
    axes[0].set_ylim(-100,100)
    axes[0].grid(True)
    axes[0].set_facecolor('black')
    axes[0].tick_params(colors = 'white')

    # Angle of side slip
    axes[1].plot(t_s, Beta_rad*180/3.14, color='magenta')
    axes[1].set_xlabel('Time [s]', color='white')
    axes[1].set_ylabel('Angle of Side Slip [deg]', color='white')
    axes[1].set_ylim(-90,90)
    axes[1].grid(True)
    axes[1].set_facecolor('black')
    axes[1].tick_params(colors = 'white')

    # Mach
    axes[2].plot(t_s, Mach, color='magenta')
    axes[2].set_xlabel('Time [s]', color='white')
    axes[2].set_ylabel('Mach Number', color='white')
    axes[2].grid(True)
    axes[2].set_facecolor('black')
    axes[2].tick_params(colors = 'white')

    # Altitude vs Mach
    axes[3].plot(Mach, -x[11,:], color='magenta')
    axes[3].set_xlabel('Mach Number []', color='white')
    axes[3].set_ylabel('Altitude [m]', color='white')
    axes[3].grid(True)
    axes[3].set_facecolor('black')
    axes[3].tick_params(colors = 'white')

    plt.tight_layout()
    if sim_config.save_air_data_plot == "on":  
        plt.savefig(
            sim_config.save_air_data_plot_dir)
    plt.show()
