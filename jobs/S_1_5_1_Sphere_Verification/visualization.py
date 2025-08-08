"""
Simple Flight Dynamics Model (FDM) example of a tumbling brick falling through the air.

History:
Adapted from "Simple FDM Loop" example by Ben Dickinson
6/24/2024
"""
import time
import numpy as np
from flightgear_python.fg_if import FDMConnection

def fdm_callback(fdm_data, event_pipe):
    
    if event_pipe.child_poll():
        
        # Unpack tuple
        alt_m_child, phi_rad_child, theta_rad_child, psi_rad_child, \
        u_b_ftps_child, v_b_ftps_child, w_b_ftps_child, alpha_rad_child, \
            beta_rad_child, long_rad_child, lat_rad_child, = event_pipe.child_recv()  
        
        # Set only the data that we need to (we can force our own values)
        fdm_data['alt_m']       = alt_m_child  
        fdm_data['phi_rad']     = phi_rad_child
        fdm_data['theta_rad']   = theta_rad_child
        fdm_data['psi_rad']     = psi_rad_child
        fdm_data['v_body_u']    = u_b_ftps_child
        fdm_data['v_body_v']    = v_b_ftps_child
        fdm_data['v_body_w']    = w_b_ftps_child
        fdm_data['alpha_rad']   = alpha_rad_child
        fdm_data['beta_rad']    = beta_rad_child
        fdm_data['lon_rad']     = long_rad_child
        fdm_data['lat_rad']     = lat_rad_child
        
    # Return the whole structure
    return fdm_data  

"""
Start FlightGear with: 
`--native-fdm=socket,out,30,localhost,5501,udp --native-fdm=socket,in,30,localhost,5502,udp`
(you probably also want `--fdm=null` and `--max-fps=30` to stop the simulation fighting with
these external commands)
"""
if __name__ == '__main__':  # NOTE: This is REQUIRED on Windows!
    
    fdm_conn = FDMConnection()
    fdm_event_pipe = fdm_conn.connect_rx('localhost', 5501, fdm_callback)
    fdm_conn.connect_tx('localhost', 5502)
    
    # Start the FDM RX/TX loop
    fdm_conn.start()  
    
    # Get Python 6-DOF simulation data
    data_01_pysim  = np.load('./jobs/S_1_5_1_Sphere_Verification/output_data/Atmos_01_Python_6DOF.npy')
    
    # Get time and other variables
    t_s   = data_01_pysim[:,0]; nt_s = t_s.size

    i = 0
    while True:
        
        # Increment time step counter
        i += 1
        
        # Get present altitude
        v_body_u_parent     = data_01_pysim[i,1]*3.28 # (converts m/s to ft/s)
        v_body_v_parent     = data_01_pysim[i,2]*3.28
        v_body_w_parent     = data_01_pysim[i,3]*3.28
        phi_rad_parent      = data_01_pysim[i,7]
        theta_rad_parent    = data_01_pysim[i,8]
        psi_rad_parent      = data_01_pysim[i,9]
        alt_m_parent        = data_01_pysim[i,13]
        alpha_rad_parent    = data_01_pysim[i,17]
        beta_rad_parent     = data_01_pysim[i,18]
        
        # Stromboli volcano 
        #lat_rad_parent  = 38.793889*0.01745329 # (converts deg to rad)
        lat_rad_parent  = 38.7971*0.01745329 # (converts deg to rad)
        #long_rad_parent = 15.211111*0.01745329
        long_rad_parent = 15.214*0.01745329
        
        # Katla volcano 
        # lat_rad_parent  = 63.633333*0.01745329 # (converts deg to rad)
        # long_rad_parent = -19.05*0.01745329
        
        # Send tuple (could also do `fdm_conn.event_pipe.parent_send` so you just need to pass around `fdm_conn`)
        fdm_event_pipe.parent_send((alt_m_parent, phi_rad_parent, theta_rad_parent, psi_rad_parent, \
            v_body_u_parent, v_body_v_parent, v_body_w_parent, alpha_rad_parent, beta_rad_parent, \
                long_rad_parent, lat_rad_parent, ))  
        
        if i == nt_s-1:
            print('\nVisualization completed.\n')
            break
        
        time.sleep(0.007)