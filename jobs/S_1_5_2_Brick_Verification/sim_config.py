# sim_config.py - This is the configuration file for the 6DOF flight simulation. It 
# contains data that defines the simulation including vehicle, initial conditions, 
# duration of simulation, integrator time step size, data and plot saving options, 
# and their directories.
#
# Written by Ben Dickinson Oct 5th, 2024
#
#----------------------------------------------------------------------------------
from vehicle_models.brick import bricks

# Define vehicle model
vmod = bricks.NASA_Atmos02_Brick()

# Empty control model
cmod = {}

# Options to save and directories of data
#save_data = 'on'
save_data = 'off'
save_dir  = './jobs/S_1_5_2_Brick_Verification/output_data/Atmos_02_Python_6DOF.npy'

#save_6dof_plot = 'on'
save_6dof_plot = 'off'
name_6dof_plot = 'S1p5_' + vmod["short_name"] + '_6dof_Plot.png'
save_6dof_plot_dir = './jobs/S_1_5_2_Brick_Verification/plots/' + name_6dof_plot

#save_euler_angle_plot = 'on'
save_euler_angle_plot = 'off'
name_euler_angle_plot = 'S1p5_' + vmod["short_name"] + '_Euler_Plot.png'
save_euler_angle_plot_dir = './jobs/S_1_5_2_Brick_Verification/plots/' + name_euler_angle_plot

#save_position_plot = 'on'
save_position_plot = 'off'
name_position_plot = 'S1p5_' + vmod["short_name"] + '_Position_Plot.png'
save_position_plot_dir = './jobs/S_1_5_2_Brick_Verification/plots/' + name_position_plot

#save_air_data_plot = 'on'
save_air_data_plot = 'off'
name_air_data_plot = 'S1p5_' + vmod["short_name"] + '_Air_Data_Plot.png'
save_air_data_plot_dir = './jobs/S_1_5_2_Brick_Verification/plots/' + name_air_data_plot


