# run.py is the main driver for this flight simulation. This run.py script was
# introduced on October 12th, 2024. The purpose is to allow all previously created
# simulation cases to be run readily and easily by the user. The user essentially
# comments/uncomments the job the desire to run, and runs this driver with 
# Spectrum Flight Simulation as the present working directory.
#
# History
# Written by Ben Dickinson, October 12th, 2024.
# www.LearnGandC.com
# Available at: https://www.patreon.com/user?u=86359827
#---------------------------------------------------------------------------------------
import subprocess
from time import time

from jobs.S_1_5_1_Sphere_Verification.main_program                   import run_job as run_job_1
from jobs.S_1_5_2_Brick_Verification.main_program                    import run_job as run_job_2
from jobs.S_1_5_3_Brick_Verification.main_program                    import run_job as run_job_3
from jobs.S_2_2_3_X15_Aileron_Doublet.main_program                   import run_job as run_job_4
from jobs.S_2_2_3_X15_Elevator_Doublet.main_program                  import run_job as run_job_5
from jobs.S_2_2_3_X15_Rudder_Doublet.main_program                    import run_job as run_job_6
from jobs.S_2_2_4_X15_Glide_Test.main_program                        import run_job as run_job_7
from jobs.S_3_1_1_X15_SAS_I.main_program                             import run_job as run_job_8
from jobs.S_3_2_1_X15_Full_Aero_Database.main_program                import run_job as run_job_9
from jobs.S_3_2_1_X15_Full_Aero_Database_Lateral_Test.main_program   import run_job as run_job_10
from jobs.S_3_2_2_X15_Dutch_Roll.main_program                        import run_job as run_job_11
from jobs.Rocket_Dev_Horizontal.main_program                         import run_job as run_job_12
from jobs.Rocket_Dev_Vertical_100m.main_program                      import run_job as run_job_13


# Step 1. Run the program
start_time = time() 
run_job_12()
# run_job_13()
end_time = time()
print(f"Total simulation time: {end_time - start_time:.4f} seconds")
#
# This is the primary flight simulation driver. It exists in each of the jobs
# and is custom to each of the jobs. Running main_program from the appropriate job 
# will load the vehicle model, simulation configuration (sim_config.py), other necessary
# data, then integrate the 6DOF simulation forward in time, and plot the results.

# Step 2. Visualize the simulation
#
# Specify the job path to the FlightGear-Python interface script, visualization.py.
#job_path = r'./jobs/S_1_5_1_Sphere_Verification/visualization.py'
#job_path = r'./jobs/S_1_5_2_Brick_Verification/visualization.py'
#job_path = r'./jobs/S_1_5_3_Brick_Verification/visualization.py'
#job_path = r'./jobs/S_2_2_3_X15_Aileron_Doublet/visualization.py'
#job_path = r'./jobs/S_2_2_3_X15_Elevator_Doublet/visualization.py'
#job_path = r'./jobs/S_2_2_3_X15_Rudder_Doublet/visualization.py'
#job_path = r'./jobs/S_2_2_4_X15_Glide_Test/visualization.py'
#job_path = r'./jobs/S_3_1_1_X15_SAS_I/visualization.py'
#job_path = r'./jobs/S_3_2_1_X15_Full_Aero_Database/visualization.py'
#job_path = r'./jobs/S_3_2_1_X15_Full_Aero_Database_Lateral_Test/visualization.py'
#job_path = r'./jobs/S_3_2_2_X15_Dutch_Roll/visualization.py'
#
# Run visualization.py
#subprocess.run(['python', job_path])
#
# This allows visualization of the flight dynamics from the simulation data run above.
# Existing job packages will be set up so that visualization.py can be run. This will
# require Flight Gear on your local machine. Before attempting to run the visualization, 
# ensure flight gear is open, the appropriate settings have been manually entered in 
# flight gear according to the flight-gear python package, the corresponding aircraft 
# has been selected in Flight Gear, and the Flight Gear is running (the "Fly!" button on 
# the bottom left of the main menu has been pressed).
#
# The data imported into flight gear for visualization already exists in each of the jobs
# folders, within output data. visualization.py specifically calls this data to produce
# the visualization. If you save new data in output_data (as specified in sim_config), 
# that data can be visualized so long as the name of the data is updated in the 
# visualization.py module.
