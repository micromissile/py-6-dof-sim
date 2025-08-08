**Python Flight Simulation - Version 0.1.1**

Developed in Python 3.11.9

PACKAGE REQUIREMENTS:
numpy
subprocess
matplotlib
ussa1976
flightgear_python (only you wish to animate flight in Flight Gear)

INSTRUCTIONS
To run PFS with existing jobs:
1. Make Python_Flight_Simulation your present working directory.
2. In run.py select the job you would like to run, such as run_job_4(), under Step 1.
3. To visualize flight in Flight Gear, uncomment subprocess.run(['python', job_path], and job_path)
    a. Open Flight Gear, select the appropriate aircraft, and hit "Fly!"
    b. The aircraft will be shown in its environment, select your desired view
4. In your terminal window type: python .\run.py (the simulation runs, plots are generated)
6. Close all plots and the visualization automatically runs in Flight Gear.

To run PFS 6DOF OR Flight Gear separately:
1. 6DOF only: Follow steps 1., 2., and 4. (leave subprocess.run(['python', job_path]) commented)
2. Flight Gear only: Comment run_job_X() in run.py, Uncomment subprocess.run(['python', job_path] 
   and job path).

To run a different job:
1. Navigate to the job folder you wish to modify (or copy and save it as a template for a new job)
2. In sim_config.py, update the initial conditions, control input, and other options.
3. If you created a new job subpackage
    a. Add import of this subpackage in run.py (similar to other packages)
    b. If you intend to visualize flight, add job_path to stored data (data store option in sim_config)
4. Follow the above instructions to run execute the job

PFS UPDATE:
This simulation is reorganized from previous versions into an improved package structure.

Python Flight Simulation/
├── __init__.py
├── control
├── jobs 
├── numerical_integrators
├── tools
├── vehicle_models
└── run.py

All simulation cases are now classified as jobs, stored in the jobs subpackage 
as uniquely named subpackages themselves. For example, the sphere simulation shown
in Section 1.5 is now captured as a job with the subpackage name,
'S_1_5_1_Sphere_Verification'. Navigating into this folder, you will see several items.
 
S_1_5_1_Sphere_Verification/
├── __init__.py
├── flat_earth_eom.py
├── sim_config.py
├── visualization.py
├── output_data/
│   ├── Atmos_01_Python_6DOF.npy
└── plots/
    ├── S1p5_Atmos01_6dof_Plot.png
    ├── S1p5_Atmos01_Air_Data_Plot.png
    ├── ...
└── verification/
    ├── Atmos_01_sim_01.csv
    ├── ...
└── main_program.py
 
The contents of each job subpackage have common core elements. However, the sphere and 
brick verification cases have the verification folder to compare simulation results to the 
NASA verification cases. Of course, this does not apply to all jobs, and so there will be 
custom subpackages and modules in each job folder to suit the jobs specific purpose.

In each job, there is the main_program (job driver), flat_earth_eom (customized for
each job e.g. Euler angles versus quaternions), sim_config (initialization data, plot, 
selections, etc.), and visualization (a Flight Gear interface for immediate visualization
of the flight dynamics). 

Note that for each job there is a specific equation of motion (EOM) module. The intent
is for the user to copy the closest EOM script for their application to make the 
necessary modifications for their purposes. It can become tedious, time-consuming, and
error-prone to make a one-sizes fits all EOM file, particularly when updates are made
and backward compatibility needs to be verified. Since each job will have its own 
EOM script, the need to check backward compatibility is completely avoided.

This rationale extends to all common files of existing jobs. They are to be used as 
templates to build from and then added to with new jobs, thereby expanding the breadth and 
utility of the simulation. By setting up the code in this way, the goal is to create an 
immediately usable archive of flight simulations while allowing the user to customize
the flight simulation for their own investigations through the creation of new jobs.

In future versions of this simulation, a wizard tool or a more generalized template is 
envisioned that will allow the user to create a job package for running simulations 
without having to review the contents of other jobs. For now, the latest job is the 
recommended template.

B.D.
October 13th, 2024