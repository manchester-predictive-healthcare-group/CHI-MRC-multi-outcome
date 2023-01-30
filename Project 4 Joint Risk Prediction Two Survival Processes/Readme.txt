-----------------------------------------------------------------------------------------------
Overview of directories
---------------------------------------------------------------------------------------------

There are three directories in the directory “Project 4”: “code”, “data” and “figures”.

- The “code” directory contains all the code that needs to be run for the simulation and clinical example.
- The “data” folder is empty, but is where all the output would be saved to.
- The “figures” folder is where all the figures would be saved too. The markdown file to create Appendix S2 is also stored here.

For the simulation, all code is ran with working directories relative to the “Project 4” folder. There is a line of code at the top of each .R file which will set the working directory, make sure this points to the “project 4” folder, wherever this is stored on your computer. This is the only change that needs to be made to the code. Alternatively, these “setwd” commands could be removed, and make the “Project 4” the home directory for R, or place an .Rproj file into this directory to make it the home for an rstudio project. Realistically to run the simulation, parallelisation is required and we used a variety of computational environments, and therefore found it more appropriate to home one line of code at the top of each program, to set the working directory manually.

For the clinical example, this requires access to CPRD data wich we cannot share. This data was stored in directories stored in the same parent directory of “Project 4”. The .R programs to run the clinical example therefore set the working directory to be one level up from the “Project 4” folder.

-----------------------------------------------------------------------------------------------
Overview of code
-----------------------------------------------------------------------------------------------

The code is organised with different prefixes and should be run in that order (although p4 – p6 is the clinical example, and can be ran separately from p1 – p3). All functions beginning with “sim_function_” contains the functions which are called in when running the simulation or clinical example.

p0: Initial miscellaneous programs to calculate minimum sample size, input parameters for the simulation, and mean risks in the CPRD cohort (used to help decide on mean risks for the simulation scenarios).
p1: Run the simulation. The .qsub files can be used to run the .R files in parallel with the required input parameters.
p2: Combine the results from the simulation
p3: Produce tables and plots
p4: Program to load data for clinical example (is called upon in p5 programs)
p5: Fit the models for the clinical example. The predicted risks for the multistate model are generated using the .qsub file.
p6: Produce calibration plots, traceplots and tables for the clinical example.
p9: Report event rates as requested in peer review.
