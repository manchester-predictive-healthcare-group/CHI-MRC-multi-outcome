Code for manuscript

All R code directories are relative to a root directory, which is the data and materials folder. Either start an R project file where the root directory is positioned here, or change the working directory in each R file. The line to do this is located at the top of each piece of code. It should be set to the directory which contains the folders /code, /data and /figures.

Code should be run in order of prefix's, p0, p1, p2, p3, p4.

p0 are preliminary pieces of code including sample size calculations and making a record of input parameters

p1 run the large sample simulation

p2 run the small sample simulation. This code was set up to run on the computational shared facility at university of manchester and requires three arguments to be read in at the command line (seed for simulation, seed for coefficients, sample size of development cohorts.

p3 processes results of small sample simulation

p4 produces all table sand figures

