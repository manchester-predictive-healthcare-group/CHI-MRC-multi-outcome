# Directory structure #

All code is ran from the 'code' directory. All data files are stored in the 'data' directory. All figures are saved into the 'Figures' directory. 

The .Rmd file to generate the supplementary figures is also stored in the Figures directory.


# Code structure #

All .R files have an initial line which clears the workspace and then sets the working directory to be the root directory where 'code', 'data' and 'Figures' are stored. All code is run relative to this root directory. Alternatively an R project can be launced in this root directory.

All files should be ran in numerical order. A number of .R files require input at the command line (normally defining scenario, sample size or number of percentiles). For these files, .sh files or .qsub files are provided to run the code with the appropriate input.

Scenario "M1C1", "M1C2" and "M1C3" correspond to NIC, WIC and SIC from the manuscript. 

Suffixes "est1", "est2" and "est3" correspond to the perfectly predicting, over predicting and under predicting estimated transition probabilities.

This should help navigate the scenarios and the Figure files.

p1.X files run the clinical example. This code assumes the CPRD dataset has already been derived. Code for deriving the CPRD dataset is available in another GitHub directory.

p2.X files generate the data for the simulation.

p3.X files run the large validation sample analysis.

p4.X files run the small validation sample analysis.

All code is commented. If any questions please contact me.