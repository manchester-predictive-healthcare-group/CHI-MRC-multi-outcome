The simulations are designed to be run relative to some parent directory. Place the four folders in the parent directory, and either place an R project in this directory, or manually set the working directory at the top of each piece of code in the 'code' folder.

Folders

Folder 'code' contains code to run to produce results.
Folder 'data' contains the results from the simulations (i.e. the table with all the differences scenarios, input data, results, etc). The 'data' folder is also where all intermediate files get saved when running the programs in the 'code' folder
Folder 'figures' contains all the figures from manuscript and supplementary material
Folder 'R' contains all the functions which are sourced and used in the .R files in the 'code' directory.

Subdirectories.

The 'code' folder contains three subdirectories:
'p1_sim_study1' contains the files to implement simulation study 1.
'p2_sim_study2' contains the files to implement simulation study 2.
'p3_applied_examples' contains the files to implement the applied examples.
'p4_sim_study3' contains the files to implement simulation study 1, using the alternative data generating mechanism outlined in section XXXX.
'p5_sim_study4' contains the files to implement simulation study 2, using the alternative data generating mechanism outlined in section XXXX..

The files within each of these folders has a similar structure and naming convention:

The files p1_run_sim_*.R in both simulation studies are designed to be run in parallel. Changing the 'set' argument, read in at the command line, will change the seed and result in the creation of different scenarios. If "corr" is present in the .R filename, this will run scenarios where there is correlation/covariance between the predictors. If "nocorr" is present in the .R filename, this runs scenarios with no correlation/covariance.

The files p2_combine_sim_data*.R will combine the simulation results into a single table, containing all the metadata about each scenario, and the results from each scenario. These tables are what all the plots and tables are derived from. These have also been placed into the folder "metadata" for easy access.

The files p3_*.R and p3_*.Rmd and above (p4, p5), produce the plots and tables presented in the manuscript and supplementary material. 