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

The files p1_run_sim_* in both simulation studies are designed to be run in parallel. Changes the 'set' argument, read in at the command line, will results in the creation of different scenarios. The simulation study 1 files were ran for set = 1:300, to create 15,000 scenarios. The simulation study 2 files were ran for set = 1:500, to create 12,500 scenarios.
