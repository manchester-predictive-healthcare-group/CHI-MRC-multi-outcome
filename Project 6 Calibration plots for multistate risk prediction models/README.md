# Directory structure #

- All code is ran from the 'code' directory. 
- All data files are stored in the 'data' directory. 
- All figures are saved into the 'Figures' directory.
- The .Rmd file to generate the html file with all Figures is also stored in the Figures directory.
- Functions used throughout the code are stored in the z_functions.R file. This is sourced is the start of every piece of code.
- Required packages are loaded by sourcing the z_load_packages.R file at the start of every piece of code.
- Items relevant to the manuscript are contained in the 'Manuscript files' section. This contains a html file with all Figures produced from the simulation
(see document aa_supplementary_material_part2.html). Figures included in the manuscript, and a validation of the formulae for estimating the transition probabilities.
- Data and individual figures have not been provided to reduce disk space. These folders are there as placeholders, and required in order to run the code provided.

# Code structure #

All .R files have an initial line which clears the workspace and then sets the working directory to be the root directory where 'code', 'data' and 'Figures' are stored. All code is run relative to this root directory.
This is the only change required to run the simulation on your own computer. Alternatively an R project can be launced in this root directory.

All files should be ran in numerical order. A number of .R files require input at the command line (normally defining scenario, sample size or number of percentiles). 
For these files, .sh files or .qsub files are provided to run the code with the appropriate input. This was done to achieve parallelisation.

Scenario "C1", "C2" and "C3" correspond to RC, WAC and SAC from the manuscript. 

## Directory p1

Files in this directory run the clinical example. This code assumes the CPRD dataset has already been derived. Code for deriving the CPRD dataset is available in another GitHub directory, however access to the underlying data cannot be provided.

Note that because calculating the pseudo-values is very computationally expensive, this is done outside of the calibmsm package to enable parallelisation. These pseudo-values
are then fed into the calibmsm package when assessing calibration using the pseudo-value approach.

1_ce_fit_csh_msm.R - fits the multistate models for the clinical example.
2_ce_calc_probtrans_noloop.R - estimates transition probabilities for individuals in the validation cohort.
3.1_ce_pv_prep.R - preps data for estimating the pseudo-values.
3.2_ce_pv_calc.R - calculates pseudo-vlaues.
3.3_ce_pv_combine.R - combines the estimated pseudo-values into a single dataframe.
4_ce_assess_calib.R - assesses calibration.
5_ce_create_plots.R - creates plots.


## Directory p2

Files in this directory prepare data for the simulation.

0_input_parameters.R - is sourced at the start of other files. This loads the input parameters for the simulation scenarios.
1_generate_population_data.R - generates the data for the simulation.
2_apply_censoring.R - applies the censoring and converts data to 'msdata' format.
3_validate_true_transition_probabilities.R and .Rmd - validate the formulae for estimating the transition probabilities.

## Directory p3

Files in this directory run the large sample simulation. 

Note that because calculating the pseudo-values is very computationally expensive, this is done outside of the calibmsm package to enable parallelisation. These pseudo-values
are then fed into the calibmsm package when assessing calibration using the pseudo-value approach.

1.1_large_sample_prep_data.R - specific data preperation steps for the large sample simulation.
1.2_large_sample_prep_data_pv.R - specific data preperation steps for the pseudo-values.
1.3_large_sample_calc_pv.R - calculates the pseudo-values.
1.4_large_sample_combine_pv.R - combines the pseudo-values into a single dataframe.
1.5_large_sample_prep_data_pv_sens.R - specific data preperation steps for calculating the pseudo-values in the sensitivity analysis.
1.6_large_sample_calc_pv_sens.R - calculates the pseudo-values for the sensitivity analysis.
1.7_large_sample_combie_pv_sens.R - combines pseudo-values for the sensitivity analysis into a single dataframe.
2.1_large_sample_mean_analysis.R - assess mean calibration.
2.2_large_sample_mean_analysis_sens.R - sensitivity analysis for mean calibration.
2.3_large_sample_mean_create_plots.R - create plots for mean calibration.
3.1_large_sample_moderate_analysis_blr.R - assess moderate calibration with BLR-IPCW
3.1_large_sample_moderate_analysis_mlr.R - assess moderate calibration with MLR-IPCW
3.1_large_sample_moderate_analysis_pv.R - assess moderate calibration with pseudo-value approach.
3.2_large_sample_moderate_analysis_sens.R - sensitivity analysis for moderate calibration with BLR-IPCW and MLR-IPCW.
3.2_large_sample_moderate_analysis_sens_pv.R - sensitivity analysis for moderate calibration with pseudo-value approach.
3.3_large_sample_moderate_estimate_true_calib_loess.R - estimate true calibration curve using loess smoother.
3.3_large_sample_moderate_estimate_true_calib_rcs.R - estimate true calibration curve using restricted cubic splines.
3.3_large_sample_moderate_estimate_true_calib_compare.R - compare true calibration curves estimated using loess and rcs.
3.4_large_sample_create_plots - create plots for moderate calibration.
3.5_large_sample_create_plots_sens - create plots for moderate calibration sensitivity analysis.


## Directory p4

Files in this directory run the large sample simulation. 

1.1_small_sample_prep_data.R - specific data preperation steps for the smallsample simulation.
2.1_small_sample_analysis.R - small sample analysis for mean calibration.
2.1_small_sample_analysis_moderate.R - small sample analysis for moderate calibration.
2.2_small_sample_analysis_sens.R - sensitivity analysis for small sample analysis.
2.3_small_sample_combine_data.R - combine data from small sample analysis for mean calibration.
2.3_small_sample_combine_data_moderate.R - combine data from small sample analysis for moderate calibration.
2.3_small_sample_combine_data_sens.R - combine data from sensitivity analysis for small sample analysis.
2.4_small_sample_create_plots.R - create plots for mean calibration and sensitivity analysis.
2.4_small_sample_create_plots_moderate.R - create plots for moderate calibration. 

All code is commented. If any questions please contact me.