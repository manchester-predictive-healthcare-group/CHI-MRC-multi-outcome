cd "/mnt/bmh01-rds/mrc-multi-outcome/Project_6/code/"
module load apps/gcc/R/4.0.2-gcc830
nohup Rscript p3.3.1_large_sample_analysis_moderate.R M1C1 200000 20 > p3.3.1_large_sample_analysis_moderate_200000_M1C1_npctls20.out &
nohup Rscript p3.3.1_large_sample_analysis_moderate.R M1C2 200000 20 > p3.3.1_large_sample_analysis_moderate_200000_M1C2_npctls20.out &
nohup Rscript p3.3.1_large_sample_analysis_moderate.R M1C3 200000 20 > p3.3.1_large_sample_analysis_moderate_200000_M1C3_npctls20.out &