cd "/mnt/bmh01-rds/mrc-multi-outcome/Project_6/code/"
module load apps/gcc/R/4.0.2-gcc830
nohup Rscript p3.2.1_large_sample_analysis_mean.R M1C1 200000 > p3.2.1_large_sample_analysis_mean_200000_M1C1.out &
nohup Rscript p3.2.1_large_sample_analysis_mean.R M1C2 200000 > p3.2.1_large_sample_analysis_mean_200000_M1C2.out &
nohup Rscript p3.2.1_large_sample_analysis_mean.R M1C3 200000 > p3.2.1_large_sample_analysis_mean_200000_M1C3.out &