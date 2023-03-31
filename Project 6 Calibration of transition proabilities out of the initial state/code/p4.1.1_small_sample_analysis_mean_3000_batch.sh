cd "/mnt/bmh01-rds/mrc-multi-outcome/Project_6/code/"
module load apps/gcc/R/4.0.2-gcc830
nohup Rscript p4.1.1_small_sample_analysis_mean.R M1C1 3000 5 > p4.1.1_small_sample_analysis_mean_M1C1_3000_npctls5.out &
nohup Rscript p4.1.1_small_sample_analysis_mean.R M1C2 3000 5 > p4.1.1_small_sample_analysis_mean_M1C2_3000_npctls5.out &
nohup Rscript p4.1.1_small_sample_analysis_mean.R M1C3 3000 5 > p4.1.1_small_sample_analysis_mean_M1C3_3000_npctls5.out &

nohup Rscript p4.1.1_small_sample_analysis_mean.R M1C1 3000 10 > p4.1.1_small_sample_analysis_mean_M1C1_3000_npctls10.out &
nohup Rscript p4.1.1_small_sample_analysis_mean.R M1C2 3000 10 > p4.1.1_small_sample_analysis_mean_M1C2_3000_npctls10.out &
nohup Rscript p4.1.1_small_sample_analysis_mean.R M1C3 3000 10 > p4.1.1_small_sample_analysis_mean_M1C3_3000_npctls10.out &

nohup Rscript p4.1.1_small_sample_analysis_mean.R M1C1 3000 20 > p4.1.1_small_sample_analysis_mean_M1C1_3000_npctls20.out &
nohup Rscript p4.1.1_small_sample_analysis_mean.R M1C2 3000 20 > p4.1.1_small_sample_analysis_mean_M1C2_3000_npctls20.out &
nohup Rscript p4.1.1_small_sample_analysis_mean.R M1C3 3000 20 > p4.1.1_small_sample_analysis_mean_M1C3_3000_npctls20.out &

