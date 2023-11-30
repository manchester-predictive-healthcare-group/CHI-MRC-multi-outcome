cd "/mnt/bmh01-rds/mrc-multi-outcome/Project_6/code/"
module load apps/gcc/R/4.0.2-gcc830
nohup Rscript p4.1.3_plots_small_sample_analysis_mean.R 1500 5 > p4.1.3_plots_small_sample_analysis_mean_N1500_npctls5.out &
nohup Rscript p4.1.3_plots_small_sample_analysis_mean.R 1500 10 > p4.1.3_plots_small_sample_analysis_mean_N1500_npctls10.out &
nohup Rscript p4.1.3_plots_small_sample_analysis_mean.R 1500 20 > p4.1.3_plots_small_sample_analysis_mean_N1500_npctls20.out &
nohup Rscript p4.1.3_plots_small_sample_analysis_mean.R 3000 5 > p4.1.3_plots_small_sample_analysis_mean_N3000_npctls5.out &
nohup Rscript p4.1.3_plots_small_sample_analysis_mean.R 3000 10 > p4.1.3_plots_small_sample_analysis_mean_N3000_npctls10.out &
nohup Rscript p4.1.3_plots_small_sample_analysis_mean.R 3000 20 > p4.1.3_plots_small_sample_analysis_mean_N3000_npctls20.out &