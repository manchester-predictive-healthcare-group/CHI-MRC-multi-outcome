module load libs/gcc/gsl/2.5 apps/gcc/R/4.3.1 libs/gcc/nlopt/2.7.1 compilers/gcc/11.2.0
nohup Rscript 2.2_large_sample_mean_analysis_sens.R C1 > 2.2_large_sample_mean_analysis_sens_C1.out &
nohup Rscript 2.2_large_sample_mean_analysis_sens.R C2 > 2.2_large_sample_mean_analysis_sens_C2.out &
nohup Rscript 2.2_large_sample_mean_analysis_sens.R C3 > 2.2_large_sample_mean_analysis_sens_C3.out &