module load libs/gcc/gsl/2.5 apps/gcc/R/4.3.1 libs/gcc/nlopt/2.7.1 compilers/gcc/11.2.0
nohup Rscript 3.2_large_sample_moderate_analysis_sens.R C1 > 3.2_large_sample_moderate_analysis_sens_C1.out &
nohup Rscript 3.2_large_sample_moderate_analysis_sens.R C2 > 3.2_large_sample_moderate_analysis_sens_C2.out &
nohup Rscript 3.2_large_sample_moderate_analysis_sens.R C3 > 3.2_large_sample_moderate_analysis_sens_C3.out &