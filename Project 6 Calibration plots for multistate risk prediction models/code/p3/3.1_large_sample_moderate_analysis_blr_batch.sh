module load libs/gcc/gsl/2.5 apps/gcc/R/4.3.1 libs/gcc/nlopt/2.7.1 compilers/gcc/11.2.0
nohup Rscript 3.1_large_sample_moderate_analysis_blr.R C1 > 3.1_large_sample_moderate_analysis_blr_C1.out &
nohup Rscript 3.1_large_sample_moderate_analysis_blr.R C2 > 3.1_large_sample_moderate_analysis_blr_C2.out &
nohup Rscript 3.1_large_sample_moderate_analysis_blr.R C3 > 3.1_large_sample_moderate_analysis_blr_C3.out &