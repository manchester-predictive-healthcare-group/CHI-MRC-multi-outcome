module load libs/gcc/gsl/2.5 apps/gcc/R/4.3.1 libs/gcc/nlopt/2.7.1 compilers/gcc/11.2.0
nohup Rscript 2.1_large_sample_mean_analysis.R C1 > 2.1_large_sample_mean_analysis_C1.out &
nohup Rscript 2.1_large_sample_mean_analysis.R C2 > 2.1_large_sample_mean_analysis_C2.out &
nohup Rscript 2.1_large_sample_mean_analysis.R C3 > 2.1_large_sample_mean_analysis_C3.out &