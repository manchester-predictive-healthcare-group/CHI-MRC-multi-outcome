module load libs/gcc/gsl/2.5 apps/gcc/R/4.3.1 libs/gcc/nlopt/2.7.1 compilers/gcc/11.2.0
nohup Rscript 1.7_large_sample_combine_pv_sens.R C1 > 1.7_large_sample_combine_pv_sens_C1.out &
nohup Rscript 1.7_large_sample_combine_pv_sens.R C2 > 1.7_large_sample_combine_pv_sens_C2.out &
nohup Rscript 1.7_large_sample_combine_pv_sens.R C3 > 1.7_large_sample_combine_pv_sens_C3.out &