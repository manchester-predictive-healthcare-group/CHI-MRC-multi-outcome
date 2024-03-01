module load libs/gcc/gsl/2.5 apps/gcc/R/4.3.1 libs/gcc/nlopt/2.7.1 compilers/gcc/11.2.0
nohup Rscript 1.1_large_sample_prep_data.R C1  > 1.1_large_sample_prep_data_C1.out &
nohup Rscript 1.1_large_sample_prep_data.R C2  > 1.1_large_sample_prep_data_C2.out &
nohup Rscript 1.1_large_sample_prep_data.R C3  > 1.1_large_sample_prep_data_C3.out &
