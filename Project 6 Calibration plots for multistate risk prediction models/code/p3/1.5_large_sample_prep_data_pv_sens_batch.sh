cd "/mnt/bmh01-rds/mrc-multi-outcome/Project_6/code/p3"
module load libs/gcc/gsl/2.5 apps/gcc/R/4.3.1 libs/gcc/nlopt/2.7.1 compilers/gcc/11.2.0
nohup Rscript 1.5_large_sample_prep_data_pv_sens.R C1 > 1.5_large_sample_prep_data_pv_sens_C1.out &
nohup Rscript 1.5_large_sample_prep_data_pv_sens.R C2 > 1.5_large_sample_prep_data_pv_sens_C2.out &
nohup Rscript 1.5_large_sample_prep_data_pv_sens.R C3 > 1.5_large_sample_prep_data_pv_sens_C3.out &