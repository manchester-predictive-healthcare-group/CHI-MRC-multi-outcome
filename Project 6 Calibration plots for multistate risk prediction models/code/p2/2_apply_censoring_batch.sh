cd "/mnt/bmh01-rds/mrc-multi-outcome/Project_6/code/p2/"
module load libs/gcc/gsl/2.5 apps/gcc/R/4.3.1 libs/gcc/nlopt/2.7.1 compilers/gcc/8.3.0
nohup Rscript 2_apply_censoring.R C1 > 2_apply_censoring_C1.out &
nohup Rscript 2_apply_censoring.R C2 > 2_apply_censoring_C2.out &
nohup Rscript 2_apply_censoring.R C3 > 2_apply_censoring_C3.out &