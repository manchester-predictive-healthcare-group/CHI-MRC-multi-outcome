module load libs/gcc/gsl/2.5 apps/gcc/R/4.3.1 libs/gcc/nlopt/2.7.1 compilers/gcc/11.2.0
nohup Rscript 1.4_large_sample_combine_pv.R C1 20 > 1.4_large_sample_combine_pv_C1_npctls20.out &
nohup Rscript 1.4_large_sample_combine_pv.R C2 20 > 1.4_large_sample_combine_pv_C2_npctls20.out &
nohup Rscript 1.4_large_sample_combine_pv.R C3 20 > 1.4_large_sample_combine_pv_C3_npctls20.out &