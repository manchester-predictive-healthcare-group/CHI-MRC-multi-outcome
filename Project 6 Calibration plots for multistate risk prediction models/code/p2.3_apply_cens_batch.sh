cd "/mnt/bmh01-rds/mrc-multi-outcome/Project_6/code/"
module load apps/gcc/R/4.0.2-gcc830
nohup Rscript p2.3_apply_cens.R M1 C1 > p2.3_apply_cens_M1C1.out &
nohup Rscript p2.3_apply_cens.R M1 C2 > p2.3_apply_cens_M1C2.out &
nohup Rscript p2.3_apply_cens.R M1 C3 > p2.3_apply_cens_M1C3.out &