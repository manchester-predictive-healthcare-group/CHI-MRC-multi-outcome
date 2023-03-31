cd "/mnt/bmh01-rds/mrc-multi-outcome/Project_6/code/"
module load apps/gcc/R/4.0.2-gcc830
nohup Rscript p3.1.2_large_sample_analysis_pv_prep_pctls.R M1C1 200000 1 > p3.1.2_large_sample_analysis_pv_prep_pctls_200000_M1C1_npctls1.out &
nohup Rscript p3.1.2_large_sample_analysis_pv_prep_pctls.R M1C2 200000 1 > p3.1.2_large_sample_analysis_pv_prep_pctls_200000_M1C2_npctls1.out &
nohup Rscript p3.1.2_large_sample_analysis_pv_prep_pctls.R M1C3 200000 1 > p3.1.2_large_sample_analysis_pv_prep_pctls_200000_M1C3_npctls1.out &
nohup Rscript p3.1.2_large_sample_analysis_pv_prep_pctls.R M1C1 200000 10 > p3.1.2_large_sample_analysis_pv_prep_pctls_200000_M1C1_npctls10.out &
nohup Rscript p3.1.2_large_sample_analysis_pv_prep_pctls.R M1C2 200000 10 > p3.1.2_large_sample_analysis_pv_prep_pctls_200000_M1C2_npctls10.out &
nohup Rscript p3.1.2_large_sample_analysis_pv_prep_pctls.R M1C3 200000 10 > p3.1.2_large_sample_analysis_pv_prep_pctls_200000_M1C3_npctls10.out &
nohup Rscript p3.1.2_large_sample_analysis_pv_prep_pctls.R M1C1 200000 20 > p3.1.2_large_sample_analysis_pv_prep_pctls_200000_M1C1_npctls20.out &
nohup Rscript p3.1.2_large_sample_analysis_pv_prep_pctls.R M1C2 200000 20 > p3.1.2_large_sample_analysis_pv_prep_pctls_200000_M1C2_npctls20.out &
nohup Rscript p3.1.2_large_sample_analysis_pv_prep_pctls.R M1C3 200000 20 > p3.1.2_large_sample_analysis_pv_prep_pctls_200000_M1C3_npctls20.out &