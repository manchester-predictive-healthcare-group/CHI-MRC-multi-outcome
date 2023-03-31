****************************************************************************;
/* This progrram will export the analysis datasets to R */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

/* R uses a different origin for time = 0 */
/* While it can be changed, i'm going to create variables which will automatically be correct for the default time origin in R */

/* Write a macro that will create a temporary dataset ready for export */
%macro create_temp(cohort);
data temp&cohort.;
	set datint.cohort_var_&cohort.;
	format study_dtindex_&cohort dtcens_r dtcens_combdeath_r death_date_r death_date_ons_r death_date_comb_r 8.2;
	study_dtindex_&cohort._r = study_dtindex_&cohort. - '01jan1970'd;
	dtcens_r = dtcens - '01jan1970'd;
	dtcens_combdeath_r = dtcens_combdeath - '01jan1970'd;
	death_date_r = death_date - '01jan1970'd;
	death_date_ons_r = death_date_ons - '01jan1970'd;
	death_date_comb_r = death_date_comb - '01jan1970'd;

run;
%mend;

%create_temp(A);
%create_temp(B);


/* Now output the datasets into csv files */

proc export data= tempA dbms = csv outfile="/mnt/bmh01-rds/mrc-multi-outcome/data_Aurum_65plus/data_intermediate/sas_cohort_var_A.csv" replace;run;
proc export data= tempB dbms = csv outfile="/mnt/bmh01-rds/mrc-multi-outcome/data_Aurum_65plus/data_intermediate/sas_cohort_var_B.csv" replace;run;


