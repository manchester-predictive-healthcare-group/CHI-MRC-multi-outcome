/* Extract history of variables from primary care for AF, CHD, MI, Stroke, TIA */
/* This includes the time since condition was developed */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

/* For each I will use the %extract_medfile_event function, which extracts time until first record, and a censoring indicator */

*************************************************************************;
/* AF */
%extract_medfile_aurum(4);

%extract_medfile_hist(all_comor_codes, AF, A);
%extract_medfile_hist(all_comor_codes, AF, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

/* Check prevalence */
proc freq data = datint.varAF_primhist_A; tables AF_primhist ;run;
proc freq data = datint.varAF_primhist_B; tables AF_primhist ;run;

/* Check distribution of time since developing condition */
proc means data = datint.varAF_primhist_A; var AF_primhist_t; where AF_primhist = 1; run;
proc means data = datint.varAF_primhist_B; var AF_primhist_t; where AF_primhist = 1; run;


*************************************************************************;
/* HF */
%extract_medfile_aurum(22);

%extract_medfile_hist(all_comor_codes, HF, A);
%extract_medfile_hist(all_comor_codes, HF, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

/* Check prevalence */
proc freq data = datint.varHF_primhist_A; tables HF_primhist ;run;
proc freq data = datint.varHF_primhist_B; tables HF_primhist ;run;

/* Check distribution of time since developing condition */
proc means data = datint.varHF_primhist_A; var HF_primhist_t; where HF_primhist = 1; run;
proc means data = datint.varHF_primhist_B; var HF_primhist_t; where HF_primhist = 1; run;


*************************************************************************;
/* CHD */
%extract_medfile_aurum(14);

%extract_medfile_hist(all_comor_codes, CHD, A);
%extract_medfile_hist(all_comor_codes, CHD, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

/* Check prevalence */
proc freq data = datint.varCHD_primhist_A; tables CHD_primhist ;run;
proc freq data = datint.varCHD_primhist_B; tables CHD_primhist ;run;

/* Check distribution of time since developing condition */
proc means data = datint.varCHD_primhist_A; var CHD_primhist_t; where CHD_primhist = 1; run;
proc means data = datint.varCHD_primhist_B; var CHD_primhist_t; where CHD_primhist = 1; run;


*************************************************************************;
/* MI */
%extract_medfile_aurum(37);

%extract_medfile_hist(all_comor_codes, MI, A);
%extract_medfile_hist(all_comor_codes, MI, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

/* Check prevalence */
proc freq data = datint.varMI_primhist_A; tables MI_primhist ;run;
proc freq data = datint.varMI_primhist_B; tables MI_primhist ;run;

/* Check distribution of time since developing condition */
proc means data = datint.varMI_primhist_A; var MI_primhist_t; where MI_primhist = 1; run;
proc means data = datint.varMI_primhist_B; var MI_primhist_t; where MI_primhist = 1; run;


*************************************************************************;
/* Stroke */
%extract_medfile_aurum(34);

%extract_medfile_hist(all_comor_codes, Stroke, A);
%extract_medfile_hist(all_comor_codes, Stroke, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

/* Check prevalence */
proc freq data = datint.varStroke_primhist_A; tables Stroke_primhist ;run;
proc freq data = datint.varStroke_primhist_B; tables Stroke_primhist ;run;

/* Check distribution of time since developing condition */
proc means data = datint.varStroke_primhist_A; var Stroke_primhist_t; where Stroke_primhist = 1; run;
proc means data = datint.varStroke_primhist_B; var Stroke_primhist_t; where Stroke_primhist = 1; run;


*************************************************************************;
/* TIA */
%extract_medfile_aurum(35);

%extract_medfile_hist(all_comor_codes, TIA, A);
%extract_medfile_hist(all_comor_codes, TIA, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

/* Check prevalence */
proc freq data = datint.varTIA_primhist_A; tables TIA_primhist ;run;
proc freq data = datint.varTIA_primhist_B; tables TIA_primhist ;run;

/* Check distribution of time since developing condition */
proc means data = datint.varTIA_primhist_A; var TIA_primhist_t; where TIA_primhist = 1; run;
proc means data = datint.varTIA_primhist_B; var TIA_primhist_t; where TIA_primhist = 1; run;










