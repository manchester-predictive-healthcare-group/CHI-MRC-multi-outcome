/* Extract outcome variables from primary care for AF, CHD, MI, Stroke, TIA */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

/* For each I will use the %extract_medfile_event function, which extracts time until first record, and a censoring indicator */

*************************************************************************;
/* AF */
%extract_medfile_aurum(4);

%extract_medfile_event(all_comor_codes, AF, A);
%extract_medfile_event(all_comor_codes, AF, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

proc means data = datint.varAF_primev_A; var AF_primev_t; run;
proc means data = datint.varAF_primev_B; var AF_primev_t; run;

proc freq data = datint.varAF_primev_A; tables AF_primev_c ;run;
proc freq data = datint.varAF_primev_B; tables AF_primev_c ;run;



*************************************************************************;
/* HF */
%extract_medfile_aurum(22);

%extract_medfile_event(all_comor_codes, HF, A);
%extract_medfile_event(all_comor_codes, HF, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

proc means data = datint.varHF_primev_A; var HF_primev_t; run;
proc means data = datint.varHF_primev_B; var HF_primev_t; run;

proc freq data = datint.varHF_primev_A; tables HF_primev_c ;run;
proc freq data = datint.varHF_primev_B; tables HF_primev_c ;run;


*************************************************************************;
/* CHD */
%extract_medfile_aurum(14);

%extract_medfile_event(all_comor_codes, CHD, A);
%extract_medfile_event(all_comor_codes, CHD, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

proc means data = datint.varCHD_primev_A; var CHD_primev_t; run;
proc means data = datint.varCHD_primev_B; var CHD_primev_t; run;

proc freq data = datint.varCHD_primev_A; tables CHD_primev_c ;run;
proc freq data = datint.varCHD_primev_B; tables CHD_primev_c ;run;


*************************************************************************;
/* MI */
%extract_medfile_aurum(37);

%extract_medfile_event(all_comor_codes, MI, A);
%extract_medfile_event(all_comor_codes, MI, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

proc means data = datint.varMI_primev_A; var MI_primev_t; run;
proc means data = datint.varMI_primev_B; var MI_primev_t; run;

proc freq data = datint.varMI_primev_A; tables MI_primev_c ;run;
proc freq data = datint.varMI_primev_B; tables MI_primev_c ;run;


*************************************************************************;
/* Stroke */
%extract_medfile_aurum(34);

%extract_medfile_event(all_comor_codes, Stroke, A);
%extract_medfile_event(all_comor_codes, Stroke, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

proc means data = datint.varStroke_primev_A; var Stroke_primev_t; run;
proc means data = datint.varStroke_primev_B; var Stroke_primev_t; run;

proc freq data = datint.varStroke_primev_A; tables Stroke_primev_c ;run;
proc freq data = datint.varStroke_primev_B; tables Stroke_primev_c ;run;


*************************************************************************;
/* TIA */
%extract_medfile_aurum(35);

%extract_medfile_event(all_comor_codes, TIA, A);
%extract_medfile_event(all_comor_codes, TIA, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

proc means data = datint.varTIA_primev_A; var TIA_primev_t; run;
proc means data = datint.varTIA_primev_B; var TIA_primev_t; run;

proc freq data = datint.varTIA_primev_A; tables TIA_primev_c ;run;
proc freq data = datint.varTIA_primev_B; tables TIA_primev_c ;run;












