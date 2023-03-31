/*************************************************************************/
/* Extract Atrial fibrillation history and event variables from HES data */
/*************************************************************************/

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

*************************************************************************;
/* Extract history of */
%extract_hes_hist(icd_primary3, icd3af., AF, A);
%extract_hes_hist(icd_primary3, icd3af., AF, B);

/* Check prevalence */
proc freq data = datint.varAF_heshist_A; tables AF_heshist;run;
proc freq data = datint.varAF_heshist_B; tables AF_heshist;run;

/* Check distribution of time since developing condition */
proc means data = datint.varAF_heshist_A; var AF_heshist_t; where AF_heshist = 1; run;
proc means data = datint.varAF_heshist_B; var AF_heshist_t; where AF_heshist = 1; run;


*************************************************************************;
/* Extract event, time until and censoring indicator */
%extract_hes_event(icd_primary3, icd3af., AF, A);
%extract_hes_event(icd_primary3, icd3af., AF, B);

/* Check distribution and prevalence of censoring indicator */
proc means data = datint.varAF_hesev_A; var AF_hesev_t;run;
proc means data = datint.varAF_hesev_B; var AF_hesev_t;run;

proc freq data = datint.varAF_hesev_A; tables AF_hesev_c;run;
proc freq data = datint.varAF_hesev_B; tables AF_hesev_c;run;




/*************************************************************************/
/* Extract CHD and MI history and event variables from HES data */
/*************************************************************************/

*************************************************************************;
/* Extract history of */
%extract_hes_hist(icd_primary3, icd3chd., CHD_MI, A);
%extract_hes_hist(icd_primary3, icd3chd., CHD_MI, B);

/* Check prevalence */
proc freq data = datint.varCHD_MI_heshist_A; tables CHD_MI_heshist;run;
proc freq data = datint.varCHD_MI_heshist_B; tables CHD_MI_heshist;run;

/* Check distribution of time since developing condition */
proc means data = datint.varCHD_MI_heshist_A; var CHD_MI_heshist_t; where CHD_MI_heshist = 1; run;
proc means data = datint.varCHD_MI_heshist_B; var CHD_MI_heshist_t; where CHD_MI_heshist = 1; run;


*************************************************************************;
/* Extract event, time until and censoring indicator */
%extract_hes_event(icd_primary3, icd3chd., CHD_MI, A);
%extract_hes_event(icd_primary3, icd3chd., CHD_MI, B);

/* Check distribution and prevalence of censoring indicator */
proc means data = datint.varCHD_MI_hesev_A; var CHD_MI_hesev_t;run;
proc means data = datint.varCHD_MI_hesev_B; var CHD_MI_hesev_t;run;

proc freq data = datint.varCHD_MI_hesev_A; tables CHD_MI_hesev_c;run;
proc freq data = datint.varCHD_MI_hesev_B; tables CHD_MI_hesev_c;run;



/*************************************************************************/
/* Extract CKD history and event variables from HES data */
/*************************************************************************/

*************************************************************************;
/* Extract history of */
%extract_hes_hist(icd5, icd5ckd., CKD, A);
%extract_hes_hist(icd5, icd5ckd., CKD, B);

/* Check prevalence */
proc freq data = datint.varCKD_heshist_A; tables CKD_heshist;run;
proc freq data = datint.varCKD_heshist_B; tables CKD_heshist;run;

/* Check distribution of time since developing condition */
proc means data = datint.varCKD_heshist_A; var CKD_heshist_t; where CKD_heshist = 1; run;
proc means data = datint.varCKD_heshist_B; var CKD_heshist_t; where CKD_heshist = 1; run;


*************************************************************************;
/* Extract event, time until and censoring indicator */
%extract_hes_event(icd5, icd5ckd., CKD, A);
%extract_hes_event(icd5, icd5ckd., CKD, B);

/* Check distribution and prevalence of censoring indicator */
proc means data = datint.varCKD_hesev_A; var CKD_hesev_t;run;
proc means data = datint.varCKD_hesev_B; var CKD_hesev_t;run;

proc freq data = datint.varCKD_hesev_A; tables CKD_hesev_c;run;
proc freq data = datint.varCKD_hesev_B; tables CKD_hesev_c;run;



/*************************************************************************/
/* Extract T2D history and event variables from HES data */
/*************************************************************************/

*************************************************************************;
/* Extract history of */
%extract_hes_hist(icd_primary3, icd3t2d., Diab_t2, A);
%extract_hes_hist(icd_primary3, icd3t2d., Diab_t2, B);

/* Check prevalence */
proc freq data = datint.varDiab_t2_heshist_A; tables Diab_t2_heshist;run;
proc freq data = datint.varDiab_t2_heshist_B; tables Diab_t2_heshist;run;

/* Check distribution of time since developing condition */
proc means data = datint.varDiab_t2_heshist_A; var Diab_t2_heshist_t; where Diab_t2_heshist = 1; run;
proc means data = datint.varDiab_t2_heshist_B; var Diab_t2_heshist_t; where Diab_t2_heshist = 1; run;



*************************************************************************;
/* Extract event, time until and censoring indicator */
%extract_hes_event(icd_primary3, icd3t2d., Diab_t2, A);
%extract_hes_event(icd_primary3, icd3t2d., Diab_t2, B);

/* Check distribution and prevalence of censoring indicator */
proc means data = datint.varDiab_t2_hesev_A; var Diab_t2_hesev_t;run;
proc means data = datint.varDiab_t2_hesev_B; var Diab_t2_hesev_t;run;

proc freq data = datint.varDiab_t2_hesev_A; tables Diab_t2_hesev_c;run;
proc freq data = datint.varDiab_t2_hesev_B; tables Diab_t2_hesev_c;run;



/*************************************************************************/
/* Extract T1D history and event variables from HES data */
/*************************************************************************/

*************************************************************************;
/* Extract history of */
%extract_hes_hist(icd_primary3, icd3t1d., Diab_t1, A);
%extract_hes_hist(icd_primary3, icd3t1d., Diab_t1, B);

/* Check prevalence */
proc freq data = datint.varDiab_t1_heshist_A; tables Diab_t1_heshist;run;
proc freq data = datint.varDiab_t1_heshist_B; tables Diab_t1_heshist;run;

/* Check distribution of time since developing condition */
proc means data = datint.varDiab_t1_heshist_A; var Diab_t1_heshist_t; where Diab_t1_heshist = 1; run;
proc means data = datint.varDiab_t1_heshist_B; var Diab_t1_heshist_t; where Diab_t1_heshist = 1; run;


*************************************************************************;
/* Extract event, time until and censoring indicator */
%extract_hes_event(icd_primary3, icd3t1d., Diab_t1, A);
%extract_hes_event(icd_primary3, icd3t1d., Diab_t1, B);

/* Check distribution and prevalence of censoring indicator */
proc means data = datint.varDiab_t1_hesev_A; var Diab_t1_hesev_t;run;
proc means data = datint.varDiab_t1_hesev_B; var Diab_t1_hesev_t;run;

proc freq data = datint.varDiab_t1_hesev_A; tables Diab_t1_hesev_c;run;
proc freq data = datint.varDiab_t1_hesev_B; tables Diab_t1_hesev_c;run;



/*************************************************************************/
/* Extract Stroke/TIA history and event variables from HES data */
/*************************************************************************/

*************************************************************************;
/* Extract history of */
%extract_hes_hist(icd_primary3, icd3str., Stroke_TIA, A);
%extract_hes_hist(icd_primary3, icd3str., Stroke_TIA, B);

/* Check prevalence */
proc freq data = datint.varStroke_TIA_heshist_A; tables Stroke_TIA_heshist;run;
proc freq data = datint.varStroke_TIA_heshist_B; tables Stroke_TIA_heshist;run;

/* Check distribution of time since developing condition */
proc means data = datint.varStroke_TIA_heshist_A; var Stroke_TIA_heshist_t; where Stroke_TIA_heshist = 1; run;
proc means data = datint.varStroke_TIA_heshist_B; var Stroke_TIA_heshist_t; where Stroke_TIA_heshist = 1; run;


*************************************************************************;
/* Extract event, time until and censoring indicator */
%extract_hes_event(icd_primary3, icd3str., Stroke_TIA, A);
%extract_hes_event(icd_primary3, icd3str., Stroke_TIA, B);

/* Check distribution and prevalence of censoring indicator */
proc means data = datint.varStroke_TIA_hesev_A; var Stroke_TIA_hesev_t;run;
proc means data = datint.varStroke_TIA_hesev_B; var Stroke_TIA_hesev_t;run;

proc freq data = datint.varStroke_TIA_hesev_A; tables Stroke_TIA_hesev_c; run;
proc freq data = datint.varStroke_TIA_hesev_B; tables Stroke_TIA_hesev_c;run;



/*************************************************************************/
/* Extract heart failure history and event variables from HES data */
/*************************************************************************/

*************************************************************************;
/* Extract history of */
%extract_hes_hist(icd_primary3, icd3hf., HF, A);
%extract_hes_hist(icd_primary3, icd3hf., HF, B);

/* Check prevalence */
proc freq data = datint.varHF_heshist_A; tables HF_heshist;run;
proc freq data = datint.varHF_heshist_B; tables HF_heshist;run;

/* Check distribution of time since developing condition */
proc means data = datint.varHF_heshist_A; var HF_heshist_t; where HF_heshist = 1; run;
proc means data = datint.varHF_heshist_B; var HF_heshist_t; where HF_heshist = 1; run;


*************************************************************************;
/* Extract event, time until and censoring indicator */
%extract_hes_event(icd_primary3, icd3hf., HF, A);
%extract_hes_event(icd_primary3, icd3hf., HF, B);

/* Check distribution and prevalence of censoring indicator */
proc means data = datint.varHF_hesev_A; var HF_hesev_t;run;
proc means data = datint.varHF_hesev_B; var HF_hesev_t;run;

proc freq data = datint.varHF_hesev_A; tables HF_hesev_c; run;
proc freq data = datint.varHF_hesev_B; tables HF_hesev_c;run;