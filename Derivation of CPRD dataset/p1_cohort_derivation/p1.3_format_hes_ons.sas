/* Format and sort ICD and HES data */

*/ file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

proc sort data = dathes.hes_diagnosis_epi; by person_id epistart; run;

proc sort data = dathes.hes_diagnosis_hosp; by person_id admidate; run;

proc sort data = dathes.hes_diagnosis_primary; by person_id admidate; run;

proc sort data = dathes.hes_episodes;by person_id discharged; run;

proc sort data = dathes.hes_hospital; by person_id admidate; run;

proc sort data = dathes.hes_patient; by person_id; run;

proc sort data = dathes.imd_patient; by person_id; run;

proc sort data = dathes.imd_practice; by care_site_id; run;

proc sort data = dathes.onsdeath; by person_id death_date_ons; run;
