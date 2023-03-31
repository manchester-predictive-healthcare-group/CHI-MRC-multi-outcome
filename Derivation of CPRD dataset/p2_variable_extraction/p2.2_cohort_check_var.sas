/* Check distributions of each extracted variable */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


/* Here I will check for missingness, min/max, do a simple univariate model */
/* List of variables is here */


%let var_cont = Age BMI SBP Cholhdl_ratio;

%let var_disc = Smoking
		Ethnicity16
		Ethnicity6
		Alcohol_misuse
		Eating_disorders
		Asthma
		Atrial_fibrillation
		Anxiety_disorders
		Depression
		Bronchiectasis
		Hepatic_failure
		Chronic_viral_hepatitis
		Chronic_sinusitis
		COPD
		Coronary_heart_disease
		Dementia
		Diverticular_disease
		Epilepsy
		Hearing_loss
		Heart_failure
		Hypertension
		Irritable_bowel_syndrome
		Intellectual_disability
		Multiple_sclerosis
		Parkinsons_disease
		Peripheral_vascular_disease
		Psoriasis;

/*
		Visual_impairment_and_blindness
		Substance_misuse
		Rheumatoid_arthritis
		Schizophrenia
		Bipolar_affective_disorder_and_mania
		Stroke_not_otherwise_specified
		Transient_ischaemic_attack
		Thyroid_disease
		Myocardial_infarction
		Peptic_ulcer
		Inflammatory_bowel_disease
		Prostate
		Diabetes_type1
		Diabetes_type2
		Diabetes_vague
		CKD
*/

/* Macro to summarise continuous variables */
%macro summary_cont(gender,cohort);
title "continuous variables summary &gender &cohort";
proc means data=datint.cohort_var_&cohort. n nmiss mean std median min max range mode;var &var_cont;where gender = &gender;run;
%mend;

/* Macro to summarise discrete variables */
%macro summary_disc(gender,cohort);
title "discretes variables summary &gender &cohort";
proc freq data=datint.cohort_var_&cohort.;table &var_disc;where gender = &gender;run;
%mend;

/* Macro to fit model for continuous variables */
%macro univ_model_cont(var,gender,cohort);
%put &var;
proc phreg data=datint.cohort_var_&cohort.;
where gender = &gender;
model CVD_time*CVD_cens(1) = &var;
run;
%mend;

/* Macro to fit model for discrete variables */
%macro univ_model_disc(var,gender,cohort);
%put &var;
%if &var = 'Ethnicity16' %then %do;
proc phreg data=datint.cohort_var_&cohort.;
	where gender = &gender;
	class &var(ref='1');
	model CVD_time*CVD_cens(1) = &var;
run;
%end;
%if &var = 'Ethnicity6' %then %do;
proc phreg data=datint.cohort_var_&cohort.;
	where gender = &gender;
	class &var(ref='1');
	model CVD_time*CVD_cens(1) = &var;
run;
%end;
%if &var ^= 'Ethnicity' %then %do;
proc phreg data=datint.cohort_var_&cohort.;
	where gender = &gender;
	class &var(ref='0');
	model CVD_time*CVD_cens(1) = &var;
run;
%end;

%mend;


/* These cycle macros are only required if doing multiple analyses per variable */
/* (for example fitting a model and seeing if effect estimate is in correct direction) */
/* Currently I am not fitting any models */

/* Cycle through summary stuff for continuous variables */
%macro cycle_summary_cont(gender1,cohort1);
%local i next_name;
%let i=1;
%do %while (%scan(&var_cont, &i) ne );
   %let next_name = %scan(&var_cont, &i);
	%summary_cont(&next_name,&gender1,&cohort1);
   %let i = %eval(&i + 1);
%end;
%mend;


/* Cycle through summary stuff for discrete variables */
%macro cycle_summary_disc(gender1,cohort1);
%local i next_name;
%let i=1;
%do %while (%scan(&var_disc, &i) ne );
   %let next_name = %scan(&var_disc, &i);
	%summary_disc(&next_name,&gender1,&cohort1);
   %let i = %eval(&i + 1);
%end;
%mend;




/* Generate discerete variables summaries - cohort A */
ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_dist_disc_male_A.pdf";
%summary_disc(0,A);
ods pdf close;

ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_dist_disc_female_A.pdf";
%summary_disc(1,A);
ods pdf close;


/* Generate continuous variables summaries - cohort A */
ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_dist_cont_male_A.pdf";
%summary_cont(0,A);
ods pdf close;

ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_dist_cont_female_A.pdf";
%summary_cont(1,A);
ods pdf close;


/* Generate discerete variables summaries - cohort B */
ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_dist_disc_male_B.pdf";
%summary_disc(0,B);
ods pdf close;

ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_dist_disc_female_B.pdf";
%summary_disc(1,B);
ods pdf close;


/* Generate continuous variables summaries - cohort B */
ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_dist_cont_male_B.pdf";
%summary_cont(0,B);
ods pdf close;

ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_dist_cont_female_B.pdf";
%summary_cont(1,B);
ods pdf close;




