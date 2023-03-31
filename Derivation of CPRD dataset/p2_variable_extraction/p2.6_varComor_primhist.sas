/* Calculate all 'history of' variables, which just require existence of one medical code prior to the index date */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


/* Next Define the macro that will searched the codes after extracting them to create the variable */
%macro extract_medhist_var(varname,codes_in,cohort);

/* I now have the medcodes or prodcodes held in all_comor_codes */
/* First sort this file */
proc sort data = &codes_in; by person_id visit_start_date;run;

/* Merge this with the cohort file */
data all_comor_codes1;
	merge &codes_in (in=ina) datint.cohort_base (in=inb);
	by person_id;
	if ina and inb then output;
run;

/* Remove all events that happen after the index date */
data all_comor_codes2;
	set all_comor_codes1;
	if visit_start_date = . then delete;
	if visit_start_date > study_dtindex_&cohort then delete;
run;

/* Make it one line per patient */
data all_comor_codes3;
	set all_comor_codes2;
	by person_id visit_start_date;
	if first.person_id=1 then output;
run;
	
/* Remerge with the cohort, creating a 0/1 variable for the presence of the comorbidity */
data datint.var&varname._&cohort (keep = person_id &varname);
	merge datint.cohort_base (in=ina) all_comor_codes3 (in=inb);
	by person_id;
	if ina and inb then &varname = 1;
	if ina and ~inb then &varname = 0;
	if ina then output;
run;

%mend;


/*
proc delete data = all_comor_codes;run;
*/

/* For each variable: */
/* 1) Extract all codes for that variable */
/* 2) Generate history of variables prior to index date for each cohort */
/* 3) Delete the all_comor_codes file, before creating a new one */

***********************************************************************;
/* Alcohol misuse */
%extract_medfile_aurum(1);

%extract_medhist_var(Alcohol_misuse_primhist, all_comor_codes, A);
%extract_medhist_var(Alcohol_misuse_primhist, all_comor_codes, B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Alcohol finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Eating_disorders,  */
%extract_medfile_aurum(2);

%extract_medhist_var(Eating_disorders_primhist, all_comor_codes,A);
%extract_medhist_var(Eating_disorders_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Eating_disorders,   finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Asthma,  */
%extract_medfile_aurum(3);

%extract_medhist_var(Asthma_primhist, all_comor_codes,A);
%extract_medhist_var(Asthma_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Asthma,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Anxiety_disorders,  */
%extract_medfile_aurum(5);

%extract_medhist_var(Anxiety_disorders_primhist, all_comor_codes,A);
%extract_medhist_var(Anxiety_disorders_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Anxiety_disorders,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Depression,  */
%extract_medfile_aurum(6);

%extract_medhist_var(Depression_primhist, all_comor_codes,A);
%extract_medhist_var(Depression_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Depression,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Visual_impairment_and_blindness,  */
%extract_medfile_aurum(7);

%extract_medhist_var(Visual_impairment_primhist, all_comor_codes,A);
%extract_medhist_var(Visual_impairment_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Visual_impairment_and_blindness,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Bronchiectasis,  */
%extract_medfile_aurum(8);

%extract_medhist_var(Bronchiectasis_primhist, all_comor_codes,A);
%extract_medhist_var(Bronchiectasis_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Bronchiectasis,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Hepatic_failure,  */
%extract_medfile_aurum(10);

%extract_medhist_var(Hepatic_failure_primhist, all_comor_codes,A);
%extract_medhist_var(Hepatic_failure_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Hepatic_failure,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Chronic_viral_hepatitis,  */
%extract_medfile_aurum(11);

%extract_medhist_var(Viral_hepatitis_primhist, all_comor_codes,A);
%extract_medhist_var(Viral_hepatitis_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Chronic_viral_hepatitis,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Chronic_sinusitis,  */
%extract_medfile_aurum(12);

%extract_medhist_var(Sinusitis_primhist, all_comor_codes,A);
%extract_medhist_var(Sinusitis_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Chronic_sinusitis,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* COPD,  */
%extract_medfile_aurum(13);

%extract_medhist_var(COPD_primhist, all_comor_codes,A);
%extract_medhist_var(COPD_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put COPD,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Dementia,  */
%extract_medfile_aurum(15);

%extract_medhist_var(Dementia_primhist, all_comor_codes,A);
%extract_medhist_var(Dementia_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Dementia,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Diverticular_disease,  */
%extract_medfile_aurum(19);

%extract_medhist_var(Diverticular_primhist, all_comor_codes,A);
%extract_medhist_var(Diverticular_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Diverticular_disease,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Epilepsy,  */
%extract_medfile_aurum(20);

%extract_medhist_var(Epilepsy_primhist, all_comor_codes,A);
%extract_medhist_var(Epilepsy_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Epilepsy,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Hearing_loss,  */
%extract_medfile_aurum(21);

%extract_medhist_var(Hearing_loss_primhist, all_comor_codes,A);
%extract_medhist_var(Hearing_loss_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Hearing_loss,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Hypertension,  */
%extract_medfile_aurum(23);

%extract_medhist_var(Hypertension_primhist, all_comor_codes,A);
%extract_medhist_var(Hypertension_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Hypertension,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Irritable_bowel_syndrome,  */
%extract_medfile_aurum(24);

%extract_medhist_var(IBS_primhist, all_comor_codes,A);
%extract_medhist_var(IBS_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Irritable_bowel_syndrome,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Intellectual_disability,  */
%extract_medfile_aurum(25);

%extract_medhist_var(Intellectual_dis_primhist, all_comor_codes,A);
%extract_medhist_var(Intellectual_dis_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Intellectual_disability,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Multiple_sclerosis,  */
%extract_medfile_aurum(26);

%extract_medhist_var(MS_primhist, all_comor_codes,A);
%extract_medhist_var(MS_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Multiple_sclerosis,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Parkinsons_disease,  */
%extract_medfile_aurum(27);

%extract_medhist_var(Parkinsons_primhist, all_comor_codes,A);
%extract_medhist_var(Parkinsons_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Parkinsons_disease,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Peripheral_vascular_disease,  */
%extract_medfile_aurum(28);

%extract_medhist_var(Perip_vascular_primhist, all_comor_codes,A);
%extract_medhist_var(Perip_vascular_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Peripheral_vascular_disease,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Psoriasis,  */
%extract_medfile_aurum(29);

%extract_medhist_var(Psoriasis_primhist, all_comor_codes,A);
%extract_medhist_var(Psoriasis_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Psoriasis,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Substance_misuse,  */
%extract_medfile_aurum(30);

%extract_medhist_var(Substance_misuse_primhist, all_comor_codes,A);
%extract_medhist_var(Substance_misuse_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Substance_misuse,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Rheumatoid_arthritis,  */
%extract_medfile_aurum(31);

%extract_medhist_var(RA_primhist, all_comor_codes,A);
%extract_medhist_var(RA_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Rheumatoid_arthritis,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Schizophrenia,  */
%extract_medfile_aurum(32);

%extract_medhist_var(Schizophrenia_primhist, all_comor_codes,A);
%extract_medhist_var(Schizophrenia_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Schizophrenia,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Bipolar_affective_disorder_and_mania,  */
%extract_medfile_aurum(33);

%extract_medhist_var(Bipolar_primhist, all_comor_codes,A);
%extract_medhist_var(Bipolar_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Bipolar_affective_disorder_and_mania,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Thyroid_disease,  */
%extract_medfile_aurum(36);

%extract_medhist_var(Thyroid_primhist, all_comor_codes,A);
%extract_medhist_var(Thyroid_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Thyroid_disease,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Peptic_ulcer */
%extract_medfile_aurum(38);

%extract_medhist_var(Peptic_ulcer_primhist, all_comor_codes,A);
%extract_medhist_var(Peptic_ulcer_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Peptic_ulcer finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Inflammatory_bowel_disease */
%extract_medfile_aurum(39);

%extract_medhist_var(IBD_primhist, all_comor_codes,A);
%extract_medhist_var(IBD_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Inflammatory_bowel_disease finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


***********************************************************************;
/* Prostate */
%extract_medfile_aurum(40);

%extract_medhist_var(Prostate_primhist, all_comor_codes,A);
%extract_medhist_var(Prostate_primhist, all_comor_codes,B);

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%put Prostate finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;

