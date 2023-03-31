/* Calculate Smoking */

options ps=60 ls=70;

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

****************************************************************************************************;


/* First extract test data for Smoking */

/* Extract all medcodes for Smoking from both medical and test file */

******************************************************************;
/* Test file */
%extract_testfile_aurum(46);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of Smoking, and the unit source values */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;
proc print data = all_comor_codes(obs=50);run;

/* Put into named file */
/* Rename variables to match that from medfile, so we can concatenate */
data smoking_dat_test;
	set all_comor_codes (rename = (measurement_date = visit_start_date measurement_concept_id = visit_concept_id));
	if visit_start_date = . then delete;
	keep person_id visit_start_date visit_concept_id;
run;

proc print data = smoking_dat_test (obs=50);run;

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


*******************************************************************;
/* Medical file */
%extract_medfile_aurum(46);

/* sort the data */
proc sort data = all_comor_codes;by person_id visit_start_date; run;

/* Put into named file */
/* Delete all observations outside of valid range */
data smoking_dat_med;
	set all_comor_codes;
	if visit_start_date = . then delete;
	keep person_id visit_start_date visit_concept_id;
run;


/* Combine medical and test files into one */
data smoking_comb;
	set smoking_dat_test smoking_dat_med;
run;	



*************************************************************************************;
/* Merge Smoking_comb with the code list, which contains info on non-smoker, ex-smoker, ex or current, or current */

/* Import code list */
%importcsv(smoking_codes ,&ffile_code.cr_smokingcodes_aurum_mcid, ddatarow = 2);

/* Change medcodeid to be visit_concept_id, create a numeric version of smokstatus, and retain only relevant variables */
data smoking_codes;
	set smoking_codes;
	visit_concept_id = input(medcodeid, best32.);
	format visit_concept_id best32.;
	if smokstatus = 'non-smoker' then smokstatus_num = 0;
	if smokstatus = 'ex-smoker' then smokstatus_num = 1;
	if smokstatus = 'current smoker' then smokstatus_num = 2;
	if smokstatus = 'current or ex-smoker' then smokstatus_num = 3;
	keep visit_concept_id smokstatus smokstatus_num;
run;


/* Sort both datasets by visit_concept_id */
proc sort data = smoking_comb;by visit_concept_id;run;
proc sort data = smoking_codes;by visit_concept_id;run;

/* Merge extracted smoking codes with the code list */
data smoking_comb_codes;
	merge smoking_comb (in = ina) smoking_codes (in = inb);
	by visit_concept_id;
	if ina then output;
run;


/* Sort this dataset by person_d, and visit_start_id */
proc sort data = smoking_comb_codes; by person_id visit_start_date;run;
proc print data = smoking_comb_codes (obs = 100);run;


/* Now create a new smoking status variable based on the following rules */
/* No = non-smoker, Ex = Ex-smoker, Cx = Current or Ex-smoker, C = current smoker */
/* If No after Ex, Cx or C, change No to Ex */
data smoking_comb_codes2;
	set smoking_comb_codes;
	retain retain_smokstatus_num;
	by person_id visit_start_date;
	if first.person_id then retain_smokstatus_num = smokstatus_num;
	else do;
		if smokstatus_num > 0 then do;
			retain_smokstatus_num = smokstatus_num;
		end;
	end;
	output;
run;

proc print data = smoking_comb_codes2(obs=500);var person_id visit_start_date smokstatus_num retain_smokstatus_num;run;


****************************************************************************************;
/* The variable retain_smokstatus_num is our variable of interest */

/* Now merge with the base cohort, and apply the following rules: */

/* Look for most recent No, Ex or C code in 5 years prior to index date */
/* If there are no codes in last 5 years, set smoking to missing */

/* Create another variable which = 1 if a patient has any history of smoking (Ex, Cx or C), not within 5 year follow up */
/* As we know these patients cannot be a "Non smoker" */


/* Write a macro to do this */
%macro create_smoking(smoking_in, cohort_in, timemeas1_in);

/* First create the variable which is precise (i.e. is there a No, Ex or C record in the last 5 years */

/* Start by merging cohort with smoking statuses */
data smokdat;
	merge datint.cohort_base (in = ina) &smoking_in (in = inb);
	by person_id;
	if ina and inb then output;
run;

/* Remove all observations where retain_smokstatus_num = 3 */
data smokdat2;
	set smokdat;
	if retain_smokstatus_num = 3 then delete;
run;

/* Now want to only output measurements within timemeas1 and timemeas2 of the index date in either direction */
data smokdat3;
	set smokdat2;
	if visit_start_date <= study_dtindex_&cohort_in and visit_start_date > study_dtindex_&cohort_in - &timemeas1_in then output;
run;

/* Now only want to output the most recent observation for each individual */
data smokdat4;
	set smokdat3 (rename = (retain_smokstatus_num = Smoking));
	by person_id visit_start_date;
	if last.person_id then output;
	keep person_id Smoking;
run;


/* Now create variable which = 1 if there is a history of smoking at some point, and 0 otherwise */
/* We can therefore distinguish people with missing status between those who had some history of smoking status */
/* We know these individuals must be Ex or C, as opposed to individuals with non-smokin history, or no recorded codes, could be No, Ex or C */
data smokdat_hist;
	merge datint.cohort_base (in = ina) &smoking_in (in = inb);
	by person_id;
	if ina and inb then output;
run;

/* Delete any observations which = 0 */
data smokdat_hist2;
	set smokdat_hist;
	if retain_smokstatus_num = 0 then delete;
run;

/* Also only retain observations that have happened prior to the index date */
data smokdat_hist3;
	set smokdat_hist2;
	if visit_start_date <= study_dtindex_&cohort_in then output;
run;

/* Reduce to 1 observation per person */
data smokdat_hist4;
	set smokdat_hist3;
	by person_id visit_start_date;
	if last.person_id then output;
	keep person_id;
run;

/* We know all these remaining individuals have a history of smoking at some point prior to the index date, so we only retain person_id */


/* Combine both these datasets with the base cohort to create variables of interest */
data datint.varSmoking_&cohort_in (keep = person_id Smoking Smoking_anyhist);
	merge datint.cohort_base (in=ina) smokdat4 (in=inb) smokdat_hist4 (in = inc);
	by person_id;
	if ina and not inb then Smoking = . ;
	if inc then Smoking_anyhist = 1;
	if not inc then Smoking_anyhist = 0;
	if ina then output;
run;

%mend;

/* Apply function */
%create_smoking(smoking_comb_codes2, A, 1826.25);
%create_smoking(smoking_comb_codes2, B, 1826.25);


/* Check and summarise data */
proc freq data = datint.varSmoking_A; tables Smoking Smoking_anyhist;run;
proc freq data = datint.varSmoking_B; tables Smoking Smoking_anyhist;run;

proc means data = datint.varSmoking_A n nmiss; var Smoking Smoking_anyhist;run;
proc means data = datint.varSmoking_B n nmiss; var Smoking Smoking_anyhist;run;


