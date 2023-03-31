/* Calculate gfr scores and export to csv's for R algorithm which will calculate CKD from the gfr values */

*********************************************************************************************************************;
/************************************************** IMPORTANT *******************************************************/
/* AGE AND ETHNICITY VARIABLES MUST HAVE ALREADY BEEN DERIVED AND STORED IN APPROPRIATE PLACE FOR THIS CODE TO WORK */
/************************************************** IMPORTANT *******************************************************/
*********************************************************************************************************************;

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


******************************************************************************;
/* Now to extract gfr scores via test data */

/* There are three possible ways to extract using test data */
/* 1) gfr score recorded in database */
/* 2) egfr score recorded in database */
/* 3) Calculate egfr using creatinine measurements in database */

/* Extract the relevant codes for each */


*************************************************************************************;
/* Extract 1) gfr scores */

/* Extract all medcodes for gfr */
%extract_testfile_aurum(51);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of gfr, and the unit source vales */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;

/* put into named file */
/* Delete extreme values */
data gfr_dat;
	set all_comor_codes;
	if value <= 0 or value > 250 or value = . then delete;
run;

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


**************************************************************************************;
/* Extract 2) egfr scores */

/* Extract all medcodes for egfr */
%extract_testfile_aurum(52);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of egfr, and the unit source vales */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;

/* put into named file */
/* Delete extreme values */
data egfr_dat;
	set all_comor_codes;
	if value <= 0 or value > 250 or value = . then delete;
run;

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


***************************************************************************************;
/* Extract 3) creatinine scores, and convert to egfr scores */

/* Extract all medcodes for creatinine */
%extract_testfile_aurum(53);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of egfr, and the unit source vales */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;

/* put into named file */
/* Delete extreme values  */
data creatinine_dat;
	set all_comor_codes;
	if value <= 0 or value > 4000 or value = . then delete;
run;


/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

/* Check creatinine */
proc univariate data = creatinine_dat; var value;where unit_source_value = 285; run;
proc univariate data = creatinine_dat; var value;where unit_source_value = 218; run;

/* Only output entries where unit = umol/L, as other units seem to have different distribution, but the conversion does not make sense */
/* For example unit_source_value = 218 is mmol/L, but the observations are on average 1/10th of the size of unit_source_value = 285 (umol/L) (see the proc univariat procedures) */
/* Given the conversion from mmol/L to umol/L, you'd expect them to be about 1/1000th */
/* Observations with unit_source_value = 285 (umol/L) is approx 87% of measurements, so we are not missing much */
/* Also do a conversion from umol/L to mg/dl, required for egfr equation */
data creatinine_dat;
	set creatinine_dat;
	value = value/88.42;
	if unit_source_value = 285 then output;
run;




***************************************************************************************************************;
/* Now write a function that will: */
/* Look for gfr or egfr measurements for the specified time prior to the index date */
/* Look for creatinine measurements for the specified time prior to the index date, and convert to egfr using published formula */
/* These formula are dependents on age at index date, and ethnicity, so Age and Ethnicity variables must have already been derived */

/* Will write out in plain text first */
%macro derive_GFR(gfr_in, egfr_in, creatinine_in, cohort);

*****************************************;
/* Put gfr measurements into a dataset */

/* Merge with cohort */
data gfr_timeperiod;
	merge &gfr_in (in = ina) datint.cohort_base (in = inb);
	by person_id;
	if ina and inb then output;
run;

/* Remove any observations not prior to index date */
/* As this is a medical history variable, we do not require 5 years */
/* Also create a variable to say which source the gfr measurement came from */
data gfr_timeperiod (keep = person_id measurement_date gfr source source_num);
	set gfr_timeperiod (rename = (value = gfr));
	length source $30.;
	source = 'gfr';
	source_num = 0;
	if dtcens_combdeath > measurement_date and measurement_date > study_dtindex_&cohort then output;
run;

/* Check sorted by measurement date */
data _null_;
	set gfr_timeperiod;
	by person_id measurement_date;
run;


*****************************************;
/* Put egfr measurements into a dataset */

/* Merge with cohort */
data egfr_gp_timeperiod;
	merge &egfr_in (in = ina) datint.cohort_base (in = inb);
	by person_id;
	if ina and inb then output;
run;

/* Remove any observations not prior to index date */
/* As this is a medical history variable, we do not require 5 years */
/* Also create a variable to say which source the gfr measurement came from */
data egfr_gp_timeperiod (keep = person_id measurement_date gfr source source_num);
	set egfr_gp_timeperiod (rename = (value = gfr));
	length source $30.;
	source = 'egfr_gp';
	source_num = 1;
	if dtcens_combdeath > measurement_date and measurement_date > study_dtindex_&cohort then output;
run;

/* Check sorted by measurement date */
data _null_;
	set egfr_gp_timeperiod;
	by person_id measurement_date;
run;


********************************************************;
/* Convert creatinine to egfr, and put into a dataset */

/* Merge with cohort */
data creatinine_timeperiod;
	merge &creatinine_in (in = ina) datint.cohort_base (in = inb);
	by person_id;
	if ina and inb then output;
run;

/* Need to read in age, ethnicity variables and gender variables */
/* Remove any observations not prior to index date */
/* As this is a medical history variable, we do not require 5 years */
data creatinine_timeperiod (keep = person_id gender measurement_date creatinine);
	set creatinine_timeperiod;
	gender = gender_concept_id - 1;
	creatinine = value;
	if dtcens_combdeath > measurement_date and measurement_date > study_dtindex_&cohort then output;
run;
proc freq data = creatinine_timeperiod; tables gender;run;


/* Check sorted by measurement date */
data _null_;
	set creatinine_timeperiod;
	by person_id measurement_date;
run;

/* Now need to read in age and ethnicity values at the appropriate index date, to estimtae egfr from creatinine */
/* Age */
data age;
	set datint.varAge_&cohort.;
	proc sort; by person_id;
run;

/* Ethnicity */
data ethnicity;
	set datint.varEthnicity_&cohort.;
	proc sort; by person_id;
run;


/* Merge with the creatinine data */
data creatinine_a_e;
	merge creatinine_timeperiod (in = ina) age (in = inb) ethnicity (in = inc);
	if ina then output;
run;


/* Remove obesrvations with missing ethnicity, as we cannot estimate egfr without it */
data creatinine_a_e;
	set creatinine_a_e;
	if Ethnicity6 = . then delete;
run;

/* Remove obesrvations without stricly male or female gender, as we cannot estimate egfr without it */
data creatinine_a_e;
	set creatinine_a_e;
	if gender not in (0,1) then delete;
run;

/* Remove obesrvations with missing age, as we cannot estimate egfr without it */
data creatinine_a_e;
	set creatinine_a_e;
	if age = . then delete;
run;

proc contents data = creatinine_a_e; run;
proc print data = creatinine_a_e (obs= 50);run;

/* Calculate egfr based on creatinine levels */
data creatinine_a_e_gfr;
	set creatinine_a_e;
	if Ethnicity6 = 4 then do;
		if gender = 1 then do;
			if creatinine <= 0.7 then gfr = 166*(creatinine/0.7)**(-0.329)*0.993**Age;
			if creatinine > 0.7 then gfr = 166*(creatinine/0.7)**(-1.209)*0.993**Age;
		end;
		if gender = 0 then do;
			if creatinine <= 0.9 then gfr = 163*(creatinine/0.9)**(-0.411)*0.993**Age;
			if creatinine > 0.9 then gfr = 163*(creatinine/0.9)**(-1.209)*0.993**Age;
		end;
	end;
	if ethnicity ^= 4 then do;
		if gender = 1 then do;
			if creatinine <= 0.7 then gfr = 144*(creatinine/0.7)**(-0.329)*0.993**Age;
			if creatinine > 0.7 then gfr = 144*(creatinine/0.7)**(-1.209)*0.993**Age;
		end;
		if gender = 0 then do;
			if creatinine <= 0.9 then gfr = 141*(creatinine/0.9)**(-0.411)*0.993**Age;
			if creatinine > 0.9 then gfr = 141*(creatinine/0.9)**(-1.209)*0.993**Age;
		end;
	end;
run;

/* Create final datset, with egfr scores estimated from creatinine, with appropriate variable  names */
data egfr_est_timeperiod (keep = person_id measurement_date gfr source source_num);
	set creatinine_a_e_gfr;
	length source $30.;
	source = 'egfr_est';
	source_num = 2;
run;
proc print data = egfr_est_timeperiod (obs=10);run;


/* Concatenate these three datasets, remove any outside valid range, remove any with missing measurement_date, also add a numbered date */
data gfr_scores_&cohort.;
	set gfr_timeperiod egfr_gp_timeperiod egfr_est_timeperiod;
	format measurement_date_num best32.;
	measurement_date_num = measurement_date;
	if gfr <= 0 or gfr > 250 or gfr = . then delete;
run;



/* Sort by person_id, measurement date and source, and remove any duplicates */
/* By the ordering they have been given, for observations on the same day, preference is given to gfr scores, followed by egfr scores recorded by gp, followed by egfr scores estimated from creatinine in the database */
proc sort data = gfr_scores_&cohort. nodupkey; by person_id measurement_date source_num;run;
proc freq data = gfr_scores_&cohort.; tables source;run;

/* Compare distribution for each derivation type */
proc univariate data = gfr_scores_&cohort.; var gfr; where source = 'gfr'; run;
proc univariate data = gfr_scores_&cohort.; var gfr; where source = 'egfr_gp'; run;
proc univariate data = gfr_scores_&cohort.; var gfr; where source = 'egfr_est'; run;

/* Now order gfr scores by descending measurement date, ready for export to R where I have written algorithm to look derive CKD from these gfr scores */
proc sort data = gfr_scores_&cohort.; by person_id descending measurement_date;run;

/* Add a number, which reports which numer is patient is (will be used i nsubsequent algorithm) */
data gfr_scores_&cohort;
	set gfr_scores_&cohort;
	retain num;
	by person_id descending measurement_date;
	if first.person_id then num + 1;
	output;
run;

%mend;

%derive_GFR(gfr_dat, egfr_dat, creatinine_dat, A);
%derive_GFR(gfr_dat, egfr_dat, creatinine_dat, B);


/* The R code looks to see if ther eare two gfr measurements over 90 days apart which are both under 60, this is classed as CKD stage 3 */
proc export data=gfr_scores_A dbms = csv outfile="/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primevent_gfr_scores_A.csv" replace;run;
proc export data=gfr_scores_A dbms = tab outfile="/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primevent_gfr_scores_A.txt" replace;run;

proc export data=gfr_scores_B dbms = csv outfile="/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primevent_gfr_scores_B.csv" replace;run;
proc export data=gfr_scores_B dbms = tab outfile="/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primevent_gfr_scores_B.txt" replace;run;

