/* Calculate BP */

options ps=60 ls=70;

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

****************************************************************************************************;

proc print data = datext.alltestfmt_1(obs=200 where=(unit_source_value = 216)); var measurement_concept_id unit_source_value value;run;
proc print data = datext.alltestfmt_1(obs=200 where=(unit_source_value = 210)); var measurement_concept_id unit_source_value value;run;
proc print data = datext.alltestfmt_1(obs=200 where=(unit_source_value = 455)); var measurement_concept_id unit_source_value value;run;

/* First extract test data for BP */

/* Extract all medcodes for BP */
%extract_testfile_aurum(50);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of SBP, and the unit source values */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;

/* Put into named file */
/* Delete all observations outside of valid range */
data bp_dat;
	set all_comor_codes;
	if value < 60 or value > 230 or value = . then delete;
run;


/* Write the function that will look for most recent BP value within specified time period */

%macro create_bp(bp_in, cohort_in, timemeas1_in);

/* Start by using function from var_extract1.sas macro container, which will extract most recent test result */
%extract_latest_test_timeperiod(&bp_in, bp_cohort, &timemeas1_in, 0, &cohort_in, BP);


/* Merge this with base cohort */
data datint.varBP_&cohort_in (keep = person_id BP);
	merge datint.cohort_base (in=ina) bp_cohort (in=inb);
	by person_id;
	if ina then output;
	if ina and not inb then BP = .;
run;

%mend;

/* Apply function */
%create_bp(bp_dat, A, 1826.25);
%create_bp(bp_dat, B, 1826.25);


/* Check and summarise data */
proc univariate data = datint.varBP_A; var BP;run;
proc univariate data = datint.varBP_B; var BP;run;

proc means data = datint.varBP_A n nmiss; var BP;run;
proc means data = datint.varBP_B n nmiss; var BP;run;

