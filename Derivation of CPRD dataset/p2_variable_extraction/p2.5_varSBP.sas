/* Calculate SBP */

options ps=60 ls=70;

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

****************************************************************************************************;


/* First extract test data for SBP */

/* Extract all medcodes for SBP */
%extract_testfile_aurum(45);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of SBP, and the unit source values */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;
proc print data = all_comor_codes(obs=200);run;

/* Put into named file */
/* Delete all observations outside of valid range */
data sbp_dat;
	set all_comor_codes;
	if value < 60 or value > 230 or value = . then delete;
run;

proc print data = all_comor_codes(obs=50); where measurement_concept_id = 3676271000006117;run;


/* Write the function that will look for most recent SBP value within specified time period */

%macro create_sbp(sbp_in, cohort_in, timemeas1_in);

/* Start by using function from var_extract1.sas macro container, which will extract most recent test result */
%extract_latest_test_timeperiod(&sbp_in, sbp_cohort, &timemeas1_in, 0, &cohort_in, SBP);


/* Merge this with base cohort */
data datint.varSBP_&cohort_in (keep = person_id SBP);
	merge datint.cohort_base (in=ina) sbp_cohort (in=inb);
	by person_id;
	if ina then output;
	if ina and not inb then SBP = .;
run;

%mend;

/* Apply function */
%create_sbp(sbp_dat, A, 1826.25);
%create_sbp(sbp_dat, B, 1826.25);


/* Check and summarise data */
proc univariate data = datint.varSBP_A; var SBP;run;
proc univariate data = datint.varSBP_B; var SBP;run;

proc means data = datint.varSBP_A n nmiss; var SBP;run;
proc means data = datint.varSBP_B n nmiss; var SBP;run;

