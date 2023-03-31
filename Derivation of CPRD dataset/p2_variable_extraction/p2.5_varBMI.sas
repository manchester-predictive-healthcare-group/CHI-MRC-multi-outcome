/* Calculate BMI */

options ps=60 ls=70;

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

****************************************************************************************************;

/* Start by creating three datasets, one for BMI measurements, one for height, and one for weight */

/* Define a function that will extract all test observations and do the appropriate conversions */


************************************************************;
/* First extract test data for BMI */

/* Extract all medcodes for bmi */
%extract_testfile_aurum(44);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of bmi, and the unit source vales */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;

/* Put into named file */
/* Delete all observations outside of valid range */
data bmi_dat;
	set all_comor_codes;
	if value < 10 or value > 70 or value = . then delete;
run;

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


************************************************;
/* Extract all medcodes for height */
%extract_testfile_aurum(42);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of bmi, and the unit source vales */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;

/* put into named file */
/* Convert value to m, when unit_source_value = cm (122) */
/* Delete extreme values that are likely not a height measurement */

data height_dat;
	set all_comor_codes;
	if unit_source_value = 122 then value = value/100;
	if value > 2.5 or value < 0.8 or value = . then delete;
run;

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


***********************************************;
/* Extract all medcodes for weight */
%extract_testfile_aurum(43);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of bmi, and the unit source vales */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;

/* put into named file */
/* Convert stone to kg */
/* Delete those outside valid range */
data weight_dat;
	set all_comor_codes;
	if unit_source_value = 6265 then value = value*6.35029;
	if value < 30 or value > 300 or value = . then delete;
run;

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

****************************************************************;

/* Write the function that will look for most recent bmi, height and weight value within specified time period */

/* If bmi measurement exists, we will use this */
/* If bmi measurement does not exist, but height and weight measurements do, we estimate bmi from height and weight */
%macro create_bmi(bmi_in, height_in, weight_in, cohort_in, timemeas1_in);

/* Start by using function from var_extract1.sas macro container, which will extract most recent test result */
%extract_latest_test_timeperiod(&bmi_in, bmi_cohort, &timemeas1_in, 0, &cohort_in, bmi_raw);
%extract_latest_test_timeperiod(&height_in, height_cohort, &timemeas1_in, 0, &cohort_in, height);
%extract_latest_test_timeperiod(&weight_in, weight_cohort, &timemeas1_in, 0, &cohort_in, weight);

/* Merge height and weight */
data height_weight;
	merge height_cohort (in = ina) weight_cohort (in = inb);
	by person_id;
	if ina and inb then output;
run;

/* Calculate bmi for individuals who have both */
data height_weight;
	set height_weight;
	bmi_hw = weight/(height**2);
run;

/* Remove any calcualted bmi measurements outside of our valid range */
data height_weight;
	set height_weight;
	if bmi_hw < 10 or bmi_hw > 70 then delete;
run;

/* Now combine the raw bmi values and the bmi values calculated from height/weight */
/* Choose bmi_raw preferentially over height/weight if both exist */
data bmi_both;
	merge bmi_cohort (in = ina) height_weight (in = inb);
	by person_id;
	if ina then BMI = bmi_raw;
	if inb and not ina then BMI = bmi_hw;
run;

/* Merge this with base cohort */
data datint.varBMI_&cohort_in (keep = person_id bmi);
	merge datint.cohort_base (in=ina) bmi_both (in=inb);
	by person_id;
	if ina then output;
	if ina and not inb then BMI = .;
run;

%mend;

/* Apply function */
%create_bmi(bmi_dat, height_dat, weight_dat, A, 1826.25);
%create_bmi(bmi_dat, height_dat, weight_dat, B, 1826.25);

/* Check and summarise data */
proc univariate data = datint.varBMI_A; var bmi;run;
proc univariate data = datint.varBMI_B; var bmi;run;

proc means data = datint.varBMI_A n nmiss; var bmi;run;
proc means data = datint.varBMI_B n nmiss; var bmi;run;


