/* Calculate Cholesterol */

options ps=60 ls=70;

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

****************************************************************************************************;

/* Start by creating three datasets, one for BMI measurements, one for height, and one for weight */

/* Define a function that will extract all test observations and do the appropriate conversions */


******************************************************************;
/* First extract test data for cholesteriol/HDL ratio */

/* Extract all medcodes for chol/hdl ratio */
%extract_testfile_aurum(49);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of bmi, and the unit source vales */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;


/* Put into named file */
/* Delete all observations outside of valid range */
data ratio_dat;
	set all_comor_codes;
	if value < 1 or value > 15 or value = . then delete;
run;

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;



*********************************************************************;
/* Extract all medcodes for hdl only */
%extract_testfile_aurum(48);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of bmi, and the unit source vales */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;


/* Put into named file */
/* Delete all observations outside of valid range */
data hdl_dat;
	set all_comor_codes;
	if value < 0.39 or value > 5 or value = . then delete;
run;

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


*********************************************************************;
/* Extract all medcodes for total chol */
%extract_testfile_aurum(47);

/* sort the data */
proc sort data = all_comor_codes;by person_id measurement_date; run;

/* look at the distribution of bmi, and the unit source vales */
proc freq data = all_comor_codes; tables unit_source_value;run;
proc univariate data = all_comor_codes; var value;run;


/* Put into named file */
/* Delete all observations outside of valid range */
data total_dat;
	set all_comor_codes;
	if value < 1.3 or value > 14 or value = . then delete;
run;

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


************************************************************************;
/* Write the function that will look for most recent chol/hdl ratio, hdl and total cholesterol value within specified time period */

/* If chol/hdl ratio measurement exists, we will use this */
/* If chol/hdl measurement does not exist, but total chol and hdl measurements do, we estimate ratio from total and hdl */
%macro create_cholhdl(ratio_in, hdl_in, total_in, cohort_in, timemeas1_in);

/* Start by using function from var_extract1.sas macro container, which will extract most recent test result */
%extract_latest_test_timeperiod(&ratio_in, ratio_cohort, &timemeas1_in, 0, &cohort_in, ratio);
%extract_latest_test_timeperiod(&hdl_in, hdl_cohort, &timemeas1_in, 0, &cohort_in, hdl);
%extract_latest_test_timeperiod(&total_in, total_cohort, &timemeas1_in, 0, &cohort_in, total);

/* Merge hdl and total */
data hdl_total;
	merge hdl_cohort (in = ina) total_cohort (in = inb);
	by person_id;
	if ina and inb then output;
run;

/* Calculate cholhdl_ratio for individuals who have both */
data hdl_total;
	set hdl_total;
	ratio_est = total/hdl;
	if ratio_est < 1 or ratio_est > 15 then delete;
run;

/* Now combine the raw bmi values and the bmi values calculated from height/weight */
/* Choose raw ratio preferentially over the ratio estimated from hdl and total if both are available */
data cholratio_both;
	merge ratio_cohort (in = ina) hdl_total (in = inb);
	by person_id;
	if ina then Cholhdl_ratio = ratio;
	if inb and not ina then Cholhdl_ratio = ratio_est;
run;

/* This checks whether the estimates ratios, and the ratios taken straight from database, have a similar distribution */
/* They have very similar median/mean, which is good */
proc univariate data = cholratio_both; var ratio ratio_est;run;

/* Merge this with base cohort */
data datint.varCholhdl_ratio_&cohort_in (keep = person_id Cholhdl_ratio);
	merge datint.cohort_base (in=ina) cholratio_both (in=inb);
	by person_id;
	if ina then output;
	if ina and not inb then Cholhdl_ratio = .;
run;

%mend;



/* Apply function */
%create_cholhdl(ratio_dat, hdl_dat, total_dat, A, 1826.25);
%create_cholhdl(ratio_dat, hdl_dat, total_dat, B, 1826.25);

/* Check and summarise data */
proc univariate data = datint.varCholhdl_ratio_A; var Cholhdl_ratio;run;
proc univariate data = datint.varCholhdl_ratio_B; var Cholhdl_ratio;run;

proc means data = datint.varCholhdl_ratio_A n nmiss; var Cholhdl_ratio;run;
proc means data = datint.varCholhdl_ratio_B n nmiss; var Cholhdl_ratio;run;

