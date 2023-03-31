/* Calculate Smoking status */

options ps=60 ls=70;

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

****************************************************************************************************;

%macro calc_miss(gender, cohort);
proc means data = datint.cohort_var_&cohort stackods n nmiss NWAY;
	where gender = &gender;
	var BMI SBP Cholhdl_ratio;
	ods output summary = summary_stats_cont;
run;

data percent_missing_cont;
	set summary_stats_cont;
	Pct_Missing = nmiss/sum(n, nmiss);
run; 


proc freq data = datint.cohort_var_&cohort;
	where gender = &gender;
	tables Ethnicity6 Smoking IMD / missing;
	ods output OneWayFreqs = summary_stats_disc;
run;

data summary_stats_disc2;
	set summary_stats_disc;
	rownum = _N_;
run;

proc sort data = summary_stats_disc2; by Table rownum;run;

data percent_missing_disc (keep = Table percent);
	set summary_stats_disc2;
	by Table rownum;
	if first.Table then output;
run;

%mend;

/* Calculate missingness for females - cohort A */
%calc_miss(1,A);
ODS pdf file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_missing_female_A.pdf";
proc print data = percent_missing_cont;
run;
proc print data = percent_missing_disc;
run;
ODS pdf close;

/* Calculate missingness for males - cohort A */
%calc_miss(0,A);
ODS pdf file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_missing_male_A.pdf";
proc print data = percent_missing_cont;
run;
proc print data = percent_missing_disc;
run;
ODS pdf close;

/* Calculate missingness for females - cohort B */
%calc_miss(1,B);
ODS pdf file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_missing_female_B.pdf";
proc print data = percent_missing_cont;
run;
proc print data = percent_missing_disc;
run;
ODS pdf close;

/* Calculate missingness for males - cohort B */
%calc_miss(0,B);
ODS pdf file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_missing_male_B.pdf";
proc print data = percent_missing_cont;
run;
proc print data = percent_missing_disc;
run;
ODS pdf close;
