/* Produce histograms for continuous variables in this program */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

%let var_cont = Age BMI SBP Cholhdl_ratio;


/* Cycle through histograms to a word document */
%macro cycle_hisotgrams(gender1,cohort1);
%local i next_name;
%let i=1;
%do %while (%scan(&var_cont, &i) ne );
   %let next_name = %scan(&var_cont, &i);
	proc univariate data=datint.cohort_var_&cohort1.;var &next_name;histogram;run;
   %let i = %eval(&i + 1);
%end;
%mend;


/*
%macro cycle_hisotgrams(gender1,cohort1);

proc univariate data=datint.cohort_var_&cohort1.;var age;where gender = &gender1;histogram;run;
proc univariate data=datint.cohort_var_&cohort1.;var BMI;where gender = &gender1;histogram;run;
proc univariate data=datint.cohort_var_&cohort1.;var SBP;where gender = &gender1;histogram;run;
proc univariate data=datint.cohort_var_&cohort1.;var Cholhdl_ratio;where gender = &gender1;histogram;run;

%mend;
*/

/* Generate histograms for continuous variables - cohort A */
ODS rtf file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_hist_male_A.rtf";
%cycle_hisotgrams(0,A);
ODS rtf close;


ODS rtf file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_hist_female_A.rtf";
%cycle_hisotgrams(1,A);
ODS rtf close;


/* Generate histograms for continuous variables - cohort B */

ODS rtf file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_hist_male_B.rtf";
%cycle_hisotgrams(0,B);
ODS rtf close;


ODS rtf file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/checks_hist_female_B.rtf";
%cycle_hisotgrams(1,B);
ODS rtf close;

