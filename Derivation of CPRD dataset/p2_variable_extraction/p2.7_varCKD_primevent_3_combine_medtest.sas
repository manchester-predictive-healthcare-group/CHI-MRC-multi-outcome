/* Extract outcome variables for CKD */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

******************************************************************************;
/* Define a macro to: */
/* 1) Create time until first event identified via medical codes */
/* 2) Create time until first event identified through test data */
/* 3) Create a combined variable which is time until first of either of these */

********************************************************************************;

/* Start with searching medical events */
%extract_medfile_aurum(9);

/* Create the time until event, and censoring indicator */
%extract_medfile_event(all_comor_codes, CKD_med, A);
%extract_medfile_event(all_comor_codes, CKD_med, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


********************************************************************************;
/* Next identify events through test data */

%macro extract_ckd_testfile_event(cohort);

/* Import individuals with CKD from test scores */
proc import datafile="/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primevent_from_test_&cohort..csv" out=varCKD_test_primev dbms=csv replace;run;

proc contents data = varCKD_test_primev;run;

/* Sort */
proc sort data = varCKD_test_primev; by person_id measurement_date_num;run;

/* Only output the first observation for each person, and retain person_id and measurement_date_num */
/* I actually think the R code only outputs the first event for each person anyway, making this step redundant, but will leave in now to double check */
data varCKD_test_primev;
	set varCKD_test_primev;
	by person_id measurement_date_num;
	if first.person_id then output;
run;

/* Merge with cohort */
/* The gfr measurements were restricted to be prior to dtcens_combdeath, therefore none of these events should be happening after dtcens_combdeath, so I don't need to do anything about this */ 
data datint.varCKD_test_primev_&cohort (keep = person_id CKD_test_primev_t CKD_test_primev_c) ;
	merge datint.cohort_base (in=ina) varCKD_test_primev (in=inb);
	by person_id;
	if ina and inb then do;
		CKD_test_primev_t = measurement_date_num - study_dtindex_&cohort.;
		CKD_test_primev_c = 1;
	end;
	if ina and ~inb then do;
		CKD_test_primev_t = dtcens_combdeath - study_dtindex_&cohort.;
		CKD_test_primev_c = 0;
	end;
	if ina then output;
run;

%mend;

%extract_ckd_testfile_event(A);
%extract_ckd_testfile_event(B);



********************************************************************************;
/* Combine the two into a single dataset */

%macro combine_med_test(cohort);

data datint.varCKD_primev_&cohort;
	merge datint.varCKD_med_primev_&cohort (in = ina) datint.varCKD_test_primev_&cohort (in = inb);
	by person_id;
	CKD_primev_t = min(CKD_med_primev_t, CKD_test_primev_t);
	CKD_primev_c = max(CKD_med_primev_c, CKD_test_primev_c);
run;

%mend;

%combine_med_test(A);
%combine_med_test(B);

proc means data = datint.varCKD_primev_A; var CKD_med_primev_t CKD_test_primev_t CKD_primev_t;run;
proc means data = datint.varCKD_primev_B; var CKD_med_primev_t CKD_test_primev_t CKD_primev_t;run;

proc freq data = datint.varCKD_primev_A; tables CKD_med_primev_c CKD_test_primev_c CKD_primev_c;run;
proc freq data = datint.varCKD_primev_B; tables CKD_med_primev_c CKD_test_primev_c CKD_primev_c;run;


/* Delete datasets put into the merge */
proc delete data = datint.varCKD_med_primev_A datint.varCKD_test_primev_A datint.varCKD_med_primev_B datint.varCKD_test_primev_B; run;
