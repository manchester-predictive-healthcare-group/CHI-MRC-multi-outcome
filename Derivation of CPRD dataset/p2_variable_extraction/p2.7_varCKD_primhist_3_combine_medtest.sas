/* Calculate CKD */
/* First extract individuals with CKD medical codes */
/* Combine this with individuals who have CKD derived through test data, from programs _1 and _2 */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


******************************************************************************;

/* Start with searching medical events */
%extract_medfile_aurum(9);

/* Create the time since event, and histry of indicator */
%extract_medfile_hist(all_comor_codes, CKD_med, A);
%extract_medfile_hist(all_comor_codes, CKD_med, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


******************************************************************************;
/* Now read in the CKD cases from the test data */
/* Note varCKD_1 and varCKD_2 files must have been run in order for this to work */

%macro extract_CKD_test(cohort);

/* Import individuals with CKD from test scores */
proc import datafile="/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p2_variable_extraction/varCKD_primhist_from_test_&cohort..csv" out=varCKD_test_primhist dbms=csv replace;run;

proc contents data = varCKD_test_primhist;run;

/* Sort */
proc sort data = varCKD_test_primhist; by person_id;run;

/* Merge with cohort */
data datint.varCKD_test_primhist_&cohort (keep = person_id CKD_test_primhist_t CKD_test_primhist) ;
	merge datint.cohort_base (in=ina) varCKD_test_primhist (in=inb);
	by person_id;
	if ina and inb then do;
		CKD_test_primhist_t = study_dtindex_&cohort. - measurement_date_num;
		CKD_test_primhist = 1;
	end;
	if ina and ~inb then do;
		CKD_test_primhist_t = -1;
		CKD_test_primhist = 0;
	end;
	if ina then output;
run;

%mend;

%extract_ckd_test(A);
%extract_ckd_test(B);


******************************************************************************;
/* Combine the two into a single dataset */

/* If _test_t and _med_t are both = -1, then CKD_t will be -1 too */
/* If either is > -1, then we pick that */
/* If both are -1, we want whichever one happened first, and was longer ago (i.e. the bigger one) */
/* Therefore we can take the max of test_t and med_t */
/* The max of the censoring variable can also be taken, this is more straight forward, as if either is ^= 0, we set CKD_c = 1 */
%macro combine_med_test(cohort);

data datint.varCKD_primhist_&cohort;
	merge datint.varCKD_med_primhist_&cohort (in=ina) datint.varCKD_test_primhist_&cohort (in=inb);
	by person_id;
	CKD_primhist = max(CKD_med_primhist, CKD_test_primhist);
	CKD_primhist_t = max(CKD_med_primhist_t, CKD_test_primhist_t);
run;

%mend;

%combine_med_test(A);
%combine_med_test(B);


/* Summarise prevalence */
proc freq data=datint.varCKD_primhist_A;tables CKD_primhist CKD_med_primhist CKD_test_primhist;run;
proc freq data=datint.varCKD_primhist_B;tables CKD_primhist CKD_med_primhist CKD_test_primhist;run;


/* Summarise distribution of time since event */
/* Check distribution of time since developing condition */
proc means data = datint.varCKD_primhist_A; var CKD_primhist_t; where CKD_primhist = 1; run;
proc means data = datint.varCKD_primhist_B; var CKD_primhist_t; where CKD_primhist = 1; run;

proc means data = datint.varCKD_med_primhist_A; var CKD_med_primhist_t; where CKD_med_primhist = 1; run;
proc means data = datint.varCKD_med_primhist_B; var CKD_med_primhist_t; where CKD_med_primhist = 1; run;

proc means data = datint.varCKD_test_primhist_A; var CKD_test_primhist_t; where CKD_test_primhist = 1; run;
proc means data = datint.varCKD_test_primhist_B; var CKD_test_primhist_t; where CKD_test_primhist = 1; run;


/* Delete datasets put into the merge */
proc delete data = datint.varCKD_med_primhist_A datint.varCKD_test_primhist_A datint.varCKD_med_primhist_B datint.varCKD_test_primhist_B; run;

