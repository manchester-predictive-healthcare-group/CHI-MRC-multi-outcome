/* Calculate Age */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


/* Define the macro */
%macro extract_age(cohort);
	
/* Remerge with the cohort, creating a 0/1 variable for the presence of the comorbidity */
data datint.varAge_&cohort (keep = person_id Age);
	set datint.cohort_base;
	Age = (study_dtindex_&cohort - dob)/365.25;
run;

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%mend;

/* These are the only variables where we simply look for a medcode at some point n the past */

/* Atrial fibrillation */
%extract_age(A);
%extract_age(B);


title 'male';
proc univariate data=datint.varAge_A;var Age;run;
proc univariate data=datint.varAge_B;var Age;run;
