/* Create a cohort (temporarily, until I've finalised how to derive it) */

*/ file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

/* Choose the cohort which we are using for extraction (multiple were created, based on different maximum age) */
data datint.cohort_base;
	set datint.cohort_base_548_85;
run;

proc contents data = datint.cohort_base;run;
proc print data = datint.cohort_base (obs = 100);var person_id death_date death_date_ons death_date_comb;run;


/* Test merging this with linkage eligibility */
/*
data temp;
	merge datint.cohort_base (in = ina) dathes.linkage_eligibility_fmt (in = inb);
	by person_id;
	if ina and inb then output;
run;

data temp2;
	set temp;
	if hes_e and death_e = 1 then output;
run;
*/