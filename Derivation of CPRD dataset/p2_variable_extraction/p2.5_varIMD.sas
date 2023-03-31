/* Calculate IMD */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


/* Define the macro */
%macro extract_imd(cohort);
	
/* Remerge with the cohort, creating a 0/1 variable for the presence of the comorbidity */
data datint.varIMD_&cohort (keep = person_id IMD);
	merge datint.cohort_base (in = ina) dathes.imd_patient_fmt (in = inb);
	by person_id;
	IMD = imd2015_5;
	if ina then output;
run;

proc delete data=work._all_;
/*Delete all the tables created in this program*/
run;

%mend;

/* These are the only variables where we simply look for a medcode at some point n the past */

/* Atrial fibrillation */
%extract_imd(A);
%extract_imd(B);


title 'male';
proc freq data=datint.varIMD_A;tables IMD;run;
proc freq data=datint.varIMD_B;tables IMD;run;



data test_imd (keep = person_id IMD);
	merge datint.cohort_base (in = ina) dathes.imd_patient_fmt (in = inb);
	by person_id;
	IMD = imd2015_5;
	if ina and inb then output;
run;


data test_linkage;
	merge datint.cohort_base (in = ina) dathes.linkage_eligibility_fmt (in = inb);
	by person_id;
	if ina and inb then output;
run;

