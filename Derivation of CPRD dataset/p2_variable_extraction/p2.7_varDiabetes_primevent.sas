/* Extract outcome variables for Diabetes type 2 */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


*************************************************************************;
/* Extract event, time until and censoring indicator, for diabetes vague, type 1, and type 2 */

***************************************************;
/* Vague diabetes codes */
%extract_medfile_aurum(16);

%extract_medfile_event(all_comor_codes, Diab_vague, A);
%extract_medfile_event(all_comor_codes, Diab_vague, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


***************************************************;
/* Diabetes type 1 */
%extract_medfile_aurum(17);

%extract_medfile_event(all_comor_codes, Diab_t1, A);
%extract_medfile_event(all_comor_codes, Diab_t1, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;


***************************************************;
/* Diabetes type 2 */
%extract_medfile_aurum(18);

%extract_medfile_event(all_comor_codes, Diab_t2_raw, A);
%extract_medfile_event(all_comor_codes, Diab_t2_raw, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;



******************************************************************************;
/* Want to create the final type 2 diabetes variable now */

/* If a patient has a T2D specific entry, we know they have T2D, and we will set the T2D to T2D_raw */
/* If a patient has a vague diabetes code, and no T1D codes, then we will assume that it is T2D */
/* If a patient has a vague diabetes code and also T1D codes, then we will assume that the vague codes are for the T1D, and not T2D */
/* Note that when doing this step, we want to look back in time for all T1D codes, as if someone already has T1D, */
/* but then only have vague codes after the index date, we will assume all these vague codes are still for the T1D */

/* Therefore create a list of patients who have a type 1 diabetes codes at any point */
%extract_medfile_aurum(17);

/* Sort data and reatin one observation per individual */
data diab_t1_anytime (keep = person_id);
	set all_comor_codes;
	proc sort nodupkey; by person_id;
run;


/* Now we want to merge these datasets */

%macro combine_diabetes(cohort);
data datint.varDiab_primev_&cohort.;
	merge datint.varDiab_vague_primev_&cohort (in = ina) datint.varDiab_t1_primev_&cohort (in = inb)
	      datint.varDiab_t2_raw_primev_&cohort (in = inc) diab_t1_anytime (in = ind);
	by person_id;
	/* Start by letting Diab_t2 = Diab_t2_raw */
	Diab_t2_primev_t = Diab_t2_raw_primev_t;
	Diab_t2_primev_c = Diab_t2_raw_primev_c;
	/* If patients have a vague event, and no T1D event ever, then let Diab_t2 = Diab_vague */
	if Diab_t2_raw_primev_c = 0 and Diab_vague_primev_c = 1 and not ind then do;
		Diab_t2_primev_t = Diab_vague_primev_t;
		Diab_t2_primev_c = Diab_vague_primev_c;
	end;
	if ina then output;
run;

%mend;

%combine_diabetes(A);
%combine_diabetes(B);

proc means data = datint.varDiab_primev_A; var Diab_vague_primev_t Diab_t1_primev_t Diab_t2_raw_primev_t Diab_t2_primev_t;run;
proc means data = datint.varDiab_primev_B; var Diab_vague_primev_t Diab_t1_primev_t Diab_t2_raw_primev_t Diab_t2_primev_t;run;

proc freq data = datint.varDiab_primev_A; tables Diab_vague_primev_c Diab_t1_primev_c Diab_t2_raw_primev_c Diab_t2_primev_c Diab_t1_primev_c*Diab_t2_primev_c;run;
proc freq data = datint.varDiab_primev_B; tables Diab_vague_primev_c Diab_t1_primev_c Diab_t2_raw_primev_c Diab_t2_primev_c Diab_t1_primev_c*Diab_t2_primev_c;run;



/* Delete datasets that were merged into one */
proc delete data = datint.varDiab_vague_primev_A datint.varDiab_t1_primev_A datint.varDiab_t2_raw_primev_A
		   datint.varDiab_vague_primev_B datint.varDiab_t1_primev_B datint.varDiab_t2_raw_primev_B; 
run;

