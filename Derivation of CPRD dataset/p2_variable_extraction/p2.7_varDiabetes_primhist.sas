/* Calculate Diabetes */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


/* Define a macro that will extract variables */
/* We extract a simple Yes/NO variable for whether there is a history of each condition */
/* For those who do have history, we extract time since condition was developed */

/* This is effectively %extract_medfile_hist macro, except we don't permanently store the output in datint library */
/* Instead we create three temporary datasets, then combine into one diabetes dataset */ 
%macro extract_medhist_var(variable_out,codenum1,cohort);

/* I now have the medcodes or prodcodes held in all_comor_codes */
/* Sort the codes */
proc sort data = all_comor_codes; by person_id visit_start_date;run;


/* Merge with cohort */
data all_comor_codes1;
	merge all_comor_codes (in=ina) datint.cohort_base (in=inb);
	by person_id;
	if ina and inb then output;
run;


/* Only output events prior to the index date */
data all_comor_codes2;
	set all_comor_codes1;
	if visit_start_date = . then delete;
	if visit_start_date > study_dtindex_&cohort then delete;
run;


/* Make it one line per patient */
/* Note that this will be the first occurence */
data all_comor_codes3;
	set all_comor_codes2;
	by person_id visit_start_date;
	if first.person_id=1 then output;
run;
	

/* Remerge with the cohort */
/* _primhist is the 'history of' variable, 0 = No, 1 = Yes */
/* _primhist_t denotes time that has elapsed since the condition was developed. */
/* For individuals who have no history of the condition, it has been set to -1 */
data var&variable_out._primhist_&cohort (keep = person_id &variable_out._primhist &variable_out._primhist_t);
	merge datint.cohort_base (in=ina) all_comor_codes3 (in=inb);
	by person_id;
	if ina and inb then do;
		&variable_out._primhist_t = study_dtindex_&cohort. - visit_start_date;
		&variable_out._primhist = 1;
	end;
	if ina and ~inb then do;
		&variable_out._primhist_t = -1;
		&variable_out._primhist = 0;
	end;
	if ina then output;
run;
proc delete data= all_comor_codes all_comor_codes1 all_comor_codes2 all_comor_codes3;
/*Delete all the tables created in this program*/
run;

%mend;

*********************************************************************;
/* Diabetes_vague,  */
%extract_medfile_aurum(16);

%extract_medfile_hist(all_comor_codes, Diab_vague, A);
%extract_medfile_hist(all_comor_codes, Diab_vague, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

%put Diabetes_vague,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


*********************************************************************;
/* Diabetes_mellitus_type1,  */
%extract_medfile_aurum(17);

%extract_medfile_hist(all_comor_codes, Diab_t1, A);
%extract_medfile_hist(all_comor_codes, Diab_t1, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

%put Diabetes_mellitus_type1,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;


*********************************************************************;
/* Diabetes_mellitus_type2,  */
%extract_medfile_aurum(18);

%extract_medfile_hist(all_comor_codes, Diab_t2_raw, A);
%extract_medfile_hist(all_comor_codes, Diab_t2_raw, B);

/* Delete all_comor_codes file */
proc delete data = all_comor_codes;run;

%put Diabetes_mellitus_type2,  finished at %sysfunc(time(),timeampm.) on %sysfunc(date(),worddate.).;



**************************************************************************************************************;
/* Want to create the final type 2 diabetes variable now */

/* If a patient has a T2D specific entry, we know they have T2D, and we will set the T2D to T2D_raw */
/* If a patient has a vague diabetes code, and no T1D codes, then we will assume that it is T2D */
/* If a patient has a vague diabetes code and also T1D codes, then we will assume that the vague codes are for the T1D, and not T2D */
/* Note that when doing this step, we want to look back in time for all T1D codes */


/* Therefore create a list of patients who have a type 1 diabetes codes at any point */
%extract_medfile_aurum(17);

/* Sort data and reatin one observation per individual */
data diab_t1_anytime (keep = person_id);
	set all_comor_codes;
	proc sort nodupkey; by person_id;
run;


%macro merge_diabetes(cohort);
data datint.varDiab_primhist_&cohort;
	merge datint.varDiab_vague_primhist_&cohort (in = ina) datint.varDiab_t1_primhist_&cohort (in = inb)
	      datint.varDiab_t2_raw_primhist_&cohort (in = inc) diab_t1_anytime (in = ind);
	by person_id;
	/* Start by letting Diab_t2 = Diab_t2_raw */
	Diab_t2_primhist = Diab_t2_raw_primhist;
	Diab_t2_primhist_t = Diab_t2_raw_primhist_t;
	/* If patients have a vague event, no T2D event, and no T1D event ever, then let Diab_t2 = Diab_vague */
	if Diab_vague_primhist = 1 and Diab_t2_primhist = 0 and not ind then do;
		Diab_t2_primhist = Diab_vague_primhist;
		Diab_t2_primhist_t = Diab_vague_primhist_t;
	end;
	if ina then output;
run;
%mend;

%merge_diabetes(A);
%merge_diabetes(B);

/* Check prevalences */
proc freq data=datint.varDiab_primhist_A;tables Diab_vague_primhist Diab_t1_primhist Diab_t2_raw_primhist 
						Diab_t2_primhist Diab_t1_primhist*Diab_t2_primhist;run;
proc freq data=datint.varDiab_primhist_B;tables Diab_vague_primhist Diab_t1_primhist Diab_t2_raw_primhist 
						Diab_t2_primhist Diab_t1_primhist*Diab_t2_primhist;run;


/* Check distribution of time since developing condition */
proc means data = datint.varDiab_primhist_A; var Diab_t1_primhist_t; where Diab_t1_primhist = 1; run;
proc means data = datint.varDiab_primhist_B; var Diab_t1_primhist_t; where Diab_t1_primhist = 1; run;

proc means data = datint.varDiab_primhist_A; var Diab_t2_primhist_t; where Diab_t2_primhist = 1; run;
proc means data = datint.varDiab_primhist_B; var Diab_t2_primhist_t; where Diab_t2_primhist = 1; run;


/* Delete datasets that were merged into one */
proc delete data = datint.varDiab_vague_primhist_A datint.varDiab_t1_primhist_A datint.varDiab_t2_raw_primhist_A
		   datint.varDiab_vague_primhist_B datint.varDiab_t1_primhist_B datint.varDiab_t2_raw_primhist_B; 
run;