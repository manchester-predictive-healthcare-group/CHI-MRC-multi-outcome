/* Calculate Ethnicity */

options ps=60 ls=70;

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

****************************************************************************************************;

/* Extract all medcodes for Smoking from the medical file */

/* Medical file */
%extract_medfile_aurum(41);

/* sort the data */
proc sort data = all_comor_codes;by person_id visit_start_date; run;

/* Output one observation per person */
data all_comor_codes;
	set all_comor_codes;
	by person_id visit_start_date;
	if last.person_id then output;
run;

/* Read in the ethnicity dataset with the description of each medical code */
%importcsv(ethnicity_codes ,&ffile_code.ethnicity_aurum_wgroups_mcid, ddatarow = 2);

/* Change medcodeid to be visit_concept_id, rename output variables, and retain only relevant variables */
data ethnicity_codes (keep = visit_concept_id Ethnicity16 Ethnicity6);
	set ethnicity_codes;
	visit_concept_id = input(medcodeid, best32.);
	Ethnicity16 = input(group_16_num, best32.);
	Ethnicity6 = input(group_6_num, best32.);
	format visit_concept_id best32.;
run;

/* Sort both datasets by visit_concept_id */
proc sort data = all_comor_codes;by visit_concept_id;run;
proc sort data = ethnicity_codes;by visit_concept_id;run;

/* Merge extracted ethnicity codes with the code list */
data ethnicity_comb_codes;
	merge all_comor_codes (in = ina) ethnicity_codes (in = inb);
	by visit_concept_id;
	if ina then output;
run;


/* Sort this dataset by person_d, and visit_start_id */
proc sort data = ethnicity_comb_codes; by person_id visit_start_date;run;
proc print data = ethnicity_comb_codes (obs = 100);run;
proc contents data = ethnicity_comb_codes;run;

%macro create_Ethnicity(ethnicity_in, cohort_in);

/* Combine both these datasets with the base cohort to create variables of interest */
data datint.varEthnicity_&cohort_in (keep = person_id Ethnicity16 Ethnicity6);
	merge datint.cohort_base (in=ina) &ethnicity_in (in=inb);
	by person_id;
	if ina and not inb then do;
		Ethnicity16 = .;
		Ethnicity6 = .;
	end;
	if ina then output;
run;

%mend;

%create_Ethnicity(ethnicity_comb_codes, A);
%create_Ethnicity(ethnicity_comb_codes, B);


/* Check and summarise data */
proc freq data = datint.varEthnicity_A; tables Ethnicity16 Ethnicity6;run;
proc freq data = datint.varEthnicity_B; tables Ethnicity16 Ethnicity6;run;

proc means data = datint.varEthnicity_A n nmiss; var Ethnicity16 Ethnicity6;run;
proc means data = datint.varEthnicity_B n nmiss; var Ethnicity16 Ethnicity6;run;
