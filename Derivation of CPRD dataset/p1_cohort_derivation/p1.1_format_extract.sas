/* Exploratory code*/

*/ file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

/* Write a macro to convert each of the files which have multiple numbers */
/* Will just write code for person file and practice file, as there is just one file */

/* Macro for convertin medical files */
%macro convert_med;

%do filenum=1 %to 8;

/* Do formatting */
data datext.allmedfmt_&filenum.;
	set datext.allmed_&filenum. (rename=(visit_concept_id = visit_concept_id_old person_id = person_id_old));
	visit_concept_id = input(visit_concept_id_old, best32.);
	person_id = input(person_id_old, best32.);
	format visit_concept_id person_id best32. visit_start_date ddmmyy10.;
	drop visit_concept_id_old person_id_old;
run;

/* Sort by person then by event date */
proc sort data = datext.allmedfmt_&filenum.; by person_id visit_start_date;run;

%end;

%mend;


/* Run the macro */
/*
%convert_med;
*/


/* Macro for converting prescription files */
%macro convert_rx;

%do filenum=1 %to 8;

/* Do formatting */
data datext.allrxfmt_&filenum.;
	set datext.allrx_&filenum. (rename=(drug_concept_id = drug_concept_id_old person_id = person_id_old));
	drug_concept_id = input(drug_concept_id_old, best32.);
	person_id = input(person_id_old, best32.);
	format drug_concept_id person_id best32. drug_exposure_start_date ddmmyy10.;
	drop drug_concept_id_old person_id_old;
run;

/* Sort by person then by event date */
proc sort data = datext.allrxfmt_&filenum.; by person_id drug_exposure_start_date;run;

%end;

%mend;

/* Run the macro */
/*
%convert_rx;
*/


/* Macro for converting prescription files */
%macro convert_test;

%do filenum=1 %to 8;

/* Do formatting */
data datext.alltestfmt_&filenum.;
	set datext.alltest_&filenum. (rename=(measurement_concept_id = measurement_concept_id_old person_id = person_id_old));
	measurement_concept_id = input(measurement_concept_id_old, best32.);
	person_id = input(person_id_old, best32.);
	format measurement_concept_id person_id best32. measurement_date ddmmyy10.;
	drop measurement_concept_id_old person_id_old;
run;

/* Sort by person then by event date */
proc sort data = datext.alltestfmt_&filenum.; by person_id measurement_date;run;

%end;

%mend;

/* Run the macro */
/*
%convert_test;
*/


/* Convert the person file */
/*
data datext.persfmt;
	set datext.pers (rename = (care_site_id = care_site_id_old person_id = person_id_old));
	care_site_id = input(care_site_id_old, best32.);
	person_id = input(person_id_old, best32.);
	format care_site_id person_id best32. death_date dtcens dtindex dtindex_hes dtvalid observation_period_end_date observation_period_end_hes
		observation_period_start_date observation_period_start_hes ddmmyy10.;
	drop person_id_old care_site_id_old;
run;
proc sort data = datext.persfmt; by person_id;run;
*/


/* Convert the practice file */
/*
data datext.practicefmt;
	set datext.practice;
	format observation_care_site_end_date observation_care_site_start_date ddmmyy10.;
run;
*/


/* Convert the ons_death file */
/*
data dathes.onsdeathfmt;
	set dathes.onsdeath (rename = (care_site_id = care_site_id_old person_id = person_id_old));
	care_site_id = input(care_site_id_old, best32.);
	person_id = input(person_id_old, best32.);
	format care_site_id person_id best32. death_date_ons ddmmyy10.;
	drop person_id_old care_site_id_old;
run;
proc sort data = dathes.onsdeathfmt; by person_id;run;
*/


/* Convert the hes primary diagnosis file */
/*
data dathes.hes_diagnosis_primary_fmt (keep = person_id admidate ICD_PRIMARY ICDx icd_primary3 icd3 icd5 description_icd3 description_icd5);
	set dathes.hes_diagnosis_primary (rename = (person_id = person_id_old));
	person_id = input(person_id_old, best32.);
	format person_id best32.;
run;
proc sort data = dathes.hes_diagnosis_primary_fmt; by person_id admidate;run;
*/


/* Convert the hes patient file */
/*
proc freq data = dathes.hes_patient; tables gen_ethnicity;run;

data dathes.hes_patient_fmt (keep = person_id gen_ethnicity);
	set dathes.hes_patient (rename = (person_id = person_id_old));
	person_id = input(person_id_old, best32.);
	format person_id best32.;
run;
proc sort data = dathes.hes_patient_fmt; by person_id;run;
*/

/* Convert the linkage_eligibilty file */
/*
data dathes.linkage_eligibility_fmt;
	set dathes.linkage_eligibility (rename = (care_site_id = care_site_id_old person_id = person_id_old));
	care_site_id = input(care_site_id_old, best32.);
	person_id = input(person_id_old, best32.);
	format care_site_id person_id best32.;
	drop care_site_id_old person_id_old;
run;
proc sort data = dathes.linkage_eligibility_fmt; by person_id;run;
*/


/* Convert the imd patient file */
data dathes.imd_patient_fmt;
	set dathes.imd_patient (rename = (care_site_id = care_site_id_old person_id = person_id_old));
	care_site_id = input(care_site_id_old, best32.);
	person_id = input(person_id_old, best32.);
	format care_site_id person_id best32.;
	drop care_site_id_old person_id_old;
run;
proc sort data = dathes.imd_patient_fmt; by person_id;run;


ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p1_prelim/exploratory.pdf";

proc contents data = datext.persfmt;run;

proc print data = datext.persfmt (obs = 50);run;

proc contents data = datext.allmedfmt_1;run;

proc print data = datext.allmedfmt_1 (obs = 50);run;

proc contents data = datext.allrxfmt_1;run;

proc print data = datext.allrxfmt_1 (obs = 50);run;

proc contents data = datext.alltestfmt_1;run;

proc print data = datext.alltestfmt_1 (obs = 50);run;

proc contents data = datext.persfmt;run;

proc print data = datext.practicefmt (obs = 50);run;

proc contents data = dathes.onsdeathfmt; run;

proc print data = dathes.onsdeathfmt (obs = 50);run;


ODS PDF close;