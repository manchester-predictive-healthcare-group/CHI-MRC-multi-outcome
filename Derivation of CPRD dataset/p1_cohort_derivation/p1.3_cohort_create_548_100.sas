/* Create a cohort which requires 548 days prior follow up, and an index date < age 100 */

*/ file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


/* Read in person file */
/* Calculate age turned 65, and age 100, create a study index date (age 65, 1 year valid follow up, study start date) */
data person;
	set datext.persfmt;
	dob = round(((year_of_birth) - 1960)*365.25, 1);
	dt_age65 = dob + round(65*365.25, 1);
	dt_age100 = round(dob + 100*365.25, 1);
	study_dtindex_A = max(regstartdate + 548, &studystartdate, dt_age65);
	format regstartdate best32.;
run;

proc print data = person(obs = 100); var person_id dt_age65 dt_age100 dob; run;


/* Merge with linkage eliiblity file, and retain a marker for if they are elgiible for linkage */
/* Read in linkage eligibility file */
data linkage (keep = person_id hes_e death_e lsoa_e);
	set dathes.linkage_eligibility_fmt;
	by person_id;
run;

/* Merge with linkage eligibility file to get linkage markers, but output all patients */
data person_linkage;
	merge person (in= ina) linkage (in = inb);
	by person_id;
	if ina and inb then linkage_elig = 1;
	if ina and not inb then linkage_elig = 0;
	if ina then output;
run;
proc freq data = person_linkage; tables linkage_elig; run;


/* Merge with ONS data and derive a new death date, minimum of CPRD death date and ONS death date */
/* Read in ONS death date and sort by person_id */
data ons_death;
	set dathes.onsdeathfmt;
run;

/* Output only one observation per person, and get rid of excess variables  */
data ons_death (keep = person_id death_date_ons);
	set ons_death;
	by person_id;
	if first.person_id then output;
run;


/* Merge ons_death_date file by person file */
data temp;
	merge person_linkage (in = ina) ons_death (in = inb);
	by person_id;
	if ina then output;
run;


/* Create a new death_date, which is the minimum of death_date from primary care and ons */
/* Create a new dtcens, which is minimum of dtcens and the new deathdate */
data temp2;
	set temp;
	format death_date_comb dtcens_combdeath ddmmyy10.;
	death_date_comb = min(death_date, death_date_ons);
	dtcens_combdeath = min(dtcens, death_date_comb);
run;



/* Remove anyone for which index_date date is not prior to dtcens_comb, or prior to them being aged 100 */
data temp3;
	set temp2;
	if study_dtindex_A < dtcens_combdeath and study_dtindex_A < dt_age100 then output;
run;


/* I actually want to define two study_dtindex's */
/* The second one is at a random point between study_dtindex and min(dtcens_combdeath, age turned 100) */
/* Also format the index dates */
data temp4;
	set temp3;
	format study_dtindex_A study_dtindex_B ddmmyy10.;
	*Create random U[0,1] variable;
	u_rand = rand("Uniform");
	*Multiply this by the length of the observation period for each individual;
	u_rand2 = u_rand*(min(dtcens_combdeath, dt_age100) - study_dtindex_A);
	*Add this onto study_dtindex, and take the floor, to get study_dtindex2;
	*We choose floor so that study_dtindex_B is never equal to dtcens_combdeath;
	study_dtindex_B = floor(study_dtindex_A + u_rand2);
	drop u_rand u_rand2;
run;

proc print data = temp4 (obs = 50);run;

/* Finally, create a dataset which is permanent (i.e. save the cohort) */
data datint.cohort_base_548_100;
	set temp4;
run;