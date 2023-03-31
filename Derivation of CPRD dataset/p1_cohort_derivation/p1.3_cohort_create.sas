/* Create a cohort (temporarily, until I've finalised how to derive it) */

*/ file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

/* Calculate age turned 65, and age 85, create a study index date (age 65, 1 year valid follow up, study start date) */
/* Remove anyone for which this date is not prior to dtcens, or prior to them being aged 85 */
data temp;
	set datext.persfmt;
	dob = ((year_of_birth) - 1960)*365.25;
	dt_age65 = dob + 65*365.25;
	dt_age85 = dob + 85*365.25;
	study_dtindex_A = max(regstartdate + 548, &studystartdate, dt_age65);
	if study_dtindex_A < dtcens and study_dtindex_A < dt_age85 then output temp;
	format regstartdate best32.;
run;



/* I actually want to define two study_dtindex's */
/* The second one is at a random point between study_dtindex and min(dtcens, age turned 85) */
/* Also format the index dates */
data temp2;
	set temp;
	format study_dtindex_A study_dtindex_B ddmmyy10.;
	*Create random U[0,1] variable;
	u_rand = rand("Uniform");
	*Multiply this by the length of the observation period for each individual;
	u_rand2 = u_rand*(min(dtcens, dt_age85) - study_dtindex_A);
	*Add this onto study_dtindex, and take the floor, to get study_dtindex2;
	*We choose floor so that study_dtindex_B is never equal to dtcens;
	study_dtindex_B = floor(study_dtindex_A + u_rand2);
	drop u_rand u_rand2;
run;

proc print data = temp2 (obs = 50);run;

/* Finally, create a dataset which is permanent (i.e. save the cohort) */
data datint.cohort_base;
	set temp2;
run;