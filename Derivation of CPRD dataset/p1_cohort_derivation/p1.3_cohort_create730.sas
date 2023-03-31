/* Create a cohort (temporarily, until I've finalised how to derive it) */

*/ file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

/* Calculate age turned 65, create a study index date (age 65, 1 year valid follow up, study start date) */
/* and remove anyone for which this date is not prior to dtcens */
data temp;
	set datext.persfmt;
	dob = ((year_of_birth) - 1960)*365.25;
	dt_age65 = dob + 65*365.25;
	study_dtindex_A = max(regstartdate + 730, &studystartdate, dt_age65);
	if study_dtindex_A < dtcens then output temp;
	format regstartdate best32.;
run;



/* I actually want to define two study_dtindex's */
/* The second one is at a random point between study_dtindex and dtcens */
data temp2;
	set temp;
	*Create random U[0,1] variable;
	u_rand = rand("Uniform");
	*Multiply this by the length of the observation period for each individual;
	u_rand2 = u_rand*(dtcens - study_dtindex_A);
	*Add this onto study_dtindex, and take the floor, to get study_dtindex2;
	study_dtindex_B = floor(study_dtindex_A + u_rand2);
	drop u_rand u_rand2;
run;

proc print data = temp2 (obs=50);run;

/* Also wnat to create date versions of study_dtindex_A and study_dtindex_B, and plot histograms */
data temp3;
	set temp2;
	format study_dtindex_A study_dtindex_B mmddyy10.;
run;

proc print data = temp3 (obs = 50);run;

/*
proc univariate data=temp3;var study_dtindex_A study_dtindex_B;histogram;run;
*/

/* Finally, create a dataset which is permanent (i.e. save the cohort) */
data datint.cohort_base730;
	set temp3;
run;


proc print data = datint.cohort_base730 (obs = 20);run;