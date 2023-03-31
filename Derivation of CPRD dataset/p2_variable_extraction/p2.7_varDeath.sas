/* Calculate time until death and a censoring indicator */
/* Also want to calculate the Nelson Aalen estimate of risk of death at either time of death, or time of censoring (this will be used in imputation) */
/* Note that if we exclude those without HES linkage, I'm going to have to re-calculate this variable */


/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


/* Write a macro to calculate time until death and a censoring indicator */
%macro extract_death(cohort);

data temp_death_&cohort (keep = person_id death_date_comb dtcens dtcens_combdeath Death_t Death_c study_dtindex_&cohort regstartdate dt_age65);
	set datint.cohort_base;
	Death_t = dtcens_combdeath - study_dtindex_&cohort;
	/* If death date is either before, or 30 days after dtcens, we assume they died and assign indicator = 1 */
	if death_date_comb = . then Death_c = 0;
	if death_date_comb ^= . then do;
		if death_date_comb <= dtcens + 30 then Death_c = 1;
		if death_date_comb > dtcens + 30 then Death_c = 0;
	end;
run;

%mend;

%extract_death(A);
%extract_death(B);

proc univariate data = temp_death_A;var Death_t;run;
proc univariate data = temp_death_B;var Death_t;run;

proc freq data = temp_death_A; tables Death_c; run;
proc freq data = temp_death_B; tables Death_c; run;



/* Now I want to calulcate the Nelson Aalen estimator of death */
/* I will have to do this twice, for the whole cohort, and for the cohort of people with linkage */

/* Do the entire cohort first */
%macro NelsonAalen(cohort);

/* Create a temp dataset to fit the model to */
data temp;
	set temp_death_&cohort;
run;

/* Calculate NelsonAalen estimator at each follow up time */
ods output ProductLimitEstimates=ProdLimEst;
proc lifetest data=temp outsurv=Outsurv nelson METHOD=PL;
	time Death_t*Death_c(0);
run;

/* Put output into a dataset */
data NelsonAalenScores (keep = Death_t CumHaz);
	set ProdLimEst;
run;

/* For time points where there is not an event, retain the cumhaz from before */
data NelsonAalenScores2;
	set NelsonAalenScores;
	retain CumHaz_retain;
	if CumHaz ^= . then CumHaz_retain = Cumhaz;
	output;
run;

proc print data = NelsonAalenScores2 (obs = 50); run; 

/* For time points with multiple rows, retain one only */
data NelsonAalenScores3;
	set NelsonAalenScores2;
	by Death_t;
	if first.Death_t then output;
run;

data datint.NelsonAalenScores3_&cohort;
	set NelsonAalenScores3;
run;

proc print data = NelsonAalenScores3 (obs = 50); run; 

/* Sort the death variable dataset by Death_t */
proc sort data=temp; by Death_t; run;

/* Merge with the Nelson Aalen estimators */ 
data datint.varDeath_&cohort (keep = person_id Death_t Death_c CumHaz_retain rename = (CumHaz_retain = Death_NelsonAalen));
	merge temp (in=ina) NelsonAalenScores3 (in=inb);
	by Death_t;
	if ina then output;
run;

/* Sort by person_id again */
proc sort data=datint.varDeath_&cohort; by person_id; run;

%mend;

%NelsonAalen(A);
%NelsonAalen(B);

proc print data = datint.varDeath_A(obs = 100); run;

proc means data = datint.varDeath_A mean median min max n nmiss; var Death_NelsonAalen; run; 
proc means data = datint.varDeath_B mean median min max n nmiss; var Death_NelsonAalen; run; 
proc means data = datint.varDeath_B mean median min max n nmiss; var Death_NelsonAalen; where Death_t > 0; run; 



/* Do the entire cohort first */
%macro NelsonAalen_link(cohort);

/* Create a temp dataset to fit the model to */
data temp;
	set temp_death_&cohort;
run;

/* Merge with the linkage eligibility file and only output those with linkge */
data temp2;
	merge temp (in = ina) dathes.linkage_eligibility_fmt (in = inb);
	by person_id;
	if ina and inb then output;
run;

/* Calculate NelsonAalen estimator at each follow up time */
ods output ProductLimitEstimates=ProdLimEst;
proc lifetest data=temp2 outsurv=Outsurv nelson METHOD=PL;
	time Death_t*Death_c(0);
run;

/* Put output into a dataset */
data NelsonAalenScores (keep = Death_t CumHaz);
	set ProdLimEst;
run;

/* For time points where there is not an event, retain the cumhaz from before */
data NelsonAalenScores2;
	set NelsonAalenScores;
	retain CumHaz_retain;
	if CumHaz ^= . then CumHaz_retain = Cumhaz;
	output;
run;

proc print data = NelsonAalenScores2 (obs = 50); run; 

/* For time points with multiple rows, retain one only */
data NelsonAalenScores3;
	set NelsonAalenScores2;
	by Death_t;
	if first.Death_t then output;
run;

data datint.NelsonAalenScores3_link_&cohort;
	set NelsonAalenScores3;
run;

proc print data = NelsonAalenScores3 (obs = 50); run; 

/* Sort the death variable dataset, which only contains those with linkage, by Death_t */
proc sort data=temp2; by Death_t; run;

/* Merge with the Nelson Aalen estimators */ 
data temp3 (keep = person_id CumHaz_retain rename = (CumHaz_retain = Death_NelsonAalen_link));
	merge temp2 (in=ina) NelsonAalenScores3 (in=inb);
	by Death_t;
	if ina then output;
run;

/* Now merge with the main dataset that contains all individuals */
/* People not in temp3 will have a missing value */
proc sort data = temp3; by person_id;run;

data datint.varDeath_&cohort (keep = person_id Death_t Death_c Death_NelsonAalen Death_NelsonAalen_link);
	merge datint.varDeath_&cohort (in=ina) temp3 (in=inb);
	by person_id;
	if ina and not inb then Death_NelsonAalen_link = .;
	if ina then output;
run;

/* Sort by person_id again */
proc sort data=datint.varDeath_&cohort; by person_id; run;

%mend;

%NelsonAalen_link(A);
%NelsonAalen_link(B);
