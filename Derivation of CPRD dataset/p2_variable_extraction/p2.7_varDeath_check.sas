/* Check the death variables, as the main programs prints > 1GB output */


/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


proc print data = datint.varDeath_A(obs = 100); run;
proc print data = datint.varDeath_B(obs = 100); run;

/* Note that for cohort B, some of these people will have a negative death_t */
proc univariate data = datint.varDeath_A;var Death_t;run;
proc univariate data = datint.varDeath_B;var Death_t;run;

proc freq data = datint.varDeath_A; tables Death_c; run;
proc freq data = datint.varDeath_B; tables Death_c; run;


proc means data = datint.varDeath_A mean median min max n nmiss; var Death_NelsonAalen; run; 
proc means data = datint.varDeath_B mean median min max n nmiss; var Death_NelsonAalen; run; 

data temp_vardeatha;
	set datint.varDeath_A;
	proc sort; by Death_t;
run;

title 'vardeatha';
proc print data = temp_vardeatha (firstobs = 4024000);run;
title 'NelsonAalenA';
proc print data = datint.NelsonAalenScores3_A;run;
title 'NelsonAalenA link';
proc print data = datint.NelsonAalenScores3_link_A;run;

data temp_vardeathb;
	set datint.varDeath_B;
	proc sort; by Death_t;
run;

title 'vardeathb';
proc print data = temp_vardeathb (firstobs = 4024000);run;
title 'NelsonAalenB';
proc print data = datint.NelsonAalenScores3_B;run;
title 'NelsonAalenB link';
proc print data = datint.NelsonAalenScores3_link_B;run;


