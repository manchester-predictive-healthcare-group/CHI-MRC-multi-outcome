/* Exploratory code*/

*/ file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 

ODS PDF file = "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/p1_prelim/p1.2_format_extract_check.pdf";

/* Print the first rows of the person file, and practice files */
ODS TEXT  = "person file";
proc print data = datext.pers(obs=50);run;

ODS TEXT  = "practice file";
proc print data = datext.practice(obs=50);run;


ODS TEXT  = "medfile";
proc contents data = datext.allmedfmt_1;run;
proc print data = datext.allmedfmt_1(obs=50);run;

ODS TEXT  = "rx";
proc contents data = datext.allrxfmt_1;run;
proc print data = datext.allrxfmt_1(obs=50);run;

ODS TEXT  = "test";
proc contents data = datext.alltestfmt_1;run;
proc print data = datext.alltestfmt_1(obs=50);run;

ODS PDF close;




