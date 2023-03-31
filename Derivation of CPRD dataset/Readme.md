This code is provided for transparency but has not been written in a reproducible format (i.e. it could not easily be ran on your system). Initial extraction of the cohort and baseline/outcome data for individuals in this cohort is predominately done in SAS, with the exception of the algorithm to identify CKD episodes from eGFR scores, which is done in R. Imputation is done in R.

p1 files derive the cohort and apply inclusion/exclusion criteria.

p2 files extract the desired variables for individuals in the cohort derived from p1. 

The p2.1 - 2.3 files create a single cohort file AFTER all the variables have been extracted. This therefore means they should be ran last. 

The p2.5 files extract baseline data which is not a 'history of' variables.

The p2.6 file extracts baseline data for all variables which are of type 'history of'.

The p2.7 files extract outcome data (time until event, and a censoring indicator).

p3 files import the cohort into R and run multiple imputation.