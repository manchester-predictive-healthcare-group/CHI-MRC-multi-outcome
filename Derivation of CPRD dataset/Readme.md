This code is provided for transparency but has not been written in a reproducible format (i.e. it could not easily be ran on your system). Initial extraction of the cohort and baseline/outcome data for individuals in this cohort is predominately done in SAS, with the exception of the algorithm to identify CKD episodes from eGFR scores, which is done in R. Imputation is done in R.

p1 files derive the cohort and apply inclusion/exclusion criteria.

p2 files extract the desired variables for individuals in the cohort derived from p1. 

The p2.1 - 2.3 files create a single cohort file AFTER all the variables have been extracted. This therefore means they should be ran last. 

The p2.5 files extract baseline data which is not a 'history of' variables.

The p2.6 file extracts baseline data for all variables which are of type 'history of'.

The p2.7 files extract outcome data (time until event, and a censoring indicator).

p3 files import the cohort into R and run multiple imputation.

More details are provided on what each of the programs do in the Data Extraction document. This document was produced for personal use and therefore may be difficult to follow, however is being provided as may provide additional information.

Codelists have been provided in the 'codelists' folder. These codelists have not been put together by our research team. A fully referenced guide to where all of these codelists have been taken from (all publically available online) can be found in the word document: Data extraction document.docx