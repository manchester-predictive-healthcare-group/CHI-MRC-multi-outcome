/* Combine all the extracted variables into a dataset */

/* file with libname statements for this study */;
%include "/mnt/bmh01-rds/mrc-multi-outcome/code_cohort_extractions/Aurum_65plus/general_sas/libb_Aurum_65plus.sas"; 


/* Create a macro which lists all the datasets to be combined */
%macro list_datasets;
	cohort_in (in=ina) 
	/* Predictors */
	datint.varAge_&cohort
	datint.varEthnicity_&cohort  
	datint.varBMI_&cohort
	datint.varCholhdl_ratio_&cohort
	datint.varSBP_&cohort
	datint.varSmoking_&cohort
	datint.varIMD_&cohort
	datint.varAlcohol_misuse_primhist_&cohort 
		datint.varEating_disorders_primhist_&cohort
		datint.varAsthma_primhist_&cohort
		datint.varAnxiety_disorders_primhist_&cohort
		datint.varDepression_primhist_&cohort
		datint.varBronchiectasis_primhist_&cohort
		datint.varVisual_impairment_primhist_&cohort
		datint.varHepatic_failure_primhist_&cohort
		datint.varViral_hepatitis_primhist_&cohort
		datint.varSinusitis_primhist_&cohort
		datint.varCOPD_primhist_&cohort
		datint.varDementia_primhist_&cohort
		datint.varDiverticular_primhist_&cohort
		datint.varEpilepsy_primhist_&cohort
		datint.varHearing_loss_primhist_&cohort
		datint.varHypertension_primhist_&cohort
		datint.varIBS_primhist_&cohort
		datint.varIntellectual_dis_primhist_&cohort
		datint.varMS_primhist_&cohort
		datint.varParkinsons_primhist_&cohort
		datint.varPerip_vascular_primhist_&cohort
		datint.varPsoriasis_primhist_&cohort
		datint.varSubstance_misuse_primhist_&cohort
		datint.varRA_primhist_&cohort
		datint.varSchizophrenia_primhist_&cohort
		datint.varBipolar_primhist_&cohort
		datint.varThyroid_primhist_&cohort
		datint.varPeptic_ulcer_primhist_&cohort
		datint.varIBD_primhist_&cohort
		datint.varProstate_primhist_&cohort
	/* Outcomes */
	/* Death */
	datint.varDeath_&cohort
	/* CKD */
	datint.varCKD_primev_&cohort datint.varCKD_primhist_&cohort datint.varCKD_hesev_&cohort datint.varCKD_heshist_&cohort
	/* Diabetes */
	datint.varDiab_primev_&cohort datint.varDiab_primhist_&cohort datint.varDiab_t1_hesev_&cohort datint.varDiab_t1_heshist_&cohort datint.varDiab_t2_hesev_&cohort datint.varDiab_t2_heshist_&cohort
	/* AF */
	datint.varAF_primev_&cohort datint.varAF_primhist_&cohort datint.varAF_hesev_&cohort datint.varAF_heshist_&cohort
	/* HF */
	datint.varHF_primev_&cohort datint.varHF_primhist_&cohort datint.varHF_hesev_&cohort datint.varHF_heshist_&cohort
	/* CHD/MI raw */
	datint.varCHD_primev_&cohort datint.varCHD_primhist_&cohort datint.varMI_primev_&cohort datint.varMI_primhist_&cohort datint.varCHD_MI_hesev_&cohort datint.varCHD_MI_heshist_&cohort
	/* Stroke/TIA raw */
	datint.varStroke_primev_&cohort datint.varStroke_primhist_&cohort datint.varTIA_primev_&cohort datint.varTIA_primhist_&cohort datint.varStroke_TIA_hesev_&cohort datint.varStroke_TIA_heshist_&cohort;
%mend;



/* Set cohort A */
title 'cohort A';
%let cohort = A;

/* Read in cohort and retain key variables */
data cohort_in;
	set datint.cohort_base;
	gender = gender_concept_id - 1;
	keep person_id gender care_site_id region_concept_id dtcens dtcens_combdeath death_date death_date_ons death_date_comb study_dtindex_&cohort.
	linkage_elig hes_e death_e lsoa_e;
run;

/* Merge this cohort with all the variables that were derived */
data datint.cohort_var_&cohort.;
	merge %list_datasets;
	by person_id;
	if ina then output;
run;

/* Create new variables, where we combine CHD/MI and Stroke/TIA extracted from primary care into one variable */
data datint.cohort_var_&cohort.;
	set datint.cohort_var_&cohort.;
	/* Combine CHD/MI and Stroke/MI from the primary care extraction, which were done seperately */
	CHD_MI_primhist = max(CHD_primhist, MI_primhist);
	CHD_MI_primhist_t = max(CHD_primhist_t, MI_primhist_t);
	CHD_MI_primev_c = max(CHD_primev_c, MI_primev_c);
	CHD_MI_primev_t = min(CHD_primev_t, MI_primev_t);
	Stroke_TIA_primhist = max(Stroke_primhist, TIA_primhist);
	Stroke_TIA_primhist_t = max(Stroke_primhist_t, TIA_primhist_t);
	Stroke_TIA_primev_c = max(Stroke_primev_c, TIA_primev_c);
	Stroke_TIA_primev_t = min(Stroke_primev_t, TIA_primev_t);
run;

proc freq data = datint.cohort_var_&cohort.; tables gender; run;

proc contents data = datint.cohort_var_&cohort.;run;

/* Repeat process for cohort B */
title 'cohort B';
%let cohort = B;

/* Read in cohort and retain key variables */
data cohort_in;
	set datint.cohort_base;
	gender = gender_concept_id - 1;
	keep person_id gender care_site_id region_concept_id dtcens dtcens_combdeath death_date death_date_ons death_date_comb study_dtindex_&cohort.
	linkage_elig hes_e death_e lsoa_e;
run;

/* Merge this cohort with all the variables that were derived */
data datint.cohort_var_&cohort.;
	merge %list_datasets;
	by person_id;
	if ina then output;
run;

/* Create new variables, where we combine CHD/MI and Stroke/TIA extracted from primary care into one variable */
data datint.cohort_var_&cohort.;
	set datint.cohort_var_&cohort.;
	/* Combine CHD/MI and Stroke/MI from the primary care extraction, which were done seperately */
	CHD_MI_primhist = max(CHD_primhist, MI_primhist);
	CHD_MI_primhist_t = max(CHD_primhist_t, MI_primhist_t);
	CHD_MI_primev_c = max(CHD_primev_c, MI_primev_c);
	CHD_MI_primev_t = min(CHD_primev_t, MI_primev_t);
	Stroke_TIA_primhist = max(Stroke_primhist, TIA_primhist);
	Stroke_TIA_primhist_t = max(Stroke_primhist_t, TIA_primhist_t);
	Stroke_TIA_primev_c = max(Stroke_primev_c, TIA_primev_c);
	Stroke_TIA_primev_t = min(Stroke_primev_t, TIA_primev_t);
run;
proc freq data = datint.cohort_var_&cohort.; tables gender; run;

proc contents data = datint.cohort_var_&cohort.;run;


