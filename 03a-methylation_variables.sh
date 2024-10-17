#!/bin/bash

source resources/setup.sh "$@"
exec &> >(tee ${section_03a_logfile})
print_version

# Remove outliers
echo "Removing outliers"
${R_directory}Rscript resources/methylation/remove_outliers.R \
		${betas} \
		${methylation_no_outliers} \
		${cohort_descriptives_commonids} \
		${methylation_summary} \
		${intersect_ids} \
		${covariates} \
		${covariates_intersect} \
		${bfile}.fam \
		${bfile}.bim

# Predict smoking
echo "Estimating smoking"
${R_directory}Rscript resources/smoking/smoking_predictor.R \
		${methylation_no_outliers} \
		${bfile}.fam \
		${smoking_pred} \
		${smoking_pred_plot} \
		${smoking_pred_SD} \
		${covariates} 

#Predict cell counts - force for all datasets to ensure consistency
echo "Predicting cell counts using epiDISH"
${R_directory}Rscript resources/cellcounts/cellcounts_epiDISH.R \
        ${methylation_no_outliers} \
        ${cellcounts_cov} \
        ${cellcounts_plot} \
        ${cellcounts_summary}


if [ "${measured_cellcounts}" != "NULL" ]
then
  echo "Comparing measured with predicted cellcounts"
    ${R_directory}Rscript resources/cellcounts/correlation.R \
	          ${cellcounts_cov} \
 	          ${measured_cellcounts} \
	          ${cor_matrix} \
	          ${cor_plot}

fi


# Estimate age accelerated residuals
echo "Estimating age"
${R_directory}Rscript resources/dnamage/dnamage.R \
		${methylation_no_outliers} \
		${covariates_intersect} \
		${bfile}.fam \
		${age_pred} \
		${age_pred_plot} \
		${age_pred_SD} \
		${age_pred_sumstats} \
		${smoking_pred}.txt \
		${cellcounts_cov} 


# Organise covariates
echo "Organising covariates"
${R_directory}Rscript resources/genetics/covariates.R \
	${covariates_intersect} \
	${pcs_all} \
	${cellcounts_cov} \
	${smoking_pred}.txt \
	${bfile}.fam \
	${covariates_combined}

echo "Successfully created methylation-related variables"
