#!/bin/bash

set -e

# Initialize variables
config_file="./config"

# Parse options using getopts
while getopts "c:" opt; do
    case $opt in
        c) config_file=$OPTARG ;;
        *) echo "Usage: $0 -c <config_file>"
           exit 1 ;;
    esac
done

# Shift option arguments, so $1 becomes the first positional argument
shift $((OPTIND - 1))

set -e
echo "-----------------------------------------------"
echo ""
echo "Using config located at:" ${config_file}
echo ""
echo "-----------------------------------------------"
	
source ${config_file}
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

#Predict cell counts
echo "Processing cell counts"
if [ "${cellcounts_required}" = "yes" ]
then

echo "Predicting cell counts using epiDISH"
${R_directory}Rscript resources/cellcounts/cellcounts_epiDISH.R \
        ${methylation_no_outliers} \
        ${cellcounts_cov} \
        ${cellcounts_plot} \
        ${cellcounts_summary}

fi

if [ "${measured_cellcounts}" != "NULL" ]
then
  echo "Comparing measured with predicted cellcounts"
    ${R_directory}Rscript resources/cellcounts/correlation.R \
	          ${cellcounts_cov} \
 	          ${measured_cellcounts} \
	          ${cor_matrix} \
	          ${cor_plot}

fi

if [ "${cellcounts_required}" = "no" ]
then
  echo "No cell counts required"
	  cellcounts_cov="NULL"

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
