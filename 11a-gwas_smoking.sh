#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_11_dir}/logs_a
touch ${section_11a_logfile}
exec &> >(tee ${section_11a_logfile})
print_version

# Step 0: epi Smoking score prediction ###################################
${R_directory}Rscript resources/smoking/smoking_predictor.R \
		${methylation_no_outliers} \
		${bfile}.fam \
		${smoking_pred} \
		${smoking_pred_plot} \
		${smoking_pred_SD} \
		${covariates} 



# Step 1: Check if cellcount is generated ###################################
if [ "${sorted_methylation}" = "yes" ]
then
  echo "No cell counts required" 
  cellcounts_cov="NULL"
fi


# Step 2: Generate qcovar for GWAS of smoking ###################################
${R_directory}Rscript resources/smoking/qcovar_gwas_smoking.R \
	      ${cellcounts_cov} \
	      ${covariates} \
          ${bfile}.fam \
	      ${smoking_pred} \
		  ${home_directory}/results/11/smoking_stats


# Step 3:GWAS of smoking ###################################
echo "Running GWAS for smoking."
${gcta} \
  		--bfile ${bfile} \
  		--grm-sparse ${grmfile_fast}  \
  		--fastGWA-mlm  \
  		--pheno ${smoking_pred}.smok.plink  \
  		--autosome \
      	--h2-limit 100 \
		--thread-num ${nthreads} \
  		--out ${section_11_dir}/gwas_smoking


${gcta} \
  		--bfile ${bfile} \
  		--grm-sparse ${grmfile_fast}  \
        --fastGWA-mlm \
  		--pheno ${smoking_pred}.smok.plink \
		--qcovar ${home_directory}/processed_data/genetic_data/gaws10_pc.eigenvec \
		--autosome \
		--h2-limit 100 \
		--thread-num ${nthreads} \
  		--out ${section_11_dir}/gwas_smoking_PCA



# Step 4: Visulization ###################################
echo "Making plots"
rm -f ${section_11_dir}/GWAlist.txt
find ${section_11_dir} -type f -name "*.fastGWA" > ${section_11_dir}/GWAlist.txt
${R_directory}Rscript resources/genetics/plot_gwas.R \
		${section_11_dir}/GWAlist.txt \
		10 \
		1 \
		3 \
		2 \
		TRUE \
		0 \
		0 \
		0 \
		0 

rm -f ${section_11_dir}/GWAlist.txt

echo "Successfully finished the GWAS on smoking!"
