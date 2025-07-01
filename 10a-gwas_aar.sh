#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_10_dir}/logs_a
touch ${section_10a_logfile}
exec &> >(tee ${section_10a_logfile})
print_version


# Step 0: re-run the epigenetic clocks for age###################################
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


# Step 1: run the GCTA-GRM: calculating the genetic relationship matrix from autosomal SNPs

# ${gcta} \
# 	--bfile ${bfile} \
# 	--autosome \
# 	--maf 0.05 \
# 	--make-grm \
# 	--out  \
# 	--thread-num 10


# Step 2: generate a sparse genetic relationship matrix (GRM) and PCA ###################################

${gcta} \
	--grm ${grmfile_all} \
	--grm-cutoff 0.05 \
	--make-grm \
	--out ${grmfile_all}_gaws10 \
	--thread-num ${nthreads}

${gcta} \
	--grm ${grmfile_all}_gaws10 \
	--pca 20 \
	--out ${home_directory}/processed_data/genetic_data/gaws10_pc

${gcta} \
	--grm ${grmfile_all}_gaws10 \
	--make-bK-sparse 0.05 \
	--autosome \
	--make-grm \
	--out ${grmfile_fast} \
	--thread-num ${nthreads}


echo 'Done on making bK sparse'

# Step 3: fastGWA ###################################
i=1
tail -n +2 ${age_pred}.txt > ${age_pred}.plink
clock_names=$(cut -d" " -f 3- ${age_pred}.txt | head -n 1)

for clock_name in $clock_names
do
    ${gcta} \
          --bfile ${bfile}  \
          --grm-sparse ${grmfile_fast} \
          --fastGWA-mlm \
          --mpheno $i \
          --pheno ${age_pred}.plink \
          --h2-limit 100 \
          --out ${section_10_dir}/${clock_name}

	
	${gcta} \
          --bfile ${bfile}  \
          --grm-sparse ${grmfile_fast}  \
          --fastGWA-mlm \
          --mpheno $i \
		  --h2-limit 100 \
		  --qcovar ${home_directory}/processed_data/genetic_data/gaws10_pc.eigenvec \
          --pheno ${age_pred}.plink \
		  --thread-num ${nthreads} \
          --out ${section_10_dir}/${clock_name}_PCA

  i=$(($i+1))
  echo "Done the GWAS on" $clock_name 
done


# Step 4: Visulization ###################################

rm -f ${section_10_dir}/GWAlist.txt
find ${section_10_dir} -type f -name "*.fastGWA" > ${section_10_dir}/GWAlist.txt
${R_directory}Rscript resources/genetics/plot_gwas.R \
	      ${section_10_dir}/GWAlist.txt \
	      10 \
	      1 \
	      3 \
	      2 \
	      TRUE \
	      0 \
	      0 \
	      0 \
	      0 

rm -f ${section_10_dir}/GWAlist.txt

echo "Successfully finished the GWAS on age accelerations!"


