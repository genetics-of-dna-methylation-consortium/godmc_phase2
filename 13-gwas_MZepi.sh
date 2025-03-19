#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_13_dir}/logs_a
touch ${section_13a_logfile}
exec &> >(tee ${section_13a_logfile})
print_version


# Step 1: compute epigenetic scores




${R_directory}Rscript resources/MZtwin/MZEpiScore.R \
    ${betas} \
    ${bfile}.fam \
    ${phenotypes_MZT}	\
    ${section_13_dir}/

echo "Computing Epi-MZ scores"

cut -d' ' -f1-12 "${pca}.eigenvec" > ${pca}_10.eigenvec

# Step 2: generate a sparse genetic relationship matrix (GRM) and PCA ###################################
# For family data, use all samples, correcting for the full (sparse) GRM.
if [ "${related}" = "yes" ]
then

${gcta} \
	--grm ${grmfile_all} \
	--make-bK-sparse 0.05 \
	--autosome \
	--make-grm \
	--out ${grmfile_fast}_gwas13 \
	--thread-num ${nthreads}


echo 'Done on making bK sparse'

# Step 3: fastGWA ###################################
${gcta} \
          --bfile ${bfile} \
	  --out ${section_13_dir}/GWASepiMZ \
          --grm-sparse ${grmfile_fast}_gwas13 \
          --fastGWA-mlm \
          --pheno ${section_13_dir}/MZEpi.pheno \
	  --qcovar ${pca}_10.eigenvec	\
	  --thread-num ${nthreads}

#For non-family data, use sparse GRM generated earlier (10a) with --grm-cutoff of 0.05
elif [ "${related}" = "no" ]
then

# Step 3: fastGWA ###################################
${gcta} \
          --bfile ${bfile} \
	  --out ${section_13_dir}/GWASepiMZ \
          --grm-sparse ${grmfile_fast}  \
          --fastGWA-mlm \
          --pheno ${section_13_dir}/MZEpi.pheno \
	  --qcovar ${pca}_10.eigenvec	\
	  --thread-num ${nthreads}

fi

# Step 4: Visualization ###################################

rm -f ${section_13_dir}/GWAlist.txt
find ${section_13_dir} -type f -name "*.fastGWA" > ${section_13_dir}/GWAlist.txt
${R_directory}Rscript resources/genetics/plot_gwas.R \
	      ${section_13_dir}/GWAlist.txt \
		10 \
		1 \
		3 \
		2 \
		TRUE \
		0 \
		0 \
		0 \
		0 

rm -f ${section_13_dir}/GWAlist.txt

echo "Successfully finished the GWAS on MZepigeneticsignature!"