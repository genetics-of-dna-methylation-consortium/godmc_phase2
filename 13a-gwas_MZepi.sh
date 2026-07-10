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
    ${section_13_dir}/ \
    ${home_directory}/processed_data/covariate_data/ \
    ${pca}_10.eigenvec \

echo "Finished computing Epi-MZ scores"

cut -d' ' -f1-12 "${pca}.eigenvec" > ${pca}_10.eigenvec

# Step 2: generate a sparse genetic relationship matrix (GRM) and PCA ###################################
# For family data, use all samples, correcting for the full (sparse) GRM.
if [ "${related}" = "yes" ]
then

echo "GWAS in related individuals is performed"

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
	  --out ${section_13_dir}/GWASepiMZ_allRelated \
          --grm-sparse ${grmfile_fast}_gwas13 \
          --fastGWA-mlm \
          --pheno ${section_13_dir}/MZEpi_all.pheno \
	  --qcovar ${home_directory}/processed_data/covariate_data/covariates_intersectids.numeric	\
	  --thread-num ${nthreads}



if [ -f ${section_13_dir}/MZEpi_MZtwins.pheno ]
then

echo "MZEpi_MZtwins.pheno is present, performing GWAS in MZ twins"

${gcta} \
          --bfile ${bfile} \
	  --out ${section_13_dir}/GWASepiMZ_MZtwins \
          --grm-sparse ${grmfile_fast}_gwas13 \
          --fastGWA-mlm \
          --pheno ${section_13_dir}/MZEpi_MZtwins.pheno \
	  --qcovar ${home_directory}/processed_data/covariate_data/covariates_intersectids.numeric	\
	  --thread-num ${nthreads}



if [ -f ${section_13_dir}/MZEpi_nontwins.pheno ]
then

echo "MZEpi_nontwins.pheno is present, performing GWAS in non-twins"


${gcta} \
 --bfile ${bfile} \
	  --out ${section_13_dir}/GWASepiMZ_nontwinsrelated \
          --grm-sparse ${grmfile_fast}_gwas13 \
          --fastGWA-mlm \
          --pheno ${section_13_dir}/MZEpi_nontwins.pheno \
	  --qcovar ${home_directory}/processed_data/covariate_data/covariates_intersectids.numeric	\
	  --thread-num ${nthreads}


fi
fi


#For non-family data, use sparse GRM generated earlier (10a) with --grm-cutoff of 0.05
elif [ "${related}" = "no" ]
then

echo "GWAS in unrelated non-twins is performed"

# Step 3: fastGWA ###################################
${gcta} \
          --bfile ${bfile} \
	  --out ${section_13_dir}/GWASepiMZ_nontwinsunrelated \
          --grm-sparse ${grmfile_fast}  \
          --fastGWA-mlm \
          --pheno ${section_13_dir}/MZEpi_nontwins.pheno \
	  --qcovar ${home_directory}/processed_data/covariate_data/covariates_intersectids.numeric	\
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
		0 \

rm -f ${section_13_dir}/GWAlist.txt

echo "Successfully finished the GWAS on MZepigeneticsignature!"