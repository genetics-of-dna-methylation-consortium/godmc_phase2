#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_14_logfile})
print_version

echo "nc886 and clustering"
${R_directory}Rscript resources/methylation/nc886_clustering.R \
    ${betas} \
    ${section_14_dir}/

if [ -f "${section_14_dir}/nc886_groups.txt" ]; then

echo "GWAS intermediated vs non-methylated"
${plink2} --bfile ${bfile} --pheno ${section_14_dir}/nc886_groups.txt --pheno-name intermediate_non --maf 0.01 --glm allow-no-covars --ci 0.95 --hardy --freq --geno-counts --no-fid --out ${section_14_dir}/intermediate_non


echo "GWAS intermediated vs imprinted"
${plink2} --bfile ${bfile} --pheno ${section_14_dir}/nc886_groups.txt --pheno-name intermediate_imprinted --maf 0.01 --glm allow-no-covars --ci 0.95 --hardy --freq --geno-counts --no-fid --out ${section_14_dir}/intermediate_imprinted


echo "GWAS non-methylated vs imprinted"
${plink2} --bfile ${bfile} --pheno ${section_14_dir}/nc886_groups.txt --pheno-name non_imprinted --maf 0.01 --glm allow-no-covars --ci 0.95 --hardy --freq --geno-counts --no-fid --out ${section_14_dir}/non_imprinted


else
	echo "Problem: There are not enough intermediatelly methylated individuals to continue the analysis"
		exit 1
  
fi

echo "Successfully completed nc866 analysis"
