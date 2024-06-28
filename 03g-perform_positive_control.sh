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
exec &> >(tee ${section_03g_logfile})
print_version

# Get the control CPG
awk -F',' -v cpg=${positive_control_cpg} '{ if(NR == 1 || $1 == cpg) { print $0 }}' ${untransformed_methylation_adjusted_pcs}.csv > ${untransformed_methylation_adjusted}.positive_control
awk -F',' -v cpg=${positive_control_cpg} '{ if(NR == 1 || $1 == cpg) { print $0 }}' ${transformed_methylation_adjusted_pcs}.csv > ${transformed_methylation_adjusted}.positive_control
#awk -v col1=IID -v col2=${positive_control_cpg} 'NR==1{for(i=1;i<=NF;i++){if($i==col1)c1=i; if ($i==col2)c2=i;} print $c1 " " $c2} NR>1{print $c1 " " $c2}' ${untransformed_methylation_adjusted}.txt > ${untransformed_methylation_adjusted}.positive_control
#awk -v col1=IID -v col2=${positive_control_cpg} 'NR==1{for(i=1;i<=NF;i++){if($i==col1)c1=i; if ($i==col2)c2=i;} print $c1 " " $c2} NR>1{print $c1 " " $c2}' ${transformed_methylation_adjusted}.txt > ${transformed_methylation_adjusted}.positive_control


nrow=`cat ${untransformed_methylation_adjusted}.positive_control | wc -l`
if [ "$nrow" -lt "2" ];
then
	echo "The positive control CPG ${positive_control_cpg} appears to be missing for mQTL analysis on untransformed methylation data. Please check."
	exit
fi

${R_directory}Rscript resources/genetics/make_control.R \
	${untransformed_methylation_adjusted}.positive_control \
	${intersect_ids_plink} \
	${untransformed_methylation_adjusted}.positive_control.plink


nrow=`cat ${untransformed_methylation_adjusted}.positive_control | wc -l`
if [ "$nrow" -lt "2" ];
then
    	echo "The positive control CPG ${positive_control_cpg} appears to be missing for mQTL analysis on transformed methylation data. Please check."
        exit
fi

${R_directory}Rscript resources/genetics/make_control.R \
	${transformed_methylation_adjusted}.positive_control \
	${intersect_ids_plink} \
	${transformed_methylation_adjusted}.positive_control.plink


echo "Perform plink"

${plink2} \
	--bfile ${bfile} \
	--pheno ${untransformed_methylation_adjusted}.positive_control.plink \
	--glm allow-no-covars \
	--out ${section_03_dir}/positive_control_untransformed_${positive_control_cpg}

tr -s " " < ${section_03_dir}/positive_control_untransformed_${positive_control_cpg}.PHENO1.glm.linear | gzip -c > ${section_03_dir}/positive_control_untransformed_${positive_control_cpg}.PHENO1.glm.linear.gz
rm ${section_03_dir}/positive_control_untransformed_${positive_control_cpg}.PHENO1.glm.linear

echo "make manhattan and qq plots"

echo ${section_03_dir}/positive_control_untransformed_${positive_control_cpg}.PHENO1.glm.linear.gz > ${section_03_dir}/positive.control.untransformed.file.txt
${R_directory}Rscript resources/genetics/plot_gwas.R \
	${section_03_dir}/positive.control.untransformed.file.txt \
	12 \
	1 \
	2 \
	3 \
	TRUE \
	${positive_control_snp_chr} \
	${positive_control_snp_pos} \
	${positive_control_snp_window} \
	${positive_control_threshold} 

echo "Perform plink"

${plink2} \
	--bfile ${bfile} \
	--pheno ${transformed_methylation_adjusted}.positive_control.plink \
	--glm allow-no-covars \
	--out ${section_03_dir}/positive_control_transformed_${positive_control_cpg}

tr -s " " < ${section_03_dir}/positive_control_transformed_${positive_control_cpg}.PHENO1.glm.linear | gzip -c > ${section_03_dir}/positive_control_transformed_${positive_control_cpg}.PHENO1.glm.linear.gz
rm ${section_03_dir}/positive_control_transformed_${positive_control_cpg}.PHENO1.glm.linear


echo "make manhattan and qq plots"
echo ${section_03_dir}/positive_control_transformed_${positive_control_cpg}.PHENO1.glm.linear.gz >${section_03_dir}/positive.control.transformed.file.txt
${R_directory}Rscript resources/genetics/plot_gwas.R \
	${section_03_dir}/positive.control.transformed.file.txt \
	12 \
	1 \
	2 \
	3 \
	TRUE \
	${positive_control_snp_chr} \
	${positive_control_snp_pos} \
	${positive_control_snp_window} \
	${positive_control_threshold}


echo "Successfully completed check"
