#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_06_dir}
mkdir -p ${section_06_dir}/logs_a
exec &> >(tee ${section_06a_logfile})
print_version

echo "Started to run 06a at $(date)"
echo "STEP1: downloading and checking if all required files are available"

sftp -P 2222 -oIdentityFile=$key ${sftp_username}@${sftp_address}:/sftp/vmeQTL_resources <<EOF
get *
EOF

mv pruned_snps_vmeQTL_phase2_chr* ${scripts_directory}/resources/genetics
mv vmeQTL* ${scripts_directory}/resources/methylation/vmeQTL

# add one step to download vmeQTL_list
if [ -f ${covariates_combined}.txt ];
then
    echo "covariate file exist"
else
    echo "candidate file DOES NOT exist. Please check if 03a has been run successfully."
    exit
fi

if [ -f ${BMI} ];
then
    echo "BMI file exist"
    n=`grep -v IID ${BMI} | wc -l`
    echo "$n samples with BMI information"
else
    echo "BMI file DOES NOT exist. Please include it following Wiki guidance."
    exit
fi

echo "STEP2: checking and processing genotype data"

count=`ls ${tabfile}.tab.*bed | wc -l`
if [ ${count} -eq ${genetic_chunks} ];
then
    echo "Genetic tab files has been generated"
else
    for chunk in $(seq 1 ${genetic_chunks})
    do

    startsnp=`cut -f 1 ${tabfile}.tab.${chunk} | grep -v snpid | head -n1`
    endsnp=`cut -f 1 ${tabfile}.tab.${chunk} | tail -n1`

    ${plink} \
        --bfile ${bfile} \
        --snps ${startsnp}-${endsnp} \
        --make-bed \
        --out ${tabfile}.tab.${chunk}
    done
fi

participate07=`grep ${study_name} ${scripts_directory}/resources/methylation/vmeQTL/vmeQTL_phase1_cohort_list.txt | wc -l`

if [ ${participate07} -eq 1 ];
then
    echo "You have participated into module 07. Start to generate the genotype data to run missing CpGs in module 07"
    for c in $(seq 1 22);
    do
    ${plink} \
        --bfile ${bfile} \
        --chr ${c} \
        --make-bed \
        --out ${tabfile}.chr${c}
    done
fi

echo "Generate genotype data with pruned SNPs"
for c in $(seq 1 22);
do
    lines=$(( ($(wc -l < ${scripts_directory}/resources/genetics/pruned_snps_vmeQTL_phase2_chr${c}) + ${prune_sub} - 1) / ${prune_sub} ))
    split -l "$lines" ${scripts_directory}/resources/genetics/pruned_snps_vmeQTL_phase2_chr${c} _tmp_subfile_ && \
    i=1 && for f in _tmp_subfile_*; do mv "$f" "${scripts_directory}/resources/genetics/pruned_snps_vmeQTL_phase2_chr${c}_sub$i"; ((i++)); done
    ${plink} \
        --bfile ${bfile} \
        --extract ${scripts_directory}/resources/genetics/pruned_snps_vmeQTL_phase2_chr${c} \
        --make-bed \
        --out ${tabfile}.prunedSNPs.chr${c}

    for chunk in $(seq 1 ${prune_sub});
	do
	${plink} \
        --bfile ${tabfile}.prunedSNPs.chr${c} \
        --extract ${scripts_directory}/resources/genetics/pruned_snps_vmeQTL_phase2_chr${c}_sub${chunk} \
        --make-bed \
        --out ${tabfile}.prunedSNPs.chr${c}.chunk${chunk}
	done
done

echo "Generate the genotype data of previously identified vQTLs"
${plink} \
    --bfile ${bfile} \
    --extract ${vmeQTL_list3} \
    --make-bed \
    --out ${tabfile}.vQTLs

echo "Generate epistasis input data"
row_num=1
IFS=$'\n'
for row in $(cat ${vmeQTL_list4}); 
do
    IFS=$'\t' read -r snp1 snp2 <<< "$row"
    
    echo "Processing Row ${row_num}: $snp1 & $snp2"
    count1=`grep -w $snp1 ${bfile}.bim | wc -l`
    count2=`grep -w $snp2 ${bfile}.bim | wc -l`
    if [ "$count1" -eq 1 ] && [ "$count2" -eq 1 ]
    then
        ${plink} \
            --bfile ${bfile} \
            --snp ${snp1} \
            --make-bed \
            --out ${tabfile}_epi_row${row_num}

        ${plink} \
            --bfile ${bfile} \
            --snp ${snp2} \
            --recode A \
            --out ${tabfile}_epi_row${row_num}
    else
        echo "Skipping row ${row_num}: One or both SNPs missing"
    fi
    ((row_num++))
done

echo "STEP3: generating environmental factor file"
${R_directory}Rscript ${scripts_directory}/resources/methylation/generate_environment_file.R \
    ${covariates_combined}.txt \
    ${BMI} \
    ${envs_input} \
    ${section_06_dir}/E_plots.pdf \
    ${section_06_dir}/E_summary.csv


echo "STEP4: generating and processing methylation data"
mkdir -p ${meth_vmeQTL_directory}/vmeQTL_phase2/

if [ ${participate07} -eq 1 ];
then
    echo "You have participated module 07. Generating subset of methylation files to run missing CpGs"
    for c in $(seq 1 22);
    do
    ${R_directory}Rscript ${scripts_directory}/resources/methylation/observe_missing_cpgs.R \
        ${c} \
        ${vmeQTL_list2} \
        ${meth_vmeQTL_input_chr}${c} \
        ${meth_vmeQTL_directory}/vmeQTL_phase2/missing_cpgs_chr${c}

    ${osca} \
        --tefile ${meth_vmeQTL_directory}/vmeQTL_phase2/missing_cpgs_chr${c} \
        --methylation-m \
        --make-bod \
        --no-fid \
        --out ${meth_vmeQTL_directory}/vmeQTL_phase2/missing_cpgs_chr${c}

    ${osca} \
        --befile ${meth_vmeQTL_directory}/vmeQTL_phase2/missing_cpgs_chr${c} \
        --update-opi ${meth_vmeQTL_annotation}.opi
    
    cp ${meth_vmeQTL_input_chr}${c}.oii ${meth_vmeQTL_directory}/vmeQTL_phase2/missing_cpgs_chr${c}.oii
    
    ${R_directory}Rscript ${scripts_directory}/resources/methylation/match_oii_plink.R \
        ${meth_vmeQTL_directory}/vmeQTL_phase2/missing_cpgs_chr${c}.oii \
        ${bfile}.fam
    done
fi

count=`ls ${meth_vmeQTL_input_chr}*[0-9].bod | wc -l`
if [ $count -eq 22 ];
then
    echo "methylation data for 22 chromosomes in OSCA format exist"
else
    ${R_directory}Rscript \
        ${scripts_directory}/resources/methylation/vmeQTL_process_tabfile.R \
        ${untransformed_methylation_adjusted_pcs}.RData \
        ${meth_vmeQTL_input_chr}

    for chr in $(seq 1 22)
    do
        echo "convert chr ${chr} methylation data to bod format"
        ${osca} \
            --tefile ${meth_vmeQTL_input_chr}${chr} \
            --methylation-m \
            --make-bod \
            --no-fid \
            --out ${meth_vmeQTL_input_chr}${chr}

        ${osca} \
            --befile ${meth_vmeQTL_input_chr}${chr} \
            --update-opi ${meth_vmeQTL_annotation}.opi
        
        ${R_directory}Rscript ${scripts_directory}/resources/methylation/match_oii_plink.R \
            ${meth_vmeQTL_input_chr}${chr}.oii \
            ${bfile}.fam
done
fi

echo "Generating methylation data of the CpGs of interest"
${R_directory}Rscript ${scripts_directory}/resources/methylation/observe_subset_cpgs.R \
    ${untransformed_methylation_adjusted_pcs}.RData \
    ${vmeQTL_list5} \
    ${meth_vmeQTL_directory}/vmeQTL_phase2/cpgs_of_interest

${osca} \
    --tefile ${meth_vmeQTL_directory}/vmeQTL_phase2/cpgs_of_interest \
    --methylation-m \
    --make-bod \
    --no-fid \
    --out ${meth_vmeQTL_directory}/vmeQTL_phase2/cpgs_of_interest

${osca} \
    --befile ${meth_vmeQTL_directory}/vmeQTL_phase2/cpgs_of_interest \
    --update-opi ${meth_vmeQTL_annotation}.opi

${R_directory}Rscript ${scripts_directory}/resources/methylation/match_oii_plink.R \
    ${meth_vmeQTL_directory}/vmeQTL_phase2/cpgs_of_interest.oii \
    ${bfile}.fam


echo "Generating methylation data for interaction analysis"
${R_directory}Rscript ${scripts_directory}/resources/methylation/vmeQTL_process_tabfile_phase2_allCpGs.R \
    ${untransformed_methylation_adjusted_pcs}.RData \
    ${meth_chunks} \
    ${meth_vmeQTL_directory}/vmeQTL_phase2/adjustcovs_cpg_phase2_allCpGs_chunk

${R_directory}Rscript ${scripts_directory}/resources/methylation/vmeQTL_process_tabfile_phase2_cisCpGs.R \
    ${untransformed_methylation_adjusted_pcs}.RData \
    ${meth_vmeQTL_directory}/vmeQTL_phase2/adjustcovs_cpg_phase2_cisCpGs_chr \
    ${vmeQTL_list1} \
    TRUE

${R_directory}Rscript ${scripts_directory}/resources/methylation/vmeQTL_process_tabfile_phase2_cisCpGs.R \
    ${untransformed_methylation_adjusted_pcs}.RData \
    ${meth_vmeQTL_directory}/vmeQTL_phase2/adjustcovs_cpg_phase2_cpg_of_interest \
    ${vmeQTL_list5} \
    FALSE

echo "06a has been done successfully at $(date)"
