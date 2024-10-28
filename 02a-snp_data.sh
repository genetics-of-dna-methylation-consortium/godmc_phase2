#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_02a_logfile})
print_version

# Provision for having dosage data and convert to best guess if the analyst doesn't have this

# Copy raw data to modifiable data


# rewrite this to use keep ids in common with meth
echo "Copying genetic data to processing folder"
# cp ${bfile_raw}.bed ${bfile}.bed
# cp ${bfile_raw}.bim ${bfile}.bim
# cp ${bfile_raw}.fam ${bfile}.fam


# For datasets using 1000G exclude indels
if [ $reference == "1000G" ] ;
then 

	echo "As 1000G reference panel used, filtering to just SNP variants"
	# create list of snps
	awk '{
		if (length($5) == 1 && length($6) == 1) {  # Check if columns 5 and 6 have single characters
			if ($5 ~ /^[ATCG]/ && $6 ~ /^[ATCG]$/) {  # Check if both columns contain A, T, C, or G (removed unnecessary $ from end of regex)
				print $2;
			} 
		}
	}' "${bfile_raw}.bim" > "${bfile}.snps.keep"

		${plink2} \
			--bfile ${bfile_raw} \
			--extract ${bfile}.snps.keep \
			--keep ${intersect_ids_plink} \
			--maf ${snp_maf} \
			--hwe ${snp_hwe} \
			--geno ${snp_miss} \
			--mind ${snp_imiss} \
			--make-bed \
			--out ${bfile} \
			--allow-extra-chr --chr-set 23 \
			--chr 1-23 \
			--threads ${nthreads}
else

	${plink2} \
		--bfile ${bfile_raw} \
		--keep ${intersect_ids_plink} \
		--maf ${snp_maf} \
		--hwe ${snp_hwe} \
		--geno ${snp_miss} \
		--mind ${snp_imiss} \
		--make-bed \
		--out ${bfile} \
		--allow-extra-chr --chr-set 23 \
		--chr 1-23 \
		--threads ${nthreads}
fi




# Sex check -note this is not implemented in PLINK2 so we are going to use PLINK1.9

n23=`grep ^23 ${bfile}.bim | wc -l`

if [ "$n23" -gt "0" ]
then
	
	${plink} \
		--bfile ${bfile} \
		--split-x b37 no-fail \
		--make-bed \
		--out ${bfile}

	${plink} \
		--bfile ${bfile} \
		--check-sex \
		--out ${section_02_dir}/data
	
	nprob=`grep "PROBLEM" ${section_02_dir}/data.sexcheck |wc -l`

    if [ "$nprob" -eq "0" ]
	then
		echo "There are ${nprob} individuals that failed the sex check."
	fi


	if [ "$nprob" -gt "0" ]
	then
		echo "There are ${nprob} individuals that failed the sex check."
		echo "They will be removed."
		echo "The summary is located here:"
		echo "${section_02_dir}/data.sexcheck"

		grep "PROBLEM" ${section_02_dir}/data.sexcheck | awk '{print $1, $2}' > ${bfile}.failed_sexcheck

		${plink} \
			--bfile ${bfile} \
			--remove ${bfile}.failed_sexcheck \
			--make-bed \
			--out ${bfile}1

		mv ${bfile}1.bed ${bfile}.bed
		mv ${bfile}1.bim ${bfile}.bim
		mv ${bfile}1.fam ${bfile}.fam	

	fi

fi

# Change SNP ids to chr:position_A1_A2 (ascii sorted order)
echo "Updating SNP ID coding"
cp ${bfile}.bim ${bfile}.bim.original

${plink2} \
	--bfile ${bfile} \
	--set-all-var-ids @:#_\$1_\$2 \
	--make-bed \
	--out ${bfile}1 \
	--threads ${nthreads}

#Recode alleles to uniform format eg. remove SNPs with alternate coding
cp ${bfile}.bim ${bfile}.bim.original1
mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam

touch ${SNPfail1}

${R_directory}Rscript resources/genetics/harmonization.R \
	${bfile}.bim \
	${SNPfail1}

# Checking for any duplicate SNPs
# 'exclude-mismatch': When unequal duplicate-ID variants are found, exclude every member of the group.

${plink2} \
	--bfile ${bfile} \
	--rm-dup exclude-mismatch \
	--make-bed \
	--out ${bfile}1 \
	--threads ${nthreads}

cp ${bfile}.bim ${bfile}.bim.original2
mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam

# Remove SNPs with low info scores
awk '$3 < 0.80 {print $1}' <${quality_scores} > ${bfile}.lowinfoSNPs.txt

cat ${SNPfail1} ${bfile}.lowinfoSNPs.txt |sort -u >${bfile}.failed.SNPs.txt

n_failedSNPs=`wc -l ${bfile}.failed.SNPs.txt | awk '{ print $1 }'`

# Remove SNPs from data
echo "Removing ${n_failedSNPs} SNPs from data"

${plink2} \
	--bfile ${bfile} \
	--exclude ${bfile}.failed.SNPs.txt \
	--make-bed \
	--out ${bfile}1 \
	--threads ${nthreads}

cp ${bfile}.bim ${bfile}.bim.original3
mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam

# Make GRMs
echo "Creating kinship matrix"
gunzip -c ${hm3_snps} > temp_hm3snps.txt

#${plink2} \
#	--bfile ${bfile}3 \
#	--extract temp_hm3snps.txt \
#	--maf ${grm_maf_cutoff} \
#	--make-rel bin4 \
#	--out ${grmfile_all}3 \
#	--threads ${nthreads} \
#	--autosome

${plink2} \
	--bfile ${bfile} \
	--extract temp_hm3snps.txt \
	--maf ${grm_maf_cutoff} \
	--make-grm-bin \
	--out ${grmfile_all} \
	--threads ${nthreads} \
	--autosome


#compare to GoDMC1
#${plink} \
#	--bfile ${bfile}3 \
#	--extract temp_hm3snps.txt \
#	--maf ${grm_maf_cutoff} \
#	--make-grm-gz no-gz \
#	--out ${grmfile_all} \
#	--threads ${nthreads} \
#	--autosome

#rm temp_hm3snps.txt

#${plink2} \
#	--bfile ${bfile}2 \
#	--extract temp_hm3snps.txt \
#	--maf ${grm_maf_cutoff} \
#	--make-king triangle bin \
#	--out ${grmfile_all} \
#	--threads ${nthreads} \
#	--autosome
#rm temp_hm3snps.txt

# Create pedigree matrix if family data, otherwise remove related individuals from existing kinship and data file
if [ "${related}" = "yes" ]
then
	echo "Creating pedigree GRM"
	${R_directory}Rscript resources/relateds/grm_relateds.R ${grmfile_all} ${grmfile_relateds} ${rel_cutoff}
elif [ "${related}" = "no" ]
then
	echo "Removing any cryptic relateds"
	#${plink} \
	#	--grm-bin ${grmfile_all}3 \
	#	--rel-cutoff ${rel_cutoff} \
	#	--make-grm-bin \
	#	--out ${grmfile_unrelateds}3 \
	#	--threads ${nthreads}


#${plink2}  \
#		--bfile ${bfile}2 \
#		--king-cutoff ${bfile} 0.0884 \
#		--make-bed \
#		--out ${bfile}3 \
#		--threads ${nthreads}

${gcta} \
		--grm ${grmfile_all} \
		--grm-cutoff ${rel_cutoff} \
		--make-grm-bin \
		--out ${grmfile_all}1

mv ${grmfile_all}1.grm.N.bin ${grmfile_all}.grm.N.bin
mv ${grmfile_all}1.grm.id ${grmfile_all}.grm.id
mv ${grmfile_all}1.grm.bin ${grmfile_all}.grm.bin

${plink2} \
		--bfile ${bfile} \
		--keep ${grmfile_all}.grm.id \
		--make-bed \
		--out ${bfile}1 \
		--threads ${nthreads}

mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam

else 
	echo "Error: Set related flag in config to yes or no"
	exit 1
fi

#Calculate PCs
#gunzip -c ${hm3_snps_no_ld} > temp_hm3snpsnold.txt

${plink2} \
	--bfile ${bfile} \
	--extract temp_hm3snps.txt \
	--indep-pairwise 10000 5 0.1 \
	--maf 0.2 \
	--out ${pca} \
	--autosome \
	--threads ${nthreads}

if [ "${related}" = "no" ]
then
	${plink2} \
		--bfile ${bfile} \
		--extract ${pca}.prune.in \
		--pca 20 \
		--out ${pca} \
		--threads ${nthreads}
else

	${plink2} \
		--bfile ${bfile} \
		--extract ${pca}.prune.in \
		--make-bed \
		--out ${bfile}_ldpruned \
		--threads ${nthreads}

	${R_directory}Rscript resources/genetics/pcs_relateds.R \
		${bfile}_ldpruned \
		${pca} \
		${n_pcs} \
		${nthreads}
fi


# Get genetic outliers
echo "Detecting genetic outliers"

${R_directory}Rscript resources/genetics/genetic_outliers.R \
	${pcs_all} \
	${pca_sd} \
	${n_pcs} \
	${genetic_outlier_ids} \
	${pcaplot}



# If there are any genetic outliers then remove them and recalculate PCs
# Otherwise don't do anything

n_outliers=`wc -l ${genetic_outlier_ids} | awk '{ print $1 }'`

if [ "${n_outliers}" = "0" ]
then
	echo "No genetic outliers detected"
else
	# Remove genetic outliers from data
	echo "Removing ${n_outliers} genetic outliers from data"
	${plink2} \
		--bfile ${bfile} \
		--remove ${genetic_outlier_ids} \
		--make-bed \
		--out ${bfile}1 \
		--threads ${nthreads}
	
	mv ${bfile}1.bed ${bfile}.bed
	mv ${bfile}1.bim ${bfile}.bim
	mv ${bfile}1.fam ${bfile}.fam

	${gcta} \
		--grm ${grmfile_all} \
		--remove ${genetic_outlier_ids} \
		--make-grm-bin \
		--out ${grmfile_all}1 \
		--thread-num ${nthreads}

mv ${grmfile_all}1.grm.N.bin ${grmfile_all}.grm.N.bin
mv ${grmfile_all}1.grm.id ${grmfile_all}.grm.id
mv ${grmfile_all}1.grm.bin ${grmfile_all}.grm.bin

fi


# Find mismatched SNPs and misaligned SNPs with EasyQC
# Get frequencis for strand check
${plink2} \
	--bfile ${bfile} \
	--freq \
	--out ${bfile}

if [ -f "processed_data/genetic_data/easyQC_hrc.ecf.out" ]
then
	echo "easyqc files present from previous run which will be removed"
	rm processed_data/genetic_data/easy*
	rm processed_data/genetic_data/*easy*
else
	echo "passed file check"
fi

    
#${R_directory}Rscript ./resources/genetics/easyQC.R ${bfile}.bim ${bfile}.frq ${easyQC} ${easyQCfile} ${easyQCscript}

#cp ${scripts_directory}/resources/genetics/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz ${home_directory}/processed_data/genetic_data/
zcat ${scripts_directory}/resources/genetics/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz |perl -pe 's/^23/X/g'|gzip >${home_directory}/processed_data/genetic_data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz

replacement_text="DEFINE --pathOut "${home_directory}"/processed_data/genetic_data"
awk 'NR==3 {$0 = "'"$replacement_text"'"} 1' ${easyQCscript} > ${easyQCscript%.ecf}_edit.ecf
${R_directory}Rscript ./resources/genetics/easyQC.R ${bfile}.afreq ${easyQC} ${easyQCfile} ${easyQCscript%.ecf}_edit.ecf

rm ${home_directory}/processed_data/genetic_data/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz

mv ${home_directory}/processed_data/genetic_data/easyQC_hrc_edit.multi.AFCHECK.png ${home_directory}/results/02/easyQC_hrc.multi.AFCHECK.png
mv ${home_directory}/processed_data/genetic_data/easyQC_hrc_edit.rep ${home_directory}/results/02/easyQC_hrc.rep

# Remove mismatched SNPs and flip misaligned SNPs
echo "Remove mismatched SNPs and NO FLIPPING"

${plink2} \
	--bfile ${bfile} \
	--exclude ${easyQC}.mismatch_afcheck.failed.SNPs.txt \
	--make-bed \
	--out ${bfile}1 \
	--threads ${nthreads}

mv ${bfile}1.bed ${bfile}.bed
mv ${bfile}1.bim ${bfile}.bim
mv ${bfile}1.fam ${bfile}.fam

# From here on, we have clean data


if [ ! "${n_outliers}" = "0" ]
then

	echo "Recalculating PCs with outliers removed"

	if [ "${related}" = "no" ]
	then
		${plink2} \
			--bfile ${bfile} \
			--extract ${pca}.prune.in \
			--pca 20 \
			--out ${pca} \
			--autosome \
			--threads ${nthreads}
	else

		${plink2} \
			--bfile ${bfile} \
			--extract ${pca}.prune.in \
			--make-bed \
			--out ${bfile}_ldpruned \
			--autosome \
			--threads ${nthreads}

		${R_directory}Rscript resources/genetics/pcs_relateds.R \
			${bfile}_ldpruned \
			${pca} \
			${n_pcs} \
			${nthreads}
	fi

fi

# Get frequencies, missingness, hwe, info scores
plink_files=${section_02_dir}/data*.gz
if [[ "${#plink_files[@]}" -gt 0 ]] ; 
then
echo "previous frequencies, missingness, hwe, info scores files present from previous run which will be removed"
	rm -f ${section_02_dir}/data.afreq.gz
	rm -f ${section_02_dir}/data.hardy.gz
	rm -f ${section_02_dir}/data.info.gz
	rm -f ${section_02_dir}/data.vmiss.gz
	rm -f ${section_02_dir}/data.smiss.gz
else
	echo "passed file check"
fi
	

${plink2} \
	--bfile ${bfile} \
	--freq \
	--hardy \
	--missing \
	--out ${section_02_dir}/data

gzip -f -c ${quality_scores} > ${section_02_dir}/data.info.gz
gzip ${section_02_dir}/data.hardy
gzip ${section_02_dir}/data.smiss
gzip ${section_02_dir}/data.vmiss
gzip ${section_02_dir}/data.afreq

# Check missingness
missingness=`zcat ${section_02_dir}/data.smiss | awk '{ sum += $6; n++ } END { if (n > 0) print sum / n; }'`

echo "Average missingness: ${missingness}"

if (( $(bc <<< "${missingness} > 0.02") ))
then
	echo ""
	echo ""
	echo ""
	echo ""
	echo "WARNING"
	echo ""
	echo ""
	echo "Your genetic data has missingness of ${missingness}"
	echo ""
	echo "This seems high considering that you should have converted to best guess format with a very high hard call threshold"
	echo ""
	echo "Please ensure that this has been done"
fi

# Update ids
awk '{print $1,$2}' < ${bfile}.fam > ${intersect_ids_plink}
awk '{print $2}' < ${bfile}.fam > ${intersect_ids}

#rm -f ${bfile}.*~
rm temp_hm3snps.txt


echo "Successfully formatted SNP data"
