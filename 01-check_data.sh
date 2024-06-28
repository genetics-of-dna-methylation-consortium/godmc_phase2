#!/bin/bash


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


exec &> >(tee ${section_01_logfile})
exec &> >(tee ${section_01_logfiles})
print_version


containsElement () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo "There is no method for ${1}."
	echo "Please run:"
	echo "./01-check_data [arg]"
	echo "where arg is an optional argument that can be one of:"
	printf '%s\n' ${@:2}
	return 1
}

arg="all"
declare -a sections=('all' 'config' 'download' 'requirements' 'genetic' 'methylation' 'covariates' 'covariates_for_PRS' 'phenotypes_for_PRS' 'summary')

vect_PRS=$(grep "PRS" ${scripts_directory}/resources/parameters | grep "weights" | awk -F"_" '{print $2}' |tr "\n" " ")

if [ -n "${1}" ]; then
	arg="${1}"
	containsElement ${1} ${sections[@]}
fi





section_message () {

	echo "-----------------------------------------------"
	echo ""
	echo "$1 section"
	echo ""
	echo "to run this part on its own type:"
	echo "$ ./01-check_data.sh $1"
	echo ""
	echo "-----------------------------------------------"
	echo ""
	echo ""

}


if [ "$arg" = "config" ] || [ "$arg" = "all" ]
then

	section_message "config"

	if ! [[ "$study_name" =~ [^a-zA-Z0-9_\] ]] && ! [ "$study_name" = "" ]
	then
		echo "The study name '${study_name}' will be used for this analysis. Change this in the config file if necessary."
		echo ""
	else
		echo "The study name '${study_name}' is invalid. Please use only alphanumeric or underscore characters, with no spaces or special characters etc."
		exit
	fi

fi


if [ "$arg" = "download" ] || [ "$arg" = "all" ]
then

	section_message "download"

	sftp -P 2222 -oIdentityFile=$key ${sftp_username}@${sftp_address}:${sftp_path} <<EOF
get HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz
get HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz.md5sum
get ref-hrc.ref.gz
get ref-hrc.ref_info.h5
EOF

	mv HRC.r1-1.GRCh37.wgs.mac5.sites.tab.cptid.maf001_recoded.gz* ${scripts_directory}/resources/genetics
	mv ref-hrc.ref.gz ${hase}/data
    mv ref-hrc.ref_info.h5 ${hase}/data
fi


if [ "$arg" = "requirements" ] || [ "$arg" = "all" ]
then

	section_message "requirements"

    ${R_directory}Rscript resources/datacheck/requirements.R $related
fi

if [ "$arg" = "genetic" ] || [ "$arg" = "all" ]
then

	section_message "genetic"


	${R_directory}Rscript resources/datacheck/genetic_data.R \
		${bfile_raw}.bim \
		${bfile_raw}.fam \
		${quality_scores} \
		${control_snps} \
		${snpchrtxt} \
		${snpchrplot} \
		${genetic_descriptives} \
		${quality_scores_plot}

	# Check missingness, there should be a small percentage of missingness (--hard-call-threshold 0.499999)
	${plink2} --bfile ${bfile_raw} --missing --out ${section_01_dir}/data
	gzip ${section_01_dir}/data.smiss
	gzip ${section_01_dir}/data.vmiss

	nrow=`zcat ${section_01_dir}/data.smiss.gz | awk 'NR>1 && $5>0.02 {print $0}'  |wc -l`

	nrow_all=`zcat ${section_01_dir}/data.smiss.gz | awk 'NR>1 {print $0}' |wc -l`
 	prop=$( printf '%.2f' $(echo "$nrow / $nrow_all" | bc -l) )
	echo "Proportion with >= 0.02 missingness: $prop"

	if [ "$prop" \> 0.01 ]
	then
		echo "Error: $prop of genotypes have more than 2% of missing values. Please don't use a genotype probability cut-off when generating best guess data."
		exit 1
	else
		echo "\n"
		echo "Best guess data appears to be correct"
	fi

fi

if [ "$arg" = "methylation" ] || [ "$arg" = "all" ]
then

	section_message "methylation"


	${R_directory}Rscript resources/datacheck/methylation_data.R \
		${betas} \
		${bfile_raw}.fam \
		${sorted_methylation} \
		${methylation_array} \
		${measured_cellcounts} \
		${meth_ids} \
		${methylation_descriptives} \
		${methylation_summary} \
    ${covariates} \
    ${sex_pred_plot} \
		${intersect_ids} \
		${intersect_ids_plink}
    
fi

if [ "$arg" = "covariates" ] || [ "$arg" = "all" ]
then

	section_message "covariates"


	${R_directory}Rscript resources/datacheck/covariates.R \
		${covariates} \
		${bfile_raw}.fam \
		${meth_ids} \
		${sorted_methylation} \
		${ageplot} \
		${covariate_descriptives}
fi

if [ "$arg" = "covariates_for_PRS" ] || [ "$arg" = "all" ]
then

	section_message "covariates_for_PRS"

  for PRS in $vect_PRS
  do

	echo "processing $PRS"
  echo ""
  
  covar_file_PRS=covariates_$PRS
  covar_desc_PRS=covariate_${PRS}_PRS_descriptives
  
  ${R_directory}Rscript resources/datacheck/covariates_for_PRS.R \
    ${!covar_file_PRS} \
    ${bfile_raw}.fam \
    ${meth_ids} \
    ${!covar_desc_PRS}
  done
fi

if [ "$arg" = "phenotypes_for_PRS" ] || [ "$arg" = "all" ]
then

  section_message "phenotypes_for_PRS"

  for PRS in $vect_PRS
  do

	echo "processing $PRS"
  echo ""
  
  pheno_file_PRS=phenotypes_$PRS
  pheno_desc_PRS=phenotype_${PRS}_PRS_descriptives
  
  ${R_directory}Rscript resources/datacheck/phenotypes_for_PRS.R \
    ${!pheno_file_PRS} \
    ${bfile_raw}.fam \
    ${!pheno_desc_PRS}
  done
fi

if [ "$arg" = "summary" ] || [ "$arg" = "all" ]
then

	section_message "summary"


		${R_directory}Rscript resources/datacheck/collect_descriptives.R \
		${genetic_descriptives} \
		${methylation_descriptives} \
		${covariate_descriptives} \
		${phenotype_descriptives} \
		${cohort_descriptives}
fi


if [ "$arg" = "all" ]
then
	echo ""
	echo ""
	echo "You successfully performed all data checks!"
fi
