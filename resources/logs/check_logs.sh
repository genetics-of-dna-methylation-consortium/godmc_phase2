#!/usr/bin/env bash
vercomp () {
	if [[ $1 == $2 ]]
	then
		echo "Correct script version"
		return 0
	fi
	local IFS=.
	local i ver1=($1) ver2=($2)
	# fill empty fields in ver1 with zeros
	for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
	do
		ver1[i]=0
	done
	for ((i=0; i<${#ver1[@]}; i++))
	do
		if [[ -z ${ver2[i]} ]]
		then
			# fill empty fields in ver2 with zeros
			ver2[i]=0
		fi
		if ((10#${ver1[i]} > 10#${ver2[i]}))
		then
			echo "Script version greater than required"
			return 0
		fi
		if ((10#${ver1[i]} < 10#${ver2[i]}))
		then
			echo ""
			echo "PROBLEM"
			echo "This analysis was performed on an outdated script."
			echo "Expecting at least version $2"
			echo "But the logs show that this was run on version $1"
			echo "Please run 'git pull' and then re-run the analysis."
			echo ""
			return 1
		fi
	done
	echo "Correct script version"
	return 0
}

compare_version () {

    logfile="section_${1}_logfile"
    version_used=$(grep "GoDMC2 version" ${!logfile}* | head -n 1 | cut -d " " -f 3)
	if [ "${version_used}" = "" ]
	then
		echo ""
		echo "WARNING"
		echo "No version number found. You are probably running an old version of git."
		echo "The scripts you used could be out of date."
		echo "Please run 'git pull' and check that no updates were made to the ${1} script you are checking."
		echo "If updates were made then please re-run this ${1} script."
		echo ""
		return 0
	fi

	version_required=$(grep "section_${1}" resources/logs/versions.txt | cut -d " " -f 2)
	echo "Version required: ${version_required}"
	echo "Version used: ${version_used}"
	vercomp ${version_used} ${version_required}

}


check_logs_01 () {

	compare_version "01"
	if grep -i -q "success" ${section_01_logfile}; then
		echo "01-check_data.sh completed successfully."
	else
		echo "Problem: 01-check_data.sh did not complete successfully"
		exit 1
	fi

}


check_logs_02 () {

	compare_version "02a"
	if grep -i -q "Successfully" ${section_02a_logfile}; then
		echo "02a-snp_data.sh completed successfully."
	else
		echo "Problem: 02a-snp_data.sh did not complete successfully"
		echo "Please ensure previous steps have been performed and rerun 02a-snp_data.sh"
		exit 1
	fi

	compare_version "02b"
	if grep -i -q "Successfully" ${section_02b_logfile}; then
		echo "02b-convert_snp_format.sh completed successfully."
	else
		echo "Problem: 02b-convert_snp_format.sh did not complete successfully"
		echo "Please ensure previous steps have been performed and rerun 02b-convert_snp_format.sh"
		exit 1
	fi

}

check_logs_03a () {

  compare_version "03a"
	if grep -i -q "success" ${section_03a_logfile}; then
		echo "03a-methylation_variables.sh completed successfully."
	else
		echo "Problem: 03a-methylation_variables.sh did not complete successfully"
		exit 1
	fi
}


check_logs_03d () {
    compare_version "03d"
    if grep -i -q "success" ${section_03d_logfile}; then
        echo "03d-non_genetic_methylation_pcs.sh completed successfully."
    else
        echo "Problem: 03d-non_genetic_methylation_pcs.sh did not complete successfully"
        exit 1
    fi
}

check_logs_03 () {

    compare_version "03a"
	if grep -i -q "success" ${section_03a_logfile}; then
		echo "03a-methylation_variables.sh completed successfully."
	else
		echo "Problem: 03a-methylation_variables.sh did not complete successfully"
		exit 1
	fi

	compare_version "03b"
	count_03b=`grep -i "success" ${section_03b_logfile}[0-9]* | wc -l`
	if [ ${count_03b} == ${meth_chunks} ]; then
		echo "03b-methylation_adjustment1.sh completed successfully."
	else
		echo "Problem: 03b-methylation_adjustment1.sh did not complete successfully"
		exit 1
	fi

	compare_version "03c"
	if grep -i -q "success" ${section_03c_logfile}; then
		echo "03c-methylation_pcs.sh completed successfully."
	else
		echo "Problem: 03c-methylation_pcs.sh did not complete successfully"
		exit 1
	fi

    compare_version "03d"
	if grep -i -q "success" ${section_03d_logfile}; then
		echo "03d-non_genetic_methylation_pcs.sh completed successfully."
	else
		echo "Problem: 03d-non_genetic_methylation_pcs.sh did not complete successfully"
		exit 1
	fi

	compare_version "03e"
	count_03e=`grep "Successfully adjusting non-genetic PCs" ${section_03e_logfile}[0-9]* | wc -l`
	if [ ${count_03e} == ${meth_chunks} ]; then
		echo "03e-methylation_adjustment2.sh completed successfully."
	else
		echo "Problem: 03e-methylation_adjustment2.sh did not complete successfully"
		exit 1
	fi

	compare_version "03f"
	if grep -i -q "success" ${section_03f_logfile}*; then
		echo "03f-convert_methylation_format.sh completed successfully."
	else
		echo "Problem: 03f-convert_methylation_format.sh did not complete successfully"
		exit 1
	fi

	compare_version "03g"
	if grep -i -q "success" ${section_03g_logfile}*; then
		echo "03g-perform_positive_control.sh completed successfully."
	else
		echo "Problem: 03g-perform_positive_control.sh did not complete successfully"
		exit 1
	fi

}

check_logs_04 () {

    compare_version "04a"
	if grep -i -q "Successfully" ${section_04a_logfile}; then
		echo "04a-convert_snp_format.sh completed successfully."
	else
		echo "Problem: 04a-convert_snp_format.sh did not complete successfully"
		exit 1
	fi

	compare_version "04b"
	if grep -i -q "Allele positions inverted" ${section_04b_logfile}; then
		echo "04b-mapper-preparation.sh completed successfully."
	else
		echo "Problem: 04b-mapper-preparation.sh did not complete successfully"
		exit 1
	fi

	compare_version "04c"
	if grep -i -q "Successfully" ${section_04c_logfile}; then
		echo "04c-create_mapper_files.sh completed successfully."
	else
		echo "Problem: 04c-create_mapper_files.sh did not complete successfully"
		exit 1
	fi

	compare_version "04d"
	if grep -i -q "Successfully" ${section_04d_logfile}; then
		echo "04d-encoding.sh completed successfully."
	else
		echo "Problem: 04d-encoding.sh did not complete successfully"
		exit 1
	fi

	compare_version "04e"
	if grep -i -q "success" ${section_04e_logfile}; then
		echo "04e-single-site-analysis.sh completed successfully."
	else
		echo "Problem: 04e-single-site-analysis.sh did not complete successfully"
		exit 1
	fi

    compare_version "04f"
    if grep -i -q "Successfully created results archives of module 04" ${section_04f_logfile}; then
        echo "04f-tar_results.sh completed successfully."
    else
        echo "Problem: 04f-tar_results.sh did not complete successfully"
        exit 1
    fi
}

check_logs_07 () {

	compare_version "07a"
	if grep -i -q "Success" ${section_07a_logfile}; then
		echo "07a-vmeQTL_data_preparation.sh completed successfully."
	else
		echo "Problem: 07a-vmeQTL_data_preparation.sh did not complete successfully"
		exit 1
	fi

	compare_version "07c"
	if grep -i -q "using drm method to detect vmeQTLs is successful" ${section_07c_logfile}_drm &
	   grep -i -q "using svlm method to detect vmeQTLs is successful" ${section_07c_logfile}_svlm &
	   grep -i -q "using BF method to detect vmeQTLs is successful" ${section_07c_logfile}_BF; then
		echo "07b-run_cis_vmeQTL.sh and 07c-check_vmeQTL_results.sh completed successfully."
	else
		echo "Problem: 07c-check_vmeQTL_results.sh did not complete successfully"
		exit 1
	fi
    
    compare_version "07d"
    count_07d=`grep -i "success" ${section_07d_logfile}_[0-9]* | wc -l`
    if [ ${count_07d} == 22 ]; then
        echo "07d-methylation_adjustment1.sh completed successfully."
    else
        echo "Problem: 07d-methylation_adjustment1.sh did not complete successfully"
        exit 1
    fi
   
}

check_logs_08 () {

	compare_version "08a"
	if grep -i -q "Success" ${section_08a_logfile}; then
		echo "08a-genotypeInversion completed successfully."
	else
		echo "Problem: 08a-genotypeInversion did not complete successfully"
		exit 1
	fi

	compare_version "08b"
	if grep -i -q "Success" ${section_08b_logfile}; then
		echo "08b-inversionmeQTL.sh completed successfully."
	else
		echo "Problem: 08b-inversionmeQTL.sh did not complete successfully"
		exit 1
	fi


}

check_logs_14 () {

	compare_version "14"
	if grep -i -q "Success" ${section_14_logfile}; then
		echo "14-nc886_gwas.sh completed successfully."
	else
		echo "Problem: 14-nc886_gwas.sh did not complete successfully"
		exit 1
	fi

}


check_logs_09 () {

  compare_version "09"
   
  vect_PRS=$(grep "PRS" ${scripts_directory}/resources/parameters | grep "weights" | awk -F"_" '{print $2}' |tr "\n" " ")
  vect_PRS_array=($vect_PRS)
  n=$((${#vect_PRS_array[*]}-1))

  for ((k=0;k<=$n;k++))
  do

    PRS=${vect_PRS_array[$k]}
    log_file=${section_09_dir}/${PRS}/logs/log.txt
 
	  if grep -i -q "run successfully" ${log_file}; then
		  echo "09 completed successfully for ${PRS}"
	  else
		  echo "Problem: 09 did not complete successfully for ${PRS}"
		  exit 1
	  fi

  done

	if grep -i -q "script finalised" ${section_09_logfile}; then
		echo "09 completed successfully for all traits"
	else
		echo "Problem: 09 did not complete successfully for all traits"
		exit 1
	fi


}


