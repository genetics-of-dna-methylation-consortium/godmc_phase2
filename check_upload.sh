#!/bin/bash

setup_args=()
section_args=()
while [ "$#" -gt 0 ]; do
	case "$1" in
		-c)
			setup_args+=("$1" "$2")
			shift 2
			;;
		*)
			section_args+=("$1")
			shift
			;;
	esac
done

source resources/setup.sh "${setup_args[@]}" "${section_args[@]}"
set -- "${section_args[@]}"

checkFirstArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $1 is not a valid section identifier"
	echo $"Need to specify a value from 01,02,03,03a,03d,04,07,08,09,10,11,14,15"
	echo $"Usage: $0 <pipeline section> {check|upload}"
	exit 1
}

checkSecondArg () {
	local e
	for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
	echo $"Error: $2 is not a valid action"
	echo $"Specify 'check', 'upload', or 'cleanup' (section 15 only)"
	echo $"Usage: $0 <pipeline section> {check|upload|cleanup}"
	exit 1
}

source resources/logs/check_logs.sh
source resources/logs/check_results.sh

require_section_15_upload_config () {
	if [ -n "${LD_SHIP_CMD:-}" ]; then
		return 0
	fi
	if [ -z "${GODMC_UPLOAD_PASSWORD:-}" ]; then
		echo "Problem: GODMC_UPLOAD_PASSWORD is required for section-15 upload" >&2
		return 1
	fi
	if ! command -v curl >/dev/null 2>&1; then
		echo "Problem: curl is required for section-15 upload" >&2
		return 1
	fi
}

ship_section_15_file () {
	local path="$1" filename http_status
	if [ -n "${LD_SHIP_CMD:-}" ]; then
		"${LD_SHIP_CMD}" "${path}"
		return $?
	fi
	filename="$(basename "${path}")"
	if ! http_status="$(curl --fail --show-error \
		--output /dev/null \
		--write-out "%{http_code}" \
		--user "uploader:${GODMC_UPLOAD_PASSWORD}" \
		--upload-file "${path}" \
		"${godmc_upload_url}/${filename}")"; then
		return 1
	fi
	if ! [[ "${http_status}" =~ ^2[0-9][0-9]$ ]]; then
		echo "Problem: section-15 upload returned HTTP ${http_status} for ${filename}" >&2
		return 1
	fi
}

cleanup_section_15 () {
	local chr chunk chunks outdir path base scaffold targets removed=0 kept=0
	local section_15_upload_dir="${section_15_dir}/upload"
	local -a section_15_files
	if ! targets="$(ld_target_chromosomes_15)"; then
		exit 1
	fi
	while IFS= read -r chr; do
		outdir="$(ld_resolve_chromosome_dir_15 "${ld_prepare_dir}" "${chr}")"
		scaffold="${study_name}_chr${chr}_15_scaffold"
		section_15_files=("${section_15_upload_dir}/${scaffold}.tgz.aes")
		if ! chunks="$(ld_manifest_chunks_15 "${outdir}" "${chr}")"; then
			echo "Problem: failed to read the section-15 chunk manifest for chr${chr}" >&2
			exit 1
		fi
		while IFS= read -r chunk; do
			chunk="${study_name}_chr${chr}_15_chr${chr}_${chunk}"
			section_15_files+=("${section_15_upload_dir}/${chunk}.tgz.aes")
		done <<< "${chunks}"
		for path in "${section_15_files[@]}"; do
			if [ ! -f "${path}" ]; then
				continue
			fi
			base="$(basename "${path%.tgz.aes}")"
			if ld_upload_receipt_present_15 "${section_15_upload_dir}" "${base}"; then
				rm -f "${path}" || exit 1
				removed=$((removed + 1))
			else
				echo "Problem: ${base}.tgz.aes has no upload record; keeping it" >&2
				kept=$((kept + 1))
			fi
		done
	done <<< "${targets}"
	echo "Section 15 cleanup removed ${removed} uploaded archive(s)"
	if [ "${kept}" -ne 0 ]; then
		echo "Problem: ${kept} section 15 archive(s) are not recorded as uploaded; re-run the upload step" >&2
		exit 1
	fi
	echo "Successfully reclaimed section 15 uploaded archives"
}

upload_section_15 () {
	local chr chunk chunks outdir path upload_receipt scaffold targets failed=0 uploaded=0 section_15_upload_dir="${section_15_dir}/upload"
	local -a section_15_files
	require_section_15_upload_config || exit 1
	if ! targets="$(ld_target_chromosomes_15)"; then
		exit 1
	fi
	mkdir -p "${section_15_upload_dir}"
	shopt -s nullglob
	while IFS= read -r chr; do
		outdir="$(ld_resolve_chromosome_dir_15 "${ld_prepare_dir}" "${chr}")"
		scaffold="${study_name}_chr${chr}_15_scaffold"
		section_15_files=(
			"${section_15_upload_dir}/${scaffold}.tgz.aes"
			"${section_15_upload_dir}/${scaffold}.md5sum"
		)
		if ! chunks="$(ld_manifest_chunks_15 "${outdir}" "${chr}")"; then
			echo "Problem: failed to read the section-15 chunk manifest for chr${chr}" >&2
			exit 1
		fi
		while IFS= read -r chunk; do
			chunk="${study_name}_chr${chr}_15_chr${chr}_${chunk}"
			section_15_files+=(
				"${section_15_upload_dir}/${chunk}.tgz.aes"
				"${section_15_upload_dir}/${chunk}.md5sum"
			)
		done <<< "${chunks}"
		for path in "${section_15_files[@]}"; do
			upload_receipt="${section_15_upload_dir}/.uploaded_$(basename "${path}")"
			if [ -f "${upload_receipt}" ]; then
				echo "Section 15 upload already recorded for $(basename "${path}")"
				uploaded=$((uploaded + 1))
				continue
			fi
			if ship_section_15_file "${path}"; then
				touch "${upload_receipt}"
				echo "Uploaded section 15 file $(basename "${path}")"
				uploaded=$((uploaded + 1))
			else
				echo "Problem: failed to upload section 15 file ${path}" >&2
				failed=1
			fi
		done
	done <<< "${targets}"
	shopt -u nullglob
	if [ "${uploaded}" -eq 0 ]; then
		echo "Problem: no section 15 upload files found in ${section_15_upload_dir}"
		exit 1
	fi
	if [ "${failed}" -ne 0 ]; then
		exit 1
	fi
	echo "Successfully uploaded section 15 staged artefacts"
}

sections=("01" "02" "03" "03a" "03d" "04" "07" "08" "09" "10" "11" "14" "15")
checkFirstArg "$1" "${sections[@]}"

actions=("check" "upload" "cleanup")
checkSecondArg "$2" "${actions[@]}"

if [[ "$2" = "cleanup" && "$1" != "15" ]]
then
	echo $"Error: the cleanup action is only supported for section 15"
	exit 1
fi

if [[ "$2" = "cleanup" && "$1" = "15" ]]
then
	cleanup_section_15
	exit 0
fi

echo ""
echo "Checking log files for $1"
eval "check_logs_$1"

echo ""
echo "Checking results for $1"
eval "check_results_$1"

echo ""
echo "Section $1 has been successfully completed!"

if [[ "$2" = "upload" && ( $1 = "01" || $1 = "02" || $1 = "03" || $1 = "03a" || $1 = "03d" || $1 = "04" || $1 = "07" || $1 = "08" ) ]]
then

	echo ""
	temp=`which sshpass 2>/dev/null | wc -l`
	port="-P 2222"
	if [[ ! "${temp}" = "0" &&  $1 != "08" ]]
	then
		echo "sshpass detected"
		sftp $port -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path} << !
bye
!
	elif [[ ! "${temp}" = "0" && $1 = "08" ]]
	then
		read -s -p "Enter SFTP password: " mypassword
	    export SSHPASS=${mypassword}
		sshpass -e sftp -P 22 -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username_inversions}@${sftp_address_inversions}:${sftp_path_inversions} << !
	
bye
!
		echo "Connection established"
	else
		echo "sshpass is not installed."
		echo "The results will now be archived, once that is done they will be uploaded to the server"
	fi

	echo ""
	echo "Tarring results, log and config files"
	mkdir -p ${home_directory}/results/config/
	if [[ "$config_file" = /* ]]; then
		cp ${config_file} ${scripts_directory}/
		config_basename=$(basename "$config_file")
		tar czf ${home_directory}/results/config/${study_name}_config.tar -C ${scripts_directory}/ ./${config_basename} ./resources/parameters
		rm ${scripts_directory}/${config_basename}
	else
		tar czf ${home_directory}/results/config/${study_name}_config.tar -C ${scripts_directory}/ ./${config_file} ./resources/parameters 
	fi	
    echo "Successfully created config archive"	

    echo "Generating md5 checksum"
    md5sum ${home_directory}/results/config/${study_name}_config.tar > ${home_directory}/results/config/${study_name}_config.md5sum
    echo "Encrypting files"
    gpg --output ${home_directory}/results/config/${study_name}_config.tar.aes --symmetric --cipher-algo AES256 ${home_directory}/results/config/${study_name}_config.tar
    echo ""

	if [ "${1}" = "17" ]
	then
		suff="tar"
		flags="cf"
	else
		suff="tgz"
		flags="czf"
	fi

    if [[ $1 = "03a" || $1 = "03d" ]]
    then
        tar ${flags} ${home_directory}/results/${study_name}_${1}.${suff} -C ${home_directory} results/03/
        echo "Successfully created results archives"
#        echo "Generating md5 checksum"
#        md5sum ${home_directory}/results/${study_name}_${1}.${suff} > ${home_directory}/results/${study_name}_${1}.md5sum
#        echo "Encrypting files"
#        gpg --output ${home_directory}/results/${study_name}_${1}.${suff}.aes --symmetric --cipher-algo AES256 ${home_directory}//results/${study_name}_${1}.${suff}
        echo ""
    elif [[ $1 = "04" || $1 = "07" ]]
    then
        echo "Tarring results have been generated"
    else
	    tar ${flags} ${home_directory}/results/${study_name}_${1}.${suff} -C ${home_directory} results/${1}
        echo "Successfully created results archives"
    fi
    
    if [[ "$1" = "07" ]]
    then
        echo "start to encrypt 07"
        for i in $(seq 1 22)
        do
        if [ -f ${home_directory}/results/${study_name}_07_chr${i}.md5sum ] && [ -f ${home_directory}/results/${study_name}_07_chr${i}.tgz.aes ]; then
            echo ""
        else
            echo "Generating md5 checksum chr${i}"
            md5sum ${home_directory}/results/${study_name}_07_chr${i}.tgz > ${home_directory}/results/${study_name}_07_chr${i}.md5sum
            echo "Encrypting files chr${i}"
            gpg --output ${home_directory}/results/${study_name}_07_chr${i}.tgz.aes --symmetric --cipher-algo AES256 ${home_directory}/results/${study_name}_07_chr${i}.tgz
            echo ""
        fi
        done
    else
        echo "Generating md5 checksum"
        md5sum ${home_directory}/results/${study_name}_${1}.${suff} > ${home_directory}/results/${study_name}_${1}.md5sum
        echo "Encrypting files"
        gpg --output ${home_directory}/results/${study_name}_${1}.${suff}.aes --symmetric --cipher-algo AES256 ${home_directory}//results/${study_name}_${1}.${suff}
        echo ""
    fi

if [[ ! "${temp}" = "0"  ]]
then
echo "Detecting sshpass"

	if [[  $1 = "08" ]]
	then
		echo "Inversions"
		export SSHPASS=${mypassword}
	    sshpass -e sftp -P 22 -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username_inversions}@${sftp_address_inversions}:${sftp_path_inversions} << !
   		dir
		cd ${sftp_username_inversions}
   		put ${home_directory}/results/${study_name}_${1}.md5sum
   		put ${home_directory}/results/${study_name}_$1.${suff}.aes
		put ${home_directory}/results/config/${study_name}_config.tar.aes
		put ${home_directory}/results/config/${study_name}_config.md5sum 
   		bye
!
    elif [[ $1 = "07" ]]
    then
        sftp $port -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path} << !
        dir
        cd ../upload
        put ${home_directory}/results/${study_name}_${1}_chr*.md5sum
        put ${home_directory}/results/${study_name}_$1_chr*.${suff}.aes
        put ${home_directory}/results/config/${study_name}_config.tar.aes
        put ${home_directory}/results/config/${study_name}_config.md5sum 
        bye
!
    else
        sftp $port -oIdentityFile=$key -oBatchMode=no -b - ${sftp_username}@${sftp_address}:${sftp_path} << !
        dir
        cd ../upload
        put ${home_directory}/results/${study_name}_${1}.md5sum
        put ${home_directory}/results/${study_name}_$1.${suff}.aes
		put ${home_directory}/results/config/${study_name}_config.tar.aes
		put ${home_directory}/results/config/${study_name}_config.md5sum 
        bye
!
	fi
else

read -s -p "Ready to upload? Press enter to continue: " anykey

if [[ $1 = "07" ]]
then
    sftp $port -oIdentityFile=$key ${sftp_username}@${sftp_address}:${sftp_path} <<EOF
    dir
    cd ../upload
    put ${home_directory}/results/${study_name}_${1}_chr*.md5sum
    put ${home_directory}/results/${study_name}_$1_chr*.${suff}.aes
    put ${home_directory}/results/config/${study_name}_config.tar.aes
    put ${home_directory}/results/config/${study_name}_config.md5sum
EOF
else
    sftp $port -oIdentityFile=$key ${sftp_username}@${sftp_address}:${sftp_path} <<EOF
    dir
    cd ../upload
    put ${home_directory}/results/${study_name}_${1}.md5sum
    put ${home_directory}/results/${study_name}_$1.${suff}.aes
	put ${home_directory}/results/config/${study_name}_config.tar.aes
	put ${home_directory}/results/config/${study_name}_config.md5sum 
EOF
fi


fi

fi

    
if [[ "$2" = "upload" && $1 = "09" ]]
then

  vect_PRS=$(grep "PRS" ${scripts_directory}/resources/parameters | grep "weights" | awk -F"_" '{print $2}' |tr "\n" " ")
	
  cd ${section_09_dir}

  for PRS in $vect_PRS
  do

    pheno_for_PRS=phenotypes_${PRS}
    pheno_desc_PRS=phenotype_${PRS}_PRS_descriptives

	  echo ""
	  echo "Tarring results and log files for ${PRS} PRS"
	  echo ""
 
    if [[ ${!pheno_for_PRS} != "NULL" ]]
    then

      tar -zcf PRS_${PRS}_${study_name}.tgz ./${PRS}/* ${!pheno_desc_PRS}
  
    else

      tar -zcf PRS_${PRS}_${study_name}.tgz ./${PRS}/*

    fi

    echo ""
    echo "Generating md5 checksum for ${PRS} PRS"
    echo ""

	  md5sum PRS_${PRS}_${study_name}.tgz > PRS_${PRS}_${study_name}.tgz.md5sum

    echo ""
	  echo "Encrypting files for ${PRS} PRS"
    echo ""

	  gpg --output PRS_${PRS}_${study_name}.tgz.gpg --symmetric --cipher-algo AES256 PRS_${PRS}_${study_name}.tgz

  done

  rm PRS*.tgz
  
  echo ""
  echo "Your results are now ready for upload in ${section_09_dir}"
  echo "Please upload .gpg and .md5sum files to google drive following the instructions in the wiki"
  echo "You have successfully run the pipeline, thank you so much! :)"
  echo ""

fi


if [[ "$2" = "upload" && $1 = "14" ]]
then
	sftp_username=${sftp_username_nc866}
	sftp_address=${sftp_address_nc866}
	sftp_path=${sftp_path_nc866}
	suff="tgz"
	flags="czf"
	
	echo ""
	echo "Tarring results and log files"

	tar ${flags} ${home_directory}/results/${study_name}_${1}.${suff} ${home_directory}/results/${1}
	echo "Successfully created results archives"
	echo "Generating md5 checksum"
	md5sum ${home_directory}/results/${study_name}_${1}.${suff} > ${home_directory}/results/${study_name}_${1}.md5sum
	echo "Encrypting files"
	gpg --output ${home_directory}/results/${study_name}_${1}.${suff}.aes --symmetric --cipher-algo AES256 ${home_directory}/results/${study_name}_${1}.${suff}

read -s -p "Ready to upload? Press enter to continue: " anykey
echo ""
read -s -p "Enter SFTP password: " mypassword

curl -k -u ${sftp_username}:${mypassword} -T ${home_directory}/results/config/${study_name}_config.tar.aes "${sftp_address}/${sftp_path}/${study_name}_${1}._config.tar.aes" 
curl -k -u ${sftp_username}:${mypassword} -T ${home_directory}/results/config/${study_name}_config.md5sum "${sftp_address}/${sftp_path}/${study_name}_${1}._config.md5sum"
curl -k -u ${sftp_username}:${mypassword} -T ${home_directory}/results/${study_name}_${1}.${suff}.aes "${sftp_address}/${sftp_path}/${study_name}_${1}.${suff}.aes" 
curl -k -u ${sftp_username}:${mypassword} -T ${home_directory}/results/${study_name}_${1}.md5sum "${sftp_address}/${sftp_path}/${study_name}_${1}.md5sum" <<EOF

EOF

fi

if [[ "$2" = "upload" && $1 = "15" ]]
then
	upload_section_15
fi
