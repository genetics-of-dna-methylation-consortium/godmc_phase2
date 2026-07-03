#!/bin/bash -l

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_05_dir}/logs_f"
exec &> >(tee "${section_05f_logfile}")
print_version

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

check_file() {
    if [ ! -f "$1" ]; then
        fail "Missing required file: $1"
    fi
}

if [ "${config_file:0:1}" = "/" ]; then
    config_to_archive="${config_file}"
else
    config_to_archive="${scripts_directory}/${config_file}"
fi

check_file "${config_to_archive}"
check_file "${scripts_directory}/resources/parameters"

if [ ! -d "${section_05_dir}/meta_inputs_female" ] && [ ! -d "${section_05_dir}/meta_inputs_male" ]; then
    fail "Missing Module 05 meta input directories. Please run 05c before 05f."
fi

cd "${home_directory}"

archive="${home_directory}/results/${study_name}_05.tgz"
checksum="${archive}.md5sum"
encrypted="${archive}.gpg"

rm -f "${archive}" "${checksum}" "${encrypted}"

echo "Compressing Module 05 results"
tar -zcf "${archive}" \
    "${config_to_archive}" \
    "${scripts_directory}/resources/parameters" \
    "results/05"

check_file "${archive}"

echo "Generating md5 checksum"
cd "${home_directory}/results"
md5sum "${study_name}_05.tgz" > "${study_name}_05.tgz.md5sum"
md5sum -c "${study_name}_05.tgz.md5sum"

# echo "Encrypting Module 05 archive"
# gpg --output "${study_name}_05.tgz.aes" \
#     --symmetric \
#     --cipher-algo AES256 \
#     "${study_name}_05.tgz"

# check_file "${checksum}"
# check_file "${encrypted}"

echo ""
echo "Module 05 archive successfully created and encrypted."
echo "Please upload these files to Google Drive: https://drive.google.com/drive/folders/1q1djBVG5ms-Ud3btZmIf_YHHnjgphxem?usp=share_link"
echo "1. ${checksum}"
echo "2. ${encrypted}"
