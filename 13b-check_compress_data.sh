#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_13_dir}/logs_b
touch ${section_13b_logfile}
exec &> >(tee ${section_13b_logfile})
print_version


#### fucntions for checking log files ############## 


if grep -i -q "Successfully finished the GWAS on MZepigeneticsignature" ${section_13a_logfile}; then
		echo "13 GWAS MZepigeneticsignature completed successfully."
	else
		echo "Problem: 13a-gwas_MZepi.sh did not complete successfully"
		exit 1


fi

#### Compress the output files from 13 ############## 
echo ""
echo "Compressing the outputs 13"
cd ${home_directory}

mkdir -p ${section_13_dir}/upload
cd ${section_13_dir}
mv  *.pdf upload
mv  *.log upload
mv  *.fastGWA upload
mv  *.RData upload
mv  *.jpeg upload
cp -r logs_* upload
cd ${home_directory}
tar -zcf results/MZepiGWAS_module13_${study_name}.tgz results/13/upload

echo "Successfully created results archives ${home_directory}/results/13/MZepiGWAS_module13_${study_name}.tgz"


# Generating md5 checksum for verify the data intensity
cd ${home_directory}/results || exit 1
md5sum MZepiGWAS_module13_${study_name}.tgz > MZepiGWAS_module13_${study_name}.tgz.md5sum
md5sum -c MZepiGWAS_module13_${study_name}.tgz.md5sum
# encryption 
gpg --output MZepiGWAS_module13_${study_name}.tgz.gpg --symmetric --cipher-algo AES256 MZepiGWAS_module13_${study_name}.tgz
echo ""
echo "Please download the following files to your own local machine and upload to https://drive.google.com/drive/folders/1Bir6C8H6zh2Li6_6SCBz4LmqaTWvhrbH?usp=drive_link"
echo "1. " ${home_directory}/results/MZepiGWAS_module13_${study_name}.tgz.md5sum
echo "2. " ${home_directory}/results/MZepiGWAS_module13_${study_name}.tgz.gpg
echo "Please share encryption passphrase to the developers by emailing the developer j.van.dongen@vu.nl mentioned in the wiki."
echo "Thank you very much for contributing to this GWAS!"

