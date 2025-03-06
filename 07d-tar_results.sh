#!/bin/bash -l

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_07d_logfile_$1})
print_version

cd $home_directory

suff="tgz"
flags="czf"

i=$1
mkdir -p results/07/chr${i}
mv ./results/07/vQTL_drm_cis_*cpgchr${i}_1_1* ./results/07/chr${i}
mv ./results/07/vQTL_svlm_cis_*cpgchr${i}_1_1* ./results/07/chr${i}
mv ./results/07/vQTL_BF_cis_*cpgchr${i}* ./results/07/chr${i}
tar ${flags} ${home_directory}/results/${study_name}_07_chr${i}.${suff} ${home_directory}/results/07/chr${i}

echo "Successfully created results archives of module 07 chr ${i}"
