#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_10_dir}/logs_b
touch ${section_10b_logfile}
exec &> >(tee ${section_10b_logfile})
print_version

#### GREML for SNP heritability##################################

tail -n +2 ${age_pred}.txt > ${age_pred}.plink
clock_names=$(cut -d" " -f 3- ${age_pred}.txt | head -n 1)

i=1
for clock_name in $clock_names
do
  echo "Runing SNP heritability for " $clock_name $i

  ${gcta} \
          --grm ${grmfile_all}_gaws10  \
          --reml \
          --mpheno $i \
          --pheno ${age_pred}.plink \
          --out ${section_10_dir}/heritability_${clock_name} \
          --thread-num ${nthreads}
  
  ${gcta} \
        --grm ${grmfile_all}_gaws10  \
        --reml \
        --mpheno $i \
        --pheno ${age_pred}.plink \
        --qcovar ${home_directory}/processed_data/genetic_data/gaws10_pc.eigenvec \
        --out ${section_10_dir}/heritability_${clock_name}_PCA \
        --thread-num ${nthreads}
  
  i=$(($i+1))
done

echo "Successfully finished the calculation on SNP heritability for age accelerations!"

