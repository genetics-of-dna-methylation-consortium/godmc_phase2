#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_11_dir}/logs_b
touch ${section_11b_logfile}
exec &> >(tee ${section_11b_logfile})
print_version

#### GREML for SNP heritability###################################

  ${gcta} \
          --grm ${grmfile_all}_gaws10  \
          --reml \
          --pheno ${smoking_pred}.smok.plink \
          --out ${section_11_dir}/heritability_smoking \
          --thread-num ${nthreads}
          
echo "Successfully finished the calculation on SNP heritability for smoking!"          
        