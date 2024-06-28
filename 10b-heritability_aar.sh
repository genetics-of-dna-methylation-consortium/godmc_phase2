#!/bin/bash

set -e

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
mkdir -p ${section_10_dir}/logs_b
touch ${section_10b_logfile}
exec &> >(tee ${section_10b_logfile})
print_version

#### GREML for SNP heritability##################################

tail -n +2 ${age_pred}.txt > age_acc.plink
clock_names=$(cut -d" " -f 3- ${age_pred}.txt | head -n 1)

i=1
for clock_name in $clock_names
do
  echo "Runing SNP heritability for " $clock_name $i
  ${gcta} \
          --grm ${grmfile_all}  \
          --reml \
          --mpheno $i \
          --pheno age_acc.plink \
          --out ${section_10_dir}/heritability_${clock_name} \
          --thread-num ${nthreads}
  
  i=$(($i+1))
done
rm age_acc.plink

echo "Successfully finished the calculation on SNP heritability for age accelerations!"

