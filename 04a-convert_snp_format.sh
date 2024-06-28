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
exec &> >(tee ${section_04a_logfile})
print_version
#Please read reasources/bin/hase/README_2.md
#An expample is also provided below

mkdir -p ${hase_dir_in}
mkdir -p ${hase_converting}
cp ${bfile}.bim	${hase_dir_in}
cp ${bfile}.fam ${hase_dir_in}
cp ${bfile}.bed ${hase_dir_in}

python ${hase}/hase.py \
    -mode converting \
    -g ${hase_dir_in} \
    -o ${hase_converting} \
    -study_name ${study_name} # the name for your study 

echo "Successfully converted genetic data"
