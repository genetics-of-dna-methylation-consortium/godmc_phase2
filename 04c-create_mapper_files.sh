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
exec &> >(tee ${section_04c_logfile})
print_version

#Please read resources/bin/hase/README_2.md
#An example is also provided below

mkdir -p ${hase_mapping}

python ${hase}/tools/mapper.py \
   -g ${hase_converting} \
   -o ${hase_mapping} \
   -study_name ${study_name} \
   -ref_name "ref-hrc"

echo "Successfully mapped the genetic data"
