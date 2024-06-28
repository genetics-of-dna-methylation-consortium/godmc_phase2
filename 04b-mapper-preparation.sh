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
exec &> >(tee ${section_04b_logfile})
print_version

#Note that this step might not be neccessary
#Please read resources/bin/hase/README_2.md
#An example is also provided below

python  ${hase}/added/invert_probes.py \
    -f  ${hase_converting}/probes \
    -n ${study_name}

echo "Allele positions inverted"
