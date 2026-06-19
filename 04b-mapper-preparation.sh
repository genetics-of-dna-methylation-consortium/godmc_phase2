#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_04b_logfile})
print_version

#Note that this step might not be neccessary
#Please read resources/bin/hase/README_2.md
#An example is also provided below

echo "Skipping old 04b because 04a aligns alleles to ref-hrc before HASE converting."

# python  ${hase}/added/invert_probes.py \
#     -f  ${hase_converting}/probes \
#     -n ${study_name}

echo "No probe inversion performed"
