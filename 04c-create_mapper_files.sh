#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

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
