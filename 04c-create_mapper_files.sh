#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_04c_logfile})
print_version

#Please read resources/bin/light_hase/README_2.md
#An example is also provided below

mkdir -p ${light_hase_mapping}

python ${light_hase}/tools/mapper.py \
   -g ${light_hase_converting} \
   -o ${light_hase_mapping} \
   -study_name ${study_name} \
   -ref_name "ref-hrc"

echo "Successfully mapped the genetic data"
