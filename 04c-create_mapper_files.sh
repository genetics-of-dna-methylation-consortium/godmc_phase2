#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_04c_logfile})
print_version

#Please read resources/bin/light_hase/README_2.md
#An example is also provided below

mkdir -p ${light_hase_mapping}

# Mapper log interpretation:
# matched means ref allele1/allele2 match the converted PLINK .bim allele1/allele2
# order; flipped means the reverse order was seen. Since HASE/light_hase decode
# PLINK .bed as allele2 dosage, a straight match means HASE beta is relative to
# ref str_allele2, while ref str_allele1 is the other allele.

"${PYTHON_RUNNER[@]}" "${light_hase}/tools/mapper.py" \
   -g ${light_hase_converting} \
   -o ${light_hase_mapping} \
   -study_name ${study_name} \
   -ref_name "ref-hrc"

echo "Successfully mapped the genetic data"
