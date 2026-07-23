#!/bin/bash

cd ..
source ./config

for i in {1..12} {16..20};
do
    sbatch --mem 32G ${scripts_directory}/06b-vmeQTL_detection_missingCpGs.sh $i
done
