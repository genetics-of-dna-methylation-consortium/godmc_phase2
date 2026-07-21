#!/bin/bash

cd ..
source ./config

for i in $(seq 1 22);
do 
    sbatch --mem 32G ${scripts_directory}/06b.vmeQTL_detection_missingCpGs.sh $i
done
