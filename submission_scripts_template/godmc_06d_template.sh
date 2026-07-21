#!/bin/bash

cd ..
source config

for i in $(seq 1 22)
do
sbatch --mem 128G 06d.vmeQTL_detection_candidateSNPs.sh BF $i
sbatch --mem 128G 06d.vmeQTL_detection_candidateSNPs.sh drm $i
sbatch --mem 128G 06d.vmeQTL_detection_candidateSNPs.sh svlm $i
done
