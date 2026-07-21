#!/bin/bash

cd ..
source config

for i in $(seq 1 22)
do
sbatch --mem 128G 06e.vmeQTL_detection_candidateCpGs.sh BF $i
sbatch --mem 128G 06e.vmeQTL_detection_candidateCpGs.sh drm $i
sbatch --mem 128G 06e.vmeQTL_detection_candidateCpGs.sh svlm $i
done
