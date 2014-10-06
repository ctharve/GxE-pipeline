#!/bin/bash
plate=$1
cellLine=$2
script=$3
ncpus=2
queue=mmtxq
jobName=${plate}-${cellLine}

echo "cd ${PWD}; Rscript ${script} ${plate} ${cellLine}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
