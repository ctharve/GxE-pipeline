#!/bin/bash
plate=$1
bamFolder=$2
ncpus=16
jobName=QC_counts

echo "cd ${PWD}; Rscript counts_QC.R ${plate} ${bamFolder}" | qsub -q masxq -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 

