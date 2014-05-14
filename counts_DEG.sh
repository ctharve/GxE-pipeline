#!/bin/bash
plate=$1
bamFolder=$2
ncpus=12
queue=mmtxq
jobName=DEG_counts_${plate}

# sh counts_QC.sh P4 /wsu/home/fl/fl97/fl9788/piquelab/charvey/system/derived_data/P4/bams

echo "cd ${PWD}; Rscript counts_DEG.R ${ncpus} ${plate} ${bamFolder}" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 

