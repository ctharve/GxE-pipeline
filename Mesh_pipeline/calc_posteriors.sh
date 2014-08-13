#!/bin/bash
ncpus=4
plate=$1
queue=mmtxq
jobName=clPst_${plate}

# sh counts_QC.sh P4 /wsu/home/fl/fl97/fl9788/piquelab/charvey/system/derived_data/P4/bams

echo "cd ${PWD}; Rscript calc_posteriors.R" | qsub -q ${queue} -l nodes=1:ppn=${ncpus} -N ${jobName} -o ${jobName}.Qsub -e ${jobName}.Qsub.e 
