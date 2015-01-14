#!/bin/bash
plate=$1
target=/wsu/home/groups/piquelab/scratch/wwwShare/charvey/system/deep/${plate}
mkdir -p ${target}
cd ../QuASAR_results_${plate}/
cp -rf ./output ${target} 
cp -rf ./plots ${target}
