#!/bin/bash
set -v
set -e
plate=$1

## results directory & QuASAR scripts
mkdir -p ../QuASAR_results_${plate}/data
cp ./QuASAR_pipeline/QuASAR_pipeline* ../QuASAR_results_${plate}

## Mesh directory & Mesh scripts
mkdir -p ../QuASAR_results_${plate}/mesh/analysis_data/
cp ./Mesh_pipeline/* ../QuASAR_results_${plate}/mesh/analysis_data/

## script for annotating scripts with gene IDs
mkdir -p ../QuASAR_results_${plate}/output/
cp ./Annotation_pipeline/add_annotations.sh ../QuASAR_results_${plate}/output/

## link to all files for analysis
cd ../QuASAR_results_${plate}/data
ls /wsu/home/groups/piquelab/charvey/GxE/derived_data/${plate}/pileups/*.clean.bed.gz | while read f; do ln -sv $f; done
