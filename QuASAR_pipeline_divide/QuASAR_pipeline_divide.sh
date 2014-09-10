#!/bin/bash
plate=$1
cutOff=$2

less ~/piquelab/charvey/GxE/derived_data/covariates/GxE_${plate}_covariates.txt \
	| awk '$13!="CellLine" {print $13}' \
	| sort \
	| uniq \
	| while read cl; do Rscript QuASAR_pipeline_divideLFC.R ${plate} $cl ${cutOff}; done 

echo $plate cutOff = ${cutOff} > runParms.txt 
