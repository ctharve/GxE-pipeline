#!/bin/bash

set -v 
set -e

cl=$1
thresh=( .5 .9 .99 )

for ii in ${thresh[@]}; do
    echo $ii;
    less ../output/QuASAR_genotypes_noPileupMin_${cl}_minCov15.txt | awk -v myii=$ii '(1-$8)>myii' | wc -l;
done;
