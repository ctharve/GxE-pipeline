#!/bin/bash
# execute within ~/wsu/home/groups/piquelab/charvey/GxE/jointGenotyping/QuASAR_results_<somePlateHere>/output

find  *_T1[36]C1_allOutput.txt | while read f; do less $f | tr ' ' '\t' | grep -v 'pos0' | intersectBed -a stdin -b ~/piquelab/data/RefTranscriptome/ensGene.hg19.v2.bed.gz -wa -wb -split | cut -f1-8,12,21- > $f.tid; done
