#!/bin/bash

find  *_allOutput.txt | while read f; do less $f | tr ' ' '\t' | grep -v 'pos0' | intersectBed -wo -a stdin -b ~/piquelab/data/RefTranscriptome/ensGene.hg19.v2.bed.gz -split | cut -f1-3,5,7-10,14,23-24 - > $f.tid; done
