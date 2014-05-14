#!/bin/bash
set -v 
set -e

plate=$1 

mkdir -p ../$plate/fastqs
mkdir -p ../$plate/bams
mkdir -p ../$plate/counts
mkdir -p ../$plate/counts/QC

cp Makefile ../$plate/bams
cp counts_DEG.* ../$plate/counts
cp counts_QC.* ../$plate/counts/QC

cd ../$plate/fastqs

for bc in {1..96}; do 
	find ${LPG}/OurData/ -name "*${plate}-HT${bc}_*R1*.fastq.gz" \
		| grep -v 'D[0-9]*${plate}' \
		| awk '{print $0,NR}' \
		| while read f r; do  
			b=${f//R1/R2}; 
			echo "${plate}-HT${bc}_S${r}";
			ln -s $f ${plate}-HT${bc}_S${r}_R1.fastq.gz; 
			ln -s $b ${plate}-HT${bc}_S${r}_R2.fastq.gz; 
		done; 
done;
