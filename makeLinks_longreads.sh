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
	#find /nfs/rprscratch/OurData/ -name "*${plate}-HT${bc}_*R1*.fastq.gz" \
	## DANGER the line below resulted from a mis-labeling of the original fastqs. They correspond to D2P1
	## If you run this script with another plate name they will be named with that prefix, so don't. Please & thank you.
	find ${LPG}/OurData/140404_SN7001329_0366_AH94AWADXX/FastQ/Project_FL/ -name "*Luca-HT${bc}_*R1*.fastq.gz" \
		| awk '{print $0,NR}' \
		| while read f r; do  
			b=${f//R1/R2}; 
			echo "${plate}-HT${bc}_L${r}";
			ln -s $f ${plate}-HT${bc}_L${r}_R1.fastq.gz; 
			ln -s $b ${plate}-HT${bc}_L${r}_R2.fastq.gz; 
		done; 
done;
