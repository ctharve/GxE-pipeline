#!/bin/bash

set -v
set -e

cl=$1
minCov="15"
echo ${cl}

less ../output/QuASAR_genotypes_noPileupMin_${cl}_minCov${minCov}.txt | awk '($8)>.5 {print $1,$3,$6}' | while read chr pos pr; do bcftools view -H -s ${cl} --types snps -Ov ../data/1KG_genotypes/ALL.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -r ${chr}:${pos} | awk '{print $10}' | sed 's/:[^\t]*//g' | awk '$1=="0|0" || $1=="1|1"'; done | wc -l
less ../output/QuASAR_genotypes_noPileupMin_${cl}_minCov${minCov}.txt | awk '($8)>.9 {print $1,$3,$6}' | while read chr pos pr; do bcftools view -H -s ${cl} --types snps -Ov ../data/1KG_genotypes/ALL.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -r ${chr}:${pos} | awk '{print $10}' | sed 's/:[^\t]*//g' | awk '$1=="0|0" || $1=="1|1"'; done | wc -l
less ../output/QuASAR_genotypes_noPileupMin_${cl}_minCov${minCov}.txt | awk '($8)>.99 {print $1,$3,$6}' | while read chr pos pr; do bcftools view -H -s ${cl} --types snps -Ov ../data/1KG_genotypes/ALL.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -r ${chr}:${pos} | awk '{print $10}' | sed 's/:[^\t]*//g' | awk '$1=="0|0" || $1=="1|1"'; done | wc -l
less ../output/QuASAR_genotypes_noPileupMin_${cl}_minCov${minCov}.txt | awk '(1-$8)>.5 {print $1,$3,$6}' | while read chr pos pr; do bcftools view -H -s ${cl} --types snps -Ov ../data/1KG_genotypes/ALL.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -r ${chr}:${pos} | awk '{print $10}' | sed 's/:[^\t]*//g' | awk '$1=="0|1" || $1=="1|0"'; done | wc -l
less ../output/QuASAR_genotypes_noPileupMin_${cl}_minCov${minCov}.txt | awk '(1-$8)>.9 {print $1,$3,$6}' | while read chr pos pr; do bcftools view -H -s ${cl} --types snps -Ov ../data/1KG_genotypes/ALL.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -r ${chr}:${pos} | awk '{print $10}' | sed 's/:[^\t]*//g' | awk '$1=="0|1" || $1=="1|0"'; done | wc -l
less ../output/QuASAR_genotypes_noPileupMin_${cl}_minCov${minCov}.txt | awk '(1-$8)>.99 {print $1,$3,$6}' | while read chr pos pr; do bcftools view -H -s ${cl} --types snps -Ov ../data/1KG_genotypes/ALL.chr${chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz -r ${chr}:${pos} | awk '{print $10}' | sed 's/:[^\t]*//g' | awk '$1=="0|1" || $1=="1|0"'; done | wc -l
