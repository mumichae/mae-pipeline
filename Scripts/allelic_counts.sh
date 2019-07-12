#!/bin/bash
chromosomeNames=$6

bcftools annotate -O b -x INFO $1 | bcftools view -s 85154  -m2 -M2 -v snps -O z > $2
bcftools index -t $2
gatk ASEReadCounter -R $3 \
-I $4 \
-V $2 \
-L "$chromosomeNames" | gzip > $5

### List of arguments
# 1: vcf
# 2: name of vcf_processed: snps
# 3: genome
# 4: bam
# 5: output: counted
#6: List of chromosome names 
 
 