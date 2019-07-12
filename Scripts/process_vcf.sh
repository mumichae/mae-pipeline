#!/bin/bash  

{ zcat /s/project/mitoMultiOmics/raw_data/helmholtz/$1/exomicout/paired-endout/stdFilenames/$1.vcf.gz | grep "^#" ; zcat /s/project/mitoMultiOmics/raw_data/helmholtz/$1/exomicout/paired-endout/stdFilenames/$1.vcf.gz | grep -v "^#" | awk -F"\t" '{$8=".";print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}';} | bgzip -c > /s/project/mitoMultiOmics/raw_data/helmholtz/$1/exomicout/paired-endout/processedData/temp_processed_$1.vcf.gz

tabix /s/project/mitoMultiOmics/raw_data/helmholtz/$1/exomicout/paired-endout/processedData/temp_processed_$1.vcf.gz

Rscript --vanilla /data/ouga/home/ag_gagneur/figueira/workspace/genePROF_tests/Scripts/VCF_Process/extract_invalid_positions.R $1 

gatk SelectVariants -V /s/project/mitoMultiOmics/raw_data/helmholtz/$1/exomicout/paired-endout/processedData/temp_processed_$1.vcf.gz -O /s/project/mitoMultiOmics/raw_data/helmholtz/$1/exomicout/paired-endout/processedData/$1_snps.vcf.gz -restrict-alleles-to BIALLELIC --select-type-to-include SNP -XL /s/project/mitoMultiOmics/raw_data/helmholtz/$1/exomicout/paired-endout/processedData/invalid.intervals

gatk ASEReadCounter -R /s/genomes/human/hg19/fasta/hg19.fa \
-I /s/project/mitoMultiOmics/raw_data/helmholtz/$2/RNAout/paired-endout/stdFilenames/$2.bam \
-V /s/project/mitoMultiOmics/raw_data/helmholtz/$1/exomicout/paired-endout/processedData/$1_snps.vcf.gz \
-L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY | gzip > /s/project/mitoMultiOmics/raw_data/helmholtz/$1/exomicout/paired-endout/processedData/count_ASE_$1.csv.gz

bcftools annotate -O b -x INFO  /s/public_webshare/public/genetic_diagnosis_demo/raw_data/vcfs-variants/SUBSET.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.exome.vcf.gz | bcftools view -s HG00189 | less -S