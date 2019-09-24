#!/bin/bash

gatk=$1 # {config[gatk]} # path to gatk
fasta=$2 # {config[genome]} # genome fasta file
bam_file=$3 # {input.bam}
vcf_file=$4 # {input.snps_filename}
sanity=$5 # {config[gatk_sanity_check]}
vcf_id=$6 # {vcf_id}
rna_id=$7 # {rna_id}
output=$8 # {output.counted}

# get chr format of BAM
bam_chr=$(samtools idxstats ${bam_file} | grep chr | wc -l)
# compare chr format and rename vcf if necessary
if [ ${bam_chr} -eq 0 ]
then 
    chr_mod="$(dirname $0)/chr_NCBI_UCSC.txt"
else
    chr_mod="$(dirname $0)/chr_UCSC_NCBI.txt"
fi

chrs=$(cut -f1 $chr_mod)
echo $chrs

#$gatk ASEReadCounter -R $fasta -I ${bam_file} -V ${vcf_file} $chrs --disable-sequence-dictionary-validation ${sanity} | \
#awk -v vcfrna="${vcf_id}--${rna_id}" -F $'\t' 'BEGIN {OFS = FS} NR==1{print $0, "mae_id"} NR>1{print $0, vcfrna}' | gzip > ${output}
