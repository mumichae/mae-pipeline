#!/bin/bash

gatk=$1 # {config[gatk]} # path to gatk
fasta=$2 # {config[genome]} # genome fasta file
bam_file=$3 # {input.bam}
vcf_file=$4 # {input.snps_filename}
sanity=$5 # {config[gatk_sanity_check]}
vcf_id=$6 # {vcf_id}
rna_id=$7 # {rna_id}
output=$8 # {output.counted}

# get chr format
vcf_chr=$(bcftools view ${vcf_file} | cut -f1 | grep -v '#' | uniq)

if [ $(echo ${vcf_chr} | grep 'chr' | wc -l) -eq 0 ]
then
    echo "use NCBI format"
    canonical="$(dirname $0)/chr_NCBI_UCSC.txt"
else
    echo "use UCSC format"
    canonical="$(dirname $0)/chr_UCSC_NCBI.txt"
fi

# subset from canonical chromosomes
chr_subset=$(comm -1  <(cut -f1 -d" " ${canonical} | sort) <(echo ${vcf_chr} | sort))
chr_subset=$(echo $chr_subset | tr ' ' '\n' | sed -e 's/^/-L /' | tr '\n' ' ')

$gatk ASEReadCounter \
    -R $fasta \
    -I ${bam_file} \
    -V ${vcf_file} \
    ${chr_subset} \
    --verbosity ERROR \
    --disable-sequence-dictionary-validation ${sanity} \
   | awk -v vcfrna="${vcf_id}--${rna_id}" \
    -F $'\t' 'BEGIN {OFS = FS} NR==1{print $0, "MAE_ID"} NR>1{print $0, vcfrna}' \
    | bgzip > ${output}

