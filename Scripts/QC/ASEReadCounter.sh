#!/bin/bash

# 1 {input.ncbi2ucsc}
# 2 {input.ucsc2ncbi}
# 3 {input.vcf_file_ucsc}
# 4 {input.vcf_file_ncbi}
# 5 {input.bam_file}
# 6 {wildcards.rna}
# 7 {config[genome]}
# 8 {config[gatk_sanity_check]}
# 9 {output.counted}

ncbi2ucsc=$1
ucsc2ncbi=$2
vcf_file_ucsc=$3
vcf_file_ncbi=$4
bam_file=$5
rna_id=$6
fasta=$7
sanity=$8
output=$9

echo $ncbi2ucsc

# get number of UCSC chromosomes in BAM
bam_chr=$(samtools idxstats ${bam_file} | grep chr | wc -l)

if [ ${bam_chr} -eq 0 ]
then
    echo "use UCSC format"
    canonical=$ncbi2ucsc
    vcf_file=$vcf_file_ucsc
else
    echo "use NCBI format"
    canonical=$ucsc2ncbi
    vcf_file=$vcf_file_ncbi
fi

# get unique chromosomes
vcf_chr=$(bcftools view ${vcf_file} | cut -f1 | grep -v '#' | uniq)
# subset from canonical chromosomes
chr_subset=$(comm -1  <(cut -f1 -d" " ${canonical} | sort) <(echo ${vcf_chr} | sort))
chr_subset=$(echo $chr_subset | tr ' ' '\n' | sed -e 's/^/-L /' | tr '\n' ' ')

gatk ASEReadCounter \
    -R $fasta \
    -I ${bam_file} \
    -V ${vcf_file} \
    ${chr_subset} \
    --verbosity ERROR \
    --disable-sequence-dictionary-validation ${sanity} \
   | awk -v rnaID="${rna_id}" \
    -F $'\t' 'BEGIN {OFS = FS} NR==1{print $0, "RNA_ID"} NR>1{print $0, rnaID}' \
    | bgzip > ${output}

