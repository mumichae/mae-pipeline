#!/bin/bash

vcf_file=$1 # {input.vcf_file}
vcf_id=$2 # {wildcards.vcf}
bam_file=$3
output=$4 # {output.snvs_filename}

#echo "vcf file:" $vcf_file
#echo "vcf ID:" $vcf_id
#echo "bam file:" $bam_file
#echo "output:" $output

tmp=${output}_tmp

bcftools annotate --force -x INFO ${vcf_file} |\
    bcftools view -s ${vcf_id} -m2 -M2 -v snps -O z -o $tmp
bcftools index -t $tmp
echo "wrote" $tmp

# compare and correct chromosome format mismatch
bam_chr=$(samtools idxstats ${bam_file} | grep chr | wc -l)
vcf_chr=$(bcftools view $tmp | cut -f1 | grep -v '#' | grep chr | wc -l)

if [ ${vcf_chr} -eq 0  ] && [ ${bam_chr} -ne 0 ]
then # rename vcf to UCSC
    chr_mod="$(dirname $0)/chr_NCBI_UCSC.txt"
    echo "converting from NCBI to UCSC format"
    bcftools annotate --rename-chrs ${chr_mod} $tmp | bgzip > ${output}
    rm ${tmp}
    rm ${tmp}.tbi
elif [ ${vcf_chr} -ne 0  ] && [ ${bam_chr} -eq 0 ]
then # rename vcf to NCBI
    chr_mod="$(dirname $0)/chr_UCSC_NCBI.txt"
    echo "converting from UCSC to NCBI format"
    bcftools annotate --rename-chrs ${chr_mod} $tmp | bgzip > ${output}
    rm ${tmp}
    rm ${tmp}.tbi
else
    mv $tmp ${output}
    rm ${tmp}.tbi
fi

bcftools index -t ${output}
