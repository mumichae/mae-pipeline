#!/bin/bash

# 1 {input.ncbi2ucsc}
# 2 {input.ucsc2ncbi}
# 3 {input.vcf_file}
# 4 {wildcards.vcf}
# 5 {input.bam_file}
# 6 {wildcards.rna}
# 7 {config[mae][genome]}
# 8 {config[mae][gatkIgnoreHeaderCheck]}
# 9 {output.counted}

ncbi2ucsc=$1
ucsc2ncbi=$2
vcf_file=$3
vcf_id=$4
bam_file=$5
rna_id=$6
fasta=$7
sanity=$8
output=$9

tmp=$(mktemp)
header="contig\tposition\tvariantID\trefAllele\taltAllele\t"
header+="refCount\taltCount\ttotalCount\tlowMAPQDepth\t"
header+="lowBaseQDepth\trawDepth\totherBases\timproperPairs"
echo -e $header >> $tmp

# get chr format
vcf_chr=$(bcftools view ${vcf_file} | cut -f1 | grep -v '#' | uniq)
if [ $(echo ${vcf_chr} | grep 'chr' | wc -l) -eq 0 ]
then
    echo "use NCBI format"
    canonical=$ncbi2ucsc
else
    echo "use UCSC format"
    canonical=$ucsc2ncbi
fi
# subset from canonical chromosomes
chr_subset=$(comm -12  <(cut -f1 -d" " ${canonical} | sort -u) <(echo ${vcf_chr} | xargs -n1 | sort -u))
#chr_subset=$(echo $chr_subset | tr ' ' '\n' | sed -e 's/^/-L /' | tr '\n' ' ')

for chr in $chr_subset
do
    echo $chr
    gatk ASEReadCounter \
    -R ${fasta} \
    -I ${bam_file} \
    -V ${vcf_file} \
    -L ${chr} \
    --verbosity ERROR \
    --disable-sequence-dictionary-validation ${sanity} \
    | tail -n+2 >> $tmp
done

cat $tmp | awk -v vcfrna="${vcf_id}--${rna_id}" \
    -F $'\t' 'BEGIN {OFS = FS} NR==1{print $0, "MAE_ID"} NR>1{print $0, vcfrna}' \
    | bgzip > ${output}

rm ${tmp}

