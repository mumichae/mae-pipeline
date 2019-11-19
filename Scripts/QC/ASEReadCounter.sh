#!/bin/bash

# 1 {input.ncbi2ucsc}
# 2 {input.ucsc2ncbi}
# 3 {input.vcf_file_ucsc}
# 4 {input.vcf_file_ncbi}
# 5 {input.bam_file}
# 6 {wildcards.rna}
# 7 {config[mae][genome]}
# 8 {config[mae][gatkIgnoreHeaderCheck]}
# 9 {output.counted}
# 10 {config[tools][bcftoolsCmd]}
# 11 {config[tools][samtoolsCmd]}

ncbi2ucsc=$1
ucsc2ncbi=$2
vcf_file_ucsc=$3
vcf_file_ncbi=$4
bam_file=$5
rna_id=$6
fasta=$7
sanity=$8
output=$9
bcftools=${10}
samtools=${11}

tmp=$(mktemp)
header="contig\tposition\tvariantID\trefAllele\taltAllele\t"
header+="refCount\taltCount\ttotalCount\tlowMAPQDepth\t"
header+="lowBaseQDepth\trawDepth\totherBases\timproperPairs"
echo -e $header >> $tmp


# get number of UCSC chromosomes in BAM
bam_chr=$($samtools idxstats ${bam_file} | grep chr | wc -l)
if [ ${bam_chr} -ne 0 ]
then
    echo "use UCSC format"
    canonical=$ucsc2ncbi
    vcf_file=$vcf_file_ucsc
else
    echo "use NCBI format"
    canonical=$ncbi2ucsc
    vcf_file=$vcf_file_ncbi
fi

# get unique chromosomes
vcf_chr=$($bcftools view ${vcf_file} | cut -f1 | grep -v '#' | uniq)
# subset from canonical chromosomes
chr_subset=$(comm -12  <(cut -f1 -d" " ${canonical} | sort -u) <(echo ${vcf_chr} | xargs -n1 | sort -u))

for chr in $chr_subset
do
    echo $chr1
    gatk ASEReadCounter \
    -R ${fasta} \
    -I ${bam_file} \
    -V ${vcf_file} \
    -L ${chr} \
    --verbosity ERROR \
    --disable-sequence-dictionary-validation ${sanity} \
    | tail -n+2 >> $tmp
done

cat $tmp | awk -v rnaID="${rna_id}" \
    -F $'\t' 'BEGIN {OFS = FS} NR==1{print $0, "RNA_ID"} NR>1{print $0, rnaID}' \
    | bgzip > ${output}

rm ${tmp}


