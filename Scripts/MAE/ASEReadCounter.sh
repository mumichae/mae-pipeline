#!/bin/bash

# $1 {config[gatk]}
# $2 {config[genome]}
# $3 {input.bam}
# $4 {input.snps_filename}
# $5 {config[gatk_sanity_check]}
# $6 {vcf}
# $7 {rna}
# $8 {output.counted}
# $9 {params.chrNames}

separator="--"
chrs=${@:9}

$1 ASEReadCounter -R $2 -I $3 -V $4 $chrs --disable-sequence-dictionary-validation $5 | awk -v vcfrna="$6$separator$7" -F $'\t' 'BEGIN {OFS = FS} NR==1{print $0, "mae_id"} NR>1{print $0, vcfrna}' | gzip > $8
