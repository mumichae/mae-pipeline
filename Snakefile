### SNAKEFILE MONOALLELIC EXPRESSION

import os
from config_parser import ConfigHelper

parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile
htmlOutputPath = config["htmlOutputPath"]

include: os.getcwd() + "/.wBuild/wBuild.snakefile"
# create temporary folder
if not os.path.exists('tmp'):
    os.makedirs('tmp')

rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("Output/all.done")
    
# overwriting wbuild rule output
rule rulegraph:
    shell: "snakemake --rulegraph | dot -Tsvg -Grankdir=TB > {config[htmlOutputPath]}/dep.svg"

rule allelic_counts:
    input:
        vcf=/s/project/mitoMultiOmics/raw_data/helmholtz/85154/exomicout/paired-endout/stdFilenames/85154.vcf.gz
        snps=/s/project/mitoMultiOmics/raw_data/helmholtz/85154/exomicout/paired-endout/processedData/85154_snps.vcf.gz
        genome=/s/genomes/human/hg19/fasta/hg19.fa
        bam=/s/project/mitoMultiOmics/raw_data/helmholtz/MUC1412/RNAout/paired-endout/stdFilenames/MUC1412.bam
        chrNames=" ".join(expand("-L {chr}", chr=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]))
    output:    
        counted=/s/project/mitoMultiOmics/raw_data/helmholtz/85154/exomicout/paired-endout/processedData/count_ASE_85154.csv.gz
    shell:
        "bcftools annotate -O b -x INFO {input.vcf} | bcftools view -s 85154  -m2 -M2 -v snps -O z > {input.snps}"
        "bcftools index -t {input.snps}"
        "gatk ASEReadCounter -R {input.genome} -I {input.bam} -V {input.snps} {input.chrNames} | gzip > {output.counted}"
"
        
        