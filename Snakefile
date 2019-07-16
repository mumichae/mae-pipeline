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


# create folders for mae results
dirs = [parser.getProcDataDir() + "/mae/snps", parser.getProcDataDir() + "/mae/allelic_counts"]
for dir in dirs:
    if not os.path.exists(dir):
        os.makedirs(dir)
        print("Created directory for MAE results: ", dir)
#print(parser.getMaeIDs())       
rule allelic_counts: 
    input:
        vcf_file=lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, isRNA=False),
        bam=lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, isRNA=True)
    params:
        snps_filename=parser.getProcDataDir() + "/mae/snps/{vcf}--{rna}.vcf.gz",
        chrNames=" ".join(expand("-L {chr}", chr=config["chr_names"]))
    output:    
        counted=parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz"
    shell:
        "bcftools annotate -O b -x INFO {input.vcf_file} | bcftools view -s {wildcards.vcf} -m2 -M2 -v snps -O z > {params.snps_filename}; "
        "bcftools index -t {params.snps_filename}; "
        "{config[gatk]} ASEReadCounter -R {config[genome]} -I {input.bam} -V {params.snps_filename} {params.chrNames} --disable-sequence-dictionary-validation {config[gatk_sanity_check]} | gzip > {output.counted}"
        
        

        
rule test_allelic_counts: 
    input:
        vcf="/s/project/mitoMultiOmics/raw_data/helmholtz/33281/exomicout/paired-endout/stdFilenames/33281.vcf.gz",
        bam="/s/project/mitoMultiOmics/raw_data/helmholtz/MUC1343/RNAout/paired-endout/stdFilenames/MUC1343.bam"
    params:
        snps_filename="tmp/33281--MUC1343.vcf.gz",
        chrNames=" ".join(expand("-L {chr}", chr=config["chr_names"]))
    output:    
        counted="tmp/counted_33281--MUC1343.csv.gz"
    shell:
        "bcftools annotate -O b -x INFO /s/project/mitoMultiOmics/raw_data/helmholtz/33281/exomicout/paired-endout/stdFilenames/33281.vcf.gz | bcftools view -s 33281  -m2 -M2 -v snps -O z > 33281--MUC1343.vcf.gz; "
        "bcftools index -t 33281--MUC1343.vcf.gz; "
        "gatk ASEReadCounter -R /s/genomes/human/hg19/fasta/hg19.fa -I /s/project/mitoMultiOmics/raw_data/helmholtz/MUC1343/RNAout/paired-endout/stdFilenames/MUC1343.bam -V 33281--MUC1343.vcf.gz -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY | gzip > counted_33281--MUC1343.csv.gz"




