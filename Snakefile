### SNAKEFILE MONOALLELIC EXPRESSION

import os
from config_parser import ConfigHelper

#print(config)
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
        vcf_file=lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, isRNA=False),
        bam=lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, isRNA=True)
    params:
        snps_filename=parser.getProcDataDir() + "/mae/snps/{vcf}--{rna}.vcf.gz",
        chrNames=" ".join(expand("-L {chr}", chr=config["chr_names"]))
    output:    
        counted=parser.getProcDataDir() + "/mae/test/{vcf}--{rna}.csv.gz"
    shell:
        "bcftools annotate -x FMT/GL -k -i'N_ALT==1' {input.vcf_file} | bcftools annotate -O b -x INFO | bcftools view -s {wildcards.vcf} -m2 -M2 -v snps -O z > {params.snps_filename}; "  #bcftools annotate -x FMT/GL -k -i'N_ALT==1' file.vcf
        "bcftools index -t {params.snps_filename}; "
        "{config[gatk]} ASEReadCounter -R {config[genome]} -I {input.bam} -V {params.snps_filename} {params.chrNames} --disable-sequence-dictionary-validation {config[gatk_sanity_check]} | gzip > {output.counted}"
        

rule test:
    input: parser.getProcDataDir() + "/mae/test/33254--33254R_GAL.csv.gz" #/s/project/prokisch/processed_data/mae/allelic_counts/33254--33254R_GAL.csv.gz



