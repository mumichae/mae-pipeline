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
        "bcftools annotate --force -x INFO {input.vcf_file} | bcftools view -s {wildcards.vcf} -m2 -M2 -v snps -O z -o {params.snps_filename}; "
        "bcftools index -t {params.snps_filename}; "
        "{config[gatk]} ASEReadCounter -R {config[genome]} -I {input.bam} -V {params.snps_filename} {params.chrNames} --disable-sequence-dictionary-validation {config[gatk_sanity_check]} | gzip > {output.counted}"
        
        

        
rule test_allelic_counts: 
    input:
        vcf_file="/s/project/sra-download/files/dbGaP-11206/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/temp_gtex",
        bam="/s/project/sra-download/bamfiles/SRR1314810.bam"
    params:
        snps_filename="/s/project/gtex_genetic_diagnosis/processed_data/mae/snps/test_GTEX-13RTJ--SRR1314810.vcf.gz",
        chrNames=" ".join(expand("-L {chr}", chr=config["chr_names"]))
    output:    
        counted="/s/project/gtex_genetic_diagnosis/processed_data/mae/allelic_counts/test_GTEX-13RTJ--SRR1314810.csv.gz"
    shell:
        "bcftools annotate --force -x INFO {input.vcf_file} | bcftools view -s GTEX-13RTJ -m2 -M2 -v snps -O z -o {params.snps_filename}; "
        "bcftools index -t {params.snps_filename}; "
        "gatk ASEReadCounter -R /s/genomes/human/hg19/fasta/Homo_sapiens_assembly19_GTEx.fasta -I {input.bam} -V {params.snps_filename} {params.chrNames} --disable-sequence-dictionary-validation {config[gatk_sanity_check]} | gzip > {output.counted}"




