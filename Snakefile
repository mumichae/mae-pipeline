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
        
rule allelic_counts: 
    input:
        vcf_file=lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, assay_name="wes_assay"),
        bam=lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay_name="rna_assay")
    params:
        snps_filename=parser.getProcDataDir() + "/mae/snps/{vcf}--{rna}.vcf.gz",
        chrNames=" ".join(expand("-L {chr}", chr=config["chr_names"]))
    output:    
        counted=parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz"
    shell:
        "bcftools annotate -O b -x INFO {input.vcf_file} | bcftools view -s {wildcards.vcf} -m2 -M2 -v snps -O z > {params.snps_filename}; "
        "bcftools index -t {params.snps_filename}; "
        "gatk ASEReadCounter -R {config[genome]} -I {input.bam} -V {params.snps_filename} {params.chrNames} | gzip > {output.counted}"
        
        
        
        
rule test_allelic_counts: 
    input:
        vcf="/s/project/mitoMultiOmics/raw_data/helmholtz/85154/exomicout/paired-endout/stdFilenames/85154.vcf.gz",
        bam="/s/project/mitoMultiOmics/raw_data/helmholtz/MUC1412/RNAout/paired-endout/stdFilenames/MUC1412.bam"
    params:
        snps_filename="/s/project/mitoMultiOmics/raw_data/helmholtz/85154/exomicout/paired-endout/processedData/v2_85154_snps.vcf.gz",
        chrNames=" ".join(expand("-L {chr}", chr=config["chr_names"]))
    output:    
        counted="/s/project/mitoMultiOmics/raw_data/helmholtz/85154/exomicout/paired-endout/processedData/cv2_ount_ASE_85154.csv.gz"
    shell:
        "bcftools annotate -O b -x INFO {input.vcf} | bcftools view -s 85154  -m2 -M2 -v snps -O z > {params.snps_filename}; "
        "bcftools index -t {params.snps_filename}; "
        "gatk ASEReadCounter -R {config[genome]} -I {input.bam} -V {params.snps_filename} {params.chrNames} | gzip > {output.counted}"



