### SNAKEFILE MONOALLELIC EXPRESSION
import os
import drop
import pathlib

## ADD tmp/ DIR
tmpdir = config["ROOT"] + '/' + config["DATASET_NAME"] + '/tmp'
config["tmpdir"] = tmpdir
if not os.path.exists(tmpdir+'/MAE'):
    os.makedirs(tmpdir+'/MAE')

parser = drop.config(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile
htmlOutputPath = config["htmlOutputPath"]
include: os.getcwd() + "/.wBuild/wBuild.snakefile" 


rule all:
    input: rules.Index.output
    output: touch(tmpdir + "/MAE.done")

# create folders for mae results for rule allelic counts
dirs = [parser.getProcDataDir() + "/mae/snps", parser.getProcDataDir() + "/mae/allelic_counts"]
for dir in dirs:
    if not os.path.exists(dir):
        os.makedirs(dir)
        print("Created directory for MAE results: ", dir)

rule create_SNVs:
    input:
        vcf_file=lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, assay=['WES_ASSAY', 'WGS_ASSAY'])
    output:
        snps_filename=parser.getProcDataDir() + "/mae/snps/{vcf}.vcf.gz",
    shell:
        "bcftools annotate --force -x INFO {input.vcf_file} | bcftools view -s {wildcards.vcf} -m2 -M2 -v snps -O z -o {output.snps_filename}; "
        "bcftools index -t {output.snps_filename}; "


rule allelic_counts: 
    input:
        snps_filename=parser.getProcDataDir() + "/mae/snps/{vcf}.vcf.gz",
        bam=lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay='RNA_ASSAY')
    params:
        chrNames=" ".join(expand("-L {chr}", chr=config["chr_names"])),
        script=pathlib.Path(drop.__file__).parent / "modules/mae-pipeline/Scripts/MAE/ASEReadCounter.sh"
    output:    
        counted=parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz"
    shell:
        "{params.script} {config[gatk]} {config[genome]} {input.bam} {input.snps_filename} {config[gatk_sanity_check]} {wildcards.vcf} {wildcards.rna} {output.counted} {params.chrNames}"

rule allelic_counts_qc: 
    input:
        snps_filename=config["qc_vcf"],
        bam=lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay='RNA_ASSAY')
    params:
        chrNames=" ".join(expand("-L {chr}", chr=config["chr_names"]))
    output:    
        counted=parser.getProcDataDir() + "/mae/allelic_counts_qc/{rna}.csv.gz"
    shell:
        "{config[gatk]} ASEReadCounter -R {config[genome]} -I {input.bam} -V {input.snps_filename} {params.chrNames} --disable-sequence-dictionary-validation {config[gatk_sanity_check]} | gzip > {output.counted}"


### RULEGRAPH  
### rulegraph only works without print statements. Call <snakemake produce_graphs> for producing output

## For rule rulegraph.. copy configfile in tmp file
import oyaml
with open(tmpdir + '/config.yaml', 'w') as yaml_file:
    oyaml.dump(config, yaml_file, default_flow_style=False)

rulegraph_filename = htmlOutputPath + "/" + os.path.basename(os.getcwd()) + "_rulegraph"
dag_filename = htmlOutputPath + "/" + os.path.basename(os.getcwd()) + "_dag"

rule produce_graphs:
    input:
        expand("{graph}.{fmt}", fmt=["svg", "png"], graph=[rulegraph_filename, dag_filename])

rule create_rulegraph:
    output:
        rulegraph_filename + ".dot"
    shell:
        "snakemake --configfile " + tmpdir + "/config.yaml --rulegraph > {output}"
        
        
rule create_dag:
    output:
        dag_filename + ".dot"
    shell:
        "snakemake --configfile " + tmpdir + "/config.yaml --dag > {output}"


rule render_dot:
    input:
        "{prefix}.dot"
    output:
        "{prefix}.{fmt,(png|svg)}"
    shell:
        "dot -T{wildcards.fmt} < {input} > {output}"
        
        



