### SNAKEFILE MONOALLELIC EXPRESSION
import os
import drop
import pathlib

tmpdir = os.path.join(config["ROOT"], 'tmp')
config["tmpdir"] = tmpdir
if not os.path.exists(tmpdir+'/MAE'):
    os.makedirs(tmpdir+'/MAE')

parser = drop.config(config)
config = parser.config

include: config['wBuildPath'] + "/wBuild.snakefile"

rule all:
    input: 
        rules.Index.output
    output: touch(tmpdir + "/MAE.done")

# create folders for mae results for rule allelic counts
dirs = [parser.getProcDataDir() + "/mae/snvs", parser.getProcDataDir() + "/mae/allelic_counts"]
for dir in dirs:
    if not os.path.exists(dir):
        os.makedirs(dir)
        print("Created directory for MAE results: ", dir)

MAE_ROOT = pathlib.Path(drop.__file__).parent / "modules/mae-pipeline"
rule create_SNVs:
    input:
        vcf_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, assay=['WES_ASSAY', 'WGS_ASSAY']),
        bam_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay='RNA_ASSAY'),
        chr_conv1 = MAE_ROOT / "Scripts/MAE/chr_NCBI_UCSC.txt", 
        chr_conv2 = MAE_ROOT / "Scripts/MAE/chr_UCSC_NCBI.txt",
        script = MAE_ROOT / "Scripts/MAE/filterSNVs.sh"
    output:
        snvs_filename=parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz",
        snvs_index=parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz.tbi"
    shell:
        "{input.script} {input.vcf_file} {wildcards.vcf} {input.bam_file} {output.snvs_filename}"

rule allelic_counts: 
    input:
        snps_filename=parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz",
        bam_file=lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay='RNA_ASSAY'),
        script=pathlib.Path(drop.__file__).parent / "modules/mae-pipeline/Scripts/MAE/ASEReadCounter.sh"
    output:    
        counted=parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz"
    shell:
        "{input.script} {config[gatk]} {config[genome]} {input.bam_file} {input.snps_filename} {config[gatk_sanity_check]} {wildcards.vcf} {wildcards.rna} {output.counted}"

rule allelic_counts_qc: 
    input:
        snps_filename=config["qc_vcf"],
        bam_file=lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay='RNA_ASSAY'),
        script=pathlib.Path(drop.__file__).parent / "modules/mae-pipeline/Scripts/MAE/ASEReadCounter.sh"
    output:    
        counted=parser.getProcDataDir() + "/mae/allelic_counts_qc/{rna}.csv.gz"
    shell:
        "{input.script} {config[gatk]} {config[genome]} {input.bam_file} {input.snps_filename} {config[gatk_sanity_check]} qc {wildcards.rna} {output.counted}"

