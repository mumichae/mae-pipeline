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
MAE_ROOT = pathlib.Path(drop.__file__).parent / "modules/mae-pipeline"

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

rule create_SNVs:
    input:
        ncbi2ucsc = MAE_ROOT / "resource/chr_NCBI_UCSC.txt",
        ucsc2ncbi = MAE_ROOT / "resource/chr_UCSC_NCBI.txt",
        vcf_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, assay=['WES_ASSAY', 'WGS_ASSAY']),
        bam_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay='RNA_ASSAY'),
        script = MAE_ROOT / "Scripts/MAE/filterSNVs.sh"
    output:
        snvs_filename=parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz",
        snvs_index=parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz.tbi"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} {input.vcf_file} {wildcards.vcf} \
        {input.bam_file} {output.snvs_filename}
        """

rule allelic_counts: 
    input:
        ncbi2ucsc = MAE_ROOT / "resource/chr_NCBI_UCSC.txt",
        ucsc2ncbi = MAE_ROOT / "resource/chr_UCSC_NCBI.txt",
        vcf_file = parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz",
        bam_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay='RNA_ASSAY'),
        script = MAE_ROOT / "Scripts/MAE/ASEReadCounter.sh"
    output:    
        counted = parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file} {wildcards.vcf} {input.bam_file} {wildcards.rna} \
        {config[genome]} {config[gatk_sanity_check]} {output.counted}
        """

rule allelic_counts_qc: 
    input:
        ncbi2ucsc = MAE_ROOT / "resource/chr_NCBI_UCSC.txt",
        ucsc2ncbi = MAE_ROOT / "resource/chr_UCSC_NCBI.txt",
        vcf_file_ucsc = config["qc_vcf"]["UCSC"],
        vcf_file_ncbi = config["qc_vcf"]["NCBI"],
        bam_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay='RNA_ASSAY'),
        script = MAE_ROOT / "Scripts/QC/ASEReadCounter.sh"
    output:    
        counted = parser.getProcDataDir() + "/mae/allelic_counts_qc/{rna}.csv.gz"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file_ucsc} {input.vcf_file_ncbi} {input.bam_file} {wildcards.rna} \
        {config[genome]} {config[gatk_sanity_check]} {output.counted}
        """

