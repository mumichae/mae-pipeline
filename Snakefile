### SNAKEFILE MONOALLELIC EXPRESSION
import os
import drop
import pathlib

parser = drop.config(config)
config = parser.config
MAE_ROOT = pathlib.Path(drop.__file__).parent / "modules/mae-pipeline"

include: config['wBuildPath'] + "/wBuild.snakefile"

rule all:
    input: 
        rules.Index.output, 
        expand(
            parser.getProcResultsDir() + "/mae/{dataset}/MAE_results_{annotation}.tsv",
            dataset=parser.mae_ids.keys(), annotation=list(config["geneAnnotation"].keys())
        ),
        parser.getProcResultsDir() + "/mae/" + config["mae"]["qcGroup"] + "/dna_rna_qc_matrix.Rds"
    output: touch(tmpdir + "/MAE.done")

# create folders for mae results for rule allelic counts
dirs = [parser.getProcDataDir() + "/mae/snvs", parser.getProcDataDir() + "/mae/allelic_counts"]
for dir_ in dirs:
    if not os.path.exists(dir_):
        os.makedirs(dirdir_)
        print("Created directory for MAE results: ", dir_)

rule create_SNVs:
    input:
        ncbi2ucsc = MAE_ROOT / "resource/chr_NCBI_UCSC.txt",
        ucsc2ncbi = MAE_ROOT / "resource/chr_UCSC_NCBI.txt",
        vcf_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, file_type='DNA_VCF_FILE'),
        bam_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, file_type='RNA_BAM_FILE'),
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
        bam_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, file_type='RNA_BAM_FILE'),
        script = MAE_ROOT / "Scripts/MAE/ASEReadCounter.sh"
    output:    
        counted = parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file} {wildcards.vcf} {input.bam_file} {wildcards.rna} \
        {config[mae][genome]} {config[mae][gatkIgnoreHeaderCheck]} {output.counted}
        """

rule allelic_counts_qc: 
    input:
        ncbi2ucsc = MAE_ROOT / "resource/chr_NCBI_UCSC.txt",
        ucsc2ncbi = MAE_ROOT / "resource/chr_UCSC_NCBI.txt",
        vcf_file_ucsc = config["mae"]["qcVcf"]["UCSC"],
        vcf_file_ncbi = config["mae"]["qcVcf"]["NCBI"],
        bam_file = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, file_type='RNA_BAM_FILE'),
        script = MAE_ROOT / "Scripts/QC/ASEReadCounter.sh"
    output:    
        counted = parser.getProcDataDir() + "/mae/allelic_counts_qc/{rna}.csv.gz"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file_ucsc} {input.vcf_file_ncbi} {input.bam_file} {wildcards.rna} \
        {config[mae][genome]} {config[mae][gatkIgnoreHeaderCheck]} {output.counted}
        """


### RULEGRAPH  
### rulegraph only works without print statements

## For rule rulegraph.. copy configfile in tmp file
import oyaml
with open(tmpdir + '/config.yaml', 'w') as yaml_file:
    oyaml.dump(config, yaml_file, default_flow_style=False)

rulegraph_filename = config["htmlOutputPath"] + "/MAE_rulegraph" # htmlOutputPath + "/" + os.path.basename(os.getcwd()) + "_rulegraph"
rule produce_rulegraph:
    input:
        expand(rulegraph_filename + ".{fmt}", fmt=["svg", "png"])

rule create_graph:
    output:
        rulegraph_filename + ".dot"
    shell:
        "snakemake --configfile " + tmpdir + "/config.yaml --rulegraph > {output}"

rule render_dot:
    input:
        "{prefix}.dot"
    output:
        "{prefix}.{fmt,(png|svg)}"
    shell:
        "dot -T{wildcards.fmt} < {input} > {output}"

