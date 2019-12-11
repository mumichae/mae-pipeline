### SNAKEFILE MONOALLELIC EXPRESSION
import os
import drop
import pathlib

METHOD = 'MAE'
SCRIPT_ROOT = os.getcwd() #drop.getMethodPath(METHOD, type_='workdir')
CONF_FILE = drop.getConfFile()

parser = drop.config(config, METHOD)
config = parser.parse()
include: config['wBuildPath'] + "/wBuild.snakefile"

rule all:
    input: 
        rules.Index.output, 
        expand(
            parser.getProcResultsDir() + "/mae/{dataset}/MAE_results_{annotation}.tsv",
            dataset=parser.mae_ids.keys(), annotation=list(config["geneAnnotation"].keys())
        ),
        parser.getProcResultsDir() + "/mae/" + config["mae"]["qcGroup"] + 
        "/dna_rna_qc_matrix.Rds"
    output: touch(drop.getMethodPath(METHOD, type_='final_file'))

rule sampleQC:
    input: rules.Scripts_QC_DNA_RNA_matrix_plot_R.output
    output: drop.getTmpDir() + "/sampleQC.done"

def fasta_dict(fasta_file):
    return ".".join(fasta_file.split('.')[0:-1]) + ".dict"

rule create_dict:
    input: config['mae']['genome']
    output: fasta_dict(config['mae']['genome'])
    shell: "gatk CreateSequenceDictionary --REFERENCE {input[0]}"
        

rule create_SNVs:
    input:
        ncbi2ucsc = os.path.join(SCRIPT_ROOT, "resource", "chr_NCBI_UCSC.txt"),
        ucsc2ncbi = os.path.join(SCRIPT_ROOT, "resource", "chr_UCSC_NCBI.txt"),
        vcf_file  = lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, 
                    file_type='DNA_VCF_FILE'),
        bam_file  = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, 
                    file_type='RNA_BAM_FILE'),
        script    = os.path.join(SCRIPT_ROOT, "Scripts", "MAE", "filterSNVs.sh")
    output:
        snvs_filename=parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz",
        snvs_index=parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz.tbi"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} {input.vcf_file} \
        {wildcards.vcf} {input.bam_file} {output.snvs_filename} \
        {config[tools][bcftoolsCmd]} {config[tools][samtoolsCmd]}
        """

rule allelic_counts: 
    input:
        ncbi2ucsc = os.path.join(SCRIPT_ROOT, "resource", "chr_NCBI_UCSC.txt"),
        ucsc2ncbi = os.path.join(SCRIPT_ROOT, "resource", "chr_UCSC_NCBI.txt"),
        vcf_file  = parser.getProcDataDir() + "/mae/snvs/{vcf}--{rna}.vcf.gz",
        bam_file  = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, 
                    file_type='RNA_BAM_FILE'),
        fasta     = config['mae']['genome'],
        dict      = fasta_dict(config['mae']['genome']),
        script    = os.path.join(SCRIPT_ROOT, "Scripts", "MAE", "ASEReadCounter.sh")
    output:    
        counted = parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file} {wildcards.vcf} {input.bam_file} {wildcards.rna} \
        {input.fasta} {config[mae][gatkIgnoreHeaderCheck]} {output.counted} \
        {config[tools][bcftoolsCmd]}
        """

rule allelic_counts_qc: 
    input:
        ncbi2ucsc     = os.path.join(SCRIPT_ROOT, "resource", "chr_NCBI_UCSC.txt"),
        ucsc2ncbi     = os.path.join(SCRIPT_ROOT, "resource", "chr_UCSC_NCBI.txt"),
        vcf_file_ucsc = config["mae"]["qcVcf"]["UCSC"],
        vcf_file_ncbi = config["mae"]["qcVcf"]["NCBI"],
        bam_file      = lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, 
                        file_type='RNA_BAM_FILE'),
        fasta         = config['mae']['genome'],
        dict          = fasta_dict(config['mae']['genome']),
        script        = os.path.join(SCRIPT_ROOT, "Scripts", "QC", "ASEReadCounter.sh")
    output:    
        counted = parser.getProcDataDir() + "/mae/allelic_counts_qc/{rna}.csv.gz"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file_ucsc} {input.vcf_file_ncbi} {input.bam_file} \
        {wildcards.rna} {input.fasta} {config[mae][gatkIgnoreHeaderCheck]} \
        {output.counted} {config[tools][bcftoolsCmd]} {config[tools][samtoolsCmd]}
        """

rulegraph_filename = f'{config["htmlOutputPath"]}/{METHOD}_rulegraph'

rule produce_rulegraph:
    input:
        expand(rulegraph_filename + ".{fmt}", fmt=["svg", "png"])

rule create_graph:
    output:
        svg = f"{rulegraph_filename}.svg",
        png = f"{rulegraph_filename}.png"
    shell:
        """
        snakemake --configfile {CONF_FILE} --rulegraph | dot -Tsvg > {output.svg}
        snakemake --configfile {CONF_FILE} --rulegraph | dot -Tpng > {output.png}
        """

rule unlock:
    output: touch(drop.getMethodPath(METHOD, type_="unlock"))
    shell: "snakemake --unlock --configfile {CONF_FILE}"

