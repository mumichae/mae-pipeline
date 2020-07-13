WORKDIR = cfg.MAE.getWorkdir(str_=False)

###### FUNCTIONS ######
def fasta_dict(fasta_file):
    return fasta_file.split('.')[0] + ".dict"

def getVcf(rna_id, vcf_id="qc"):
    if vcf_id == "qc":
        return config["mae"]["qcVcf"]
    else:
        return cfg.getProcessedDataDir() + f"/mae/snvs/{vcf_id}--{rna_id}.vcf.gz"
        
def getQC(format):
    if format == "UCSC":
        return config["mae"]["qcVcf"]
    elif format == "NCBI":
        return cfg.processedDataDir / "mae" / "qc_vcf_ncbi.vcf.gz"
    else:
        raise ValueError(f"getQC: {format} is an invalid chromosome format")

def getChrMap(WORKDIR , conversion):
    if conversion == 'ncbi2ucsc':
        return WORKDIR /"resource"/"chr_NCBI_UCSC.txt"
    elif conversion == 'ucsc2ncbi':
        return WORKDIR /"resource"/"chr_UCSC_NCBI.txt"
    else:
        raise ValueError(f"getChrMap: {conversion} is an invalid conversion option")
        
######

rule mae:
    input:
        cfg.getHtmlFromScript(WORKDIR / "MAE" / "Datasets.R"),
        cfg.getHtmlFromScript(WORKDIR / "QC" / "Datasets.R")

rule sampleQC:
    input: cfg.getHtmlFromScript(WORKDIR / "QC" / "Datasets.R")

rule create_dict:
    input: config['mae']['genome']
    output: fasta_dict(config['mae']['genome'])
    shell: "gatk CreateSequenceDictionary --REFERENCE {input[0]}"
        
## MAE
rule mae_createSNVs:
    input:
        ncbi2ucsc = getChrMap(WORKDIR, "ncbi2ucsc"),
        ucsc2ncbi = getChrMap(WORKDIR, "ucsc2ncbi"),
        vcf_file  = lambda w: sa.getFilePath(w.vcf, 'DNA_VCF_FILE'),
        bam_file  = lambda w: sa.getFilePath(w.rna, 'RNA_BAM_FILE'),
        script    = WORKDIR / "MAE" / "filterSNVs.sh"
    output:
        snvs_filename = cfg.processedDataDir / "mae" / "snvs" / "{vcf}--{rna}.vcf.gz",
        snvs_index = cfg.processedDataDir / "mae" / "snvs" / "{vcf}--{rna}.vcf.gz.tbi"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} {input.vcf_file} \
        {wildcards.vcf} {input.bam_file} {output.snvs_filename} \
        {config[tools][bcftoolsCmd]} {config[tools][samtoolsCmd]}
        """

rule mae_allelicCounts:
    input:
        ncbi2ucsc = getChrMap(WORKDIR, "ncbi2ucsc"),
        ucsc2ncbi = getChrMap(WORKDIR, "ucsc2ncbi"),
        vcf_file  = lambda w: getVcf(w.rna, w.vcf),
        bam_file  = lambda w: sa.getFilePath(w.rna, 'RNA_BAM_FILE'),
        fasta     = config['mae']['genome'],
        dict      = fasta_dict(config['mae']['genome']),
        script    = WORKDIR / "MAE" / "ASEReadCounter.sh"
    output:    
        counted = cfg.processedDataDir / "mae" / "allelic_counts" / "{vcf}--{rna}.csv.gz"
    shell:
        """
        {input.script} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file} {input.bam_file} {wildcards.vcf}--{wildcards.rna} \
        {input.fasta} {config[mae][gatkIgnoreHeaderCheck]} {output.counted} \
        {config[tools][bcftoolsCmd]}
        """
## QC
rule mae_renameChrQC:
    input:
        ucsc2ncbi = getChrMap(WORKDIR, "ucsc2ncbi"),
        ncbi_vcf = getQC(format="UCSC")
    output:
        ncbi_vcf = getQC(format="NCBI")
    shell:
        """
        bcftools={config[tools][bcftoolsCmd]}
        echo 'converting from UCSC to NCBI format'
        $bcftools annotate --rename-chrs {input.ucsc2ncbi} {input.ncbi_vcf} \
            | bgzip > {output.ncbi_vcf}
        $bcftools index -t {output.ncbi_vcf}
        """

rule mae_allelicCountsQC:
    input:
        ncbi2ucsc = getChrMap(WORKDIR, "ncbi2ucsc"),
        ucsc2ncbi = getChrMap(WORKDIR, "ucsc2ncbi"),
        vcf_file_ucsc = getQC(format="UCSC"),
        vcf_file_ncbi = getQC(format="NCBI"),
        bam_file      = lambda w: sa.getFilePath(w.rna, 'RNA_BAM_FILE'),
        fasta         = config['mae']['genome'],
        dict          = fasta_dict(config['mae']['genome']),
        script_qc = WORKDIR / "QC" / "ASEReadCounter.sh",
        script_mae = WORKDIR / "MAE" / "ASEReadCounter.sh"
    output:    
        counted = cfg.processedDataDir / "mae" / "allelic_counts" / "qc_{rna}.csv.gz"
    shell:
        """
        {input.script_qc} {input.ncbi2ucsc} {input.ucsc2ncbi} \
        {input.vcf_file_ucsc} {input.vcf_file_ncbi} {input.bam_file} \
        {wildcards.rna} {input.fasta} {config[mae][gatkIgnoreHeaderCheck]} \
        {output.counted} {config[tools][bcftoolsCmd]} \
        {config[tools][samtoolsCmd]} {input.script_mae}
        """

####
rulegraph_filename = f'{config["htmlOutputPath"]}/MAE_rulegraph'
rule mae_rulegraph:
    output:
        svg = f"{rulegraph_filename}.svg",
        png = f"{rulegraph_filename}.png"
    shell:
        """
        snakemake mae --rulegraph | dot -Tsvg > {output.svg}
        snakemake mae --rulegraph | dot -Tpng > {output.png}
        """
