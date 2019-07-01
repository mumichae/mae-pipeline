#'---
#' title: Run MAE for a sample
#' author: mumichae
#' wb:
#'  input:
#'   - vcf: '`sm lambda wildcards: parser.getFilePath(sampleId=wildcards.vcf, assay_name="dna_assay") `'
#'   - rna: '`sm lambda wildcards: parser.getFilePath(sampleId=wildcards.rna, assay_name="rna_assay")`'
#'  output:
#'   - mae: '`sm parser.getProcDataDir() + "/mae/{vcf}--{rna}.Rds"`'
#'  threads: 1
#'  type: script
#'---

# #'   - vcf: '`sm standardFileNames("Data/helmholtz/{vcf}/exomicout/paired-endout/stdFilenames/{vcf}.vcf.gz")`'
# #'   - rna: '`sm standardFileNames("Data/helmholtz/{rna}/RNAout/paired-endout/stdFilenames/{rna}.bam")`'

saveRDS(snakemake, 'tmp/mae.Rds')
snakemake <-  readRDS('tmp/mae.Rds')

suppressPackageStartupMessages({
    devtools::load_all("mae/") ########### ADD THIS TO genetic-diagnosis-tools
    library(VariantAnnotation)
    devtools::load_all("../genetic-diagnosis-tools")
})

#source("Scripts/_functions/filter_sets.R")
vcfs <- snakemake@input$vcf
rnas <- snakemake@input$rna


BPPARAM = MulticoreParam(snakemake@threads, snakemake@threads, progressbar=TRUE)
#gr <- countMAEReads(vcfs, rnas, filter_function = filter_vcf_quality, BPPARAM=BPPARAM, genome=snakemake@config$gene_assembly)  # already filters for quality
gr <- countMAEReads(vcfs, rnas, BPPARAM=BPPARAM, genome=snakemake@config$gene_assembly, seqlevelsStyle_bam=snakemake@config$CHROMOSOME_FORMAT_rna, seqlevelsStyle_vcf=snakemake@config$CHROMOSOME_FORMAT_dna ) ### filter_vcf_quality hard-coded


# TODO: add vcf and rna as input. It can be that the vcf file doesn't have the id imbedded.
mcols(gr)$EXOME_ID <- snakemake@wildcards$vcf
mcols(gr)$RNA_ID <- snakemake@wildcards$rna
mcols(gr)$sample <- paste(mcols(gr)$EXOME_ID, mcols(gr)$RNA_ID, sep = "-")

saveRDS(gr, snakemake@output$mae)
message(paste("Finished with allelic counts for sample", snakemake@wildcards$vcf))