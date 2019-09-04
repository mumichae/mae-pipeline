#'---
#' title: Create QC matrix
#' author: vyepez
#' wb:
#'  py:
#'  - |
#'   config["rna_ids_qc"] = parser.all_rna_ids
#'   config["wes_ids_qc"] = parser.getSampleIDs(experiment=config["wes_assay"])
#'  input: 
#'    - mae_res: '`sm lambda wildcards: expand(parser.getProcDataDir() + "/mae/RNA_GT/{rna}.Rds", rna=parser.getRNAByGroup({wildcards.dataset}))`'
#'    - vcf: '`sm parser.getVCFsFilePaths(assay="wes_assay")`'
#'  output:
#'    - mat_qc: '`sm parser.getProcResultsDir() + "/mae/{dataset}/dna_rna_qc_matrix.Rds"`'
#'  threads: 50
#'  type: script
#'---

saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/MAE/qc_matrix.snakemake"))
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/MAE/qc_matrix.snakemake"))

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(magrittr)
  library(BiocParallel)
  library(data.table)
})

register(MulticoreParam(snakemake@threads))

# Read the test vcf as GRanges
gr_test <- readVcf("/s/project/genetic_diagnosis/resource/qc_ucsc.vcf.gz") %>% granges()
mcols(gr_test)$GT <- "0/0"

# Read the vcf and rna files
input_vcf <- snakemake@input$vcf
wes_samples <- snakemake@config$wes_ids_qc


rna_samples <- snakemake@config$rna_ids_qc[[snakemake@wildcards$dataset]]
mae_res <- snakemake@input$mae_res

N <- length(input_vcf)
lp <- bplapply(1:N, function(i){
  
  # Read sample vcf file
  sample <- wes_samples[i]
  param <-  ScanVcfParam(fixed=NA, info=NA, geno='GT', samples=sample, trimEmpty=TRUE) 
  vcf_sample <- readVcf(input_vcf[i], param = param, row.names = FALSE)
  # Get GRanges and add Genotype
  gr_sample <- granges(vcf_sample)
  mcols(gr_sample)$GT <- VariantAnnotation::geno(vcf_sample)$GT
  
  # Find overlaps between test and sample
  ov <- findOverlaps(gr_test, gr_sample, type = 'equal')
  gr_res <- copy(gr_test)
  mcols(gr_res)[from(ov),]$GT <- mcols(gr_sample)[to(ov),]$GT
  
  # Find simmilarity between DNA sample and RNA sample
  x <- vapply(mae_res, function(m){
    gr_rna <- readRDS(m)
    ov <- findOverlaps(gr_res, gr_rna, type = 'equal')
    gt_dna <- gr_res[from(ov)]$GT
    gt_rna <- gr_rna[to(ov)]$RNA_GT
    sum(gt_dna == gt_rna) / length(gt_dna)
  }, 1.0)
  return(x)
})

# Create a matrix
mat <- do.call(rbind, lp)
row.names(mat) <- wes_samples
colnames(mat) <- rna_samples


saveRDS(mat, snakemake@output$mat_qc)