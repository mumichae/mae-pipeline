#'---
#' title: MAE test on qc variants
#' author: vyepez
#' wb:
#'  input:
#'   - qc_counts: '`sm parser.getProcDataDir() + "/mae/allelic_counts_qc/{rna}.csv.gz" `'
#'  output:
#'   - mae_res: '`sm parser.getProcDataDir() + "/mae/RNA_GT/{rna}.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, paste0(snakemake@config$tmpdir,'/MAE/deseq_qc.snakemake'))
# snakemake <- readRDS('paste0(snakemake@config$tmpdir,'MAE/deseq_qc.snakemake'))

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(GenomicRanges)
  devtools::load_all("tMAE")
})

# Read MA counts for qc
qc_counts <- fread(snakemake@input$qc_counts, fill=TRUE)
qc_counts <- qc_counts[!is.na(position)]

# Run DESeq
rmae <- run_deseq_all_mae(qc_counts, min_cov = 10)

# Convert to granges
qc_gr <- GRanges(seqnames = rmae$chr, 
                 ranges = IRanges(start = rmae$pos, end = rmae$pos), 
                 strand = '*')
mcols(qc_gr) = DataFrame(RNA_GT = rmae$RNA_GT)

saveRDS(qc_gr, snakemake@output$mae_res)
