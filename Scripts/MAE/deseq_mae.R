#'---
#' title: Get MAE results
#' author: vyepez
#' wb:
#'  input:
#'   - mae_counts: '`sm parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz" `'
#'   - v29_dt: '`sm parser.getProcDataDir() + "/mae/v29/gene_name_mapping.Rds" `'
#'  output:
#'   - mae_res: '`sm parser.getProcResultsDir() + "/mae/samples/{vcf}--{rna}_res.Rds"`'
#'  type: script
#'---

print("Started with deseq")
saveRDS(snakemake, paste0(snakemake@config$tmpdir,'/MAE/deseq_mae.snakemake'))
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, '/MAE/deseq_mae.snakemake'))

suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(magrittr)
    library(tidyr)
    library(GenomicRanges)
    library(SummarizedExperiment)
    ## LOAD tMAE package
    devtools::load_all("tMAE")
})



mae_raw <- fread(snakemake@input$mae_counts, fill=TRUE)
# mae_raw[, sample := paste(snakemake@wildcards$vcf, snakemake@wildcards$rna, sep = "--")]

# Function from MAE pkg
rmae <- run_deseq_all_mae(mae_raw) ## build test for counting REF and ALT in MAE
rmae[, sample := paste(snakemake@wildcards$vcf, snakemake@wildcards$rna, sep = "--")]


setorderv(rmae, c('chr', 'pos'))

print("Done with deseq")
saveRDS(rmae, snakemake@output$mae_res)

