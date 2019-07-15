#'---
#' title: Get MAE results
#' author: vyepez
#' wb:
#'  input:
#'   - mae_counts: '`sm parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz" `'
#'  output:
#'   - mae_res: '`sm parser.getProcResultsDir() + "/mae/samples/{vcf}--{rna}_res.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, 'tmp/res_mae.Rds')
# snakemake <- readRDS('tmp/res_mae.Rds')

suppressPackageStartupMessages({
    devtools::load_all("mae/")
    library(dplyr)
    library(data.table)
    library(magrittr)
    library(tidyr)
})


mae_raw <- fread(snakemake@input$mae_counts)
# Function from MAE pkg
#rmae <- run_deseq_all_mae(mae_raw) ## build test for counting REF and ALT in MAE
#print("Done with deseq")
rmae <- mae_raw
saveRDS(rmae, snakemake@output$rmae)


