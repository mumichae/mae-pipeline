#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  input:
#'   - mae_res: '`sm expand(parser.getProcResultsDir() + "/mae/samplesNEW/{id}_res.Rds", id = parser.getMaeIDs())`'
#'  output:
#'   - res_signif_all: '`sm parser.getProcResultsDir() + "/mae/MAE_results.Rds"`' 
#'  type: script
#'---

# #' output: 
# #'   html_document:
# #'    code_folding: show
# #'    code_download: TRUE
# #'---

print("BUILD RESULTS FOR MAE")
#+ echo=F
saveRDS(snakemake, 'tmp/mae_res_all.Rds')
# snakemake <- readRDS('tmp/mae_res_all.Rds')
suppressPackageStartupMessages({
    library(data.table)
    library(magrittr)
    library(ggplot2)
    library(cowplot)
    library(tidyr)
    devtools::load_all("../genetic-diagnosis-tools")
})

#source("Scripts/_functions/gene_annotation/add_gene_info_cols.R")

## Read all mae files
res <- lapply(snakemake@input$mae_res, function(m){
    rt <- readRDS(m)
    return(rt)
}) %>% rbindlist()

res <- separate(res, 'sample', into = c('EXOME_ID', 'RNA_ID'), sep = "--", remove = FALSE)
res[, c('GT', 'as_gt') := NULL] ## Do we need this?

# Add gene info
# setnames(res, "hgncid", "gene_name")
res <- add_all_gene_info(res, dis_genes = F)

# Bring gene_name column front
res <- cbind(res[, .(gene_name)], res[, -"gene_name"])

saveRDS(res, snakemake@output$res_signif_all)
