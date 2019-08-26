#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  input:
#'   - mae_res: '`sm expand(parser.getProcResultsDir() + "/mae/samplesNEW/{id}_res.Rds", id = parser.getMaeIDs())`'
#'  output:
#'   - res_signif_all: '`sm parser.getProcResultsDir() + "/mae/MAE_results.Rds"`' 
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---


#+ echo=F
saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/MAE/mae_res_all.Rds") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/MAE/mae_res_all.Rds"))

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(ggplot2)
  library(cowplot)
  library(tidyr)
})
   

#' ### Read all mae files
res <- lapply(snakemake@input$mae_res, function(m){
  rt <- readRDS(m)
  return(rt)
}) %>% rbindlist()

res <- separate(res, 'sample', into = c('EXOME_ID', 'RNA_ID'), sep = "--", remove = FALSE)
res[, c('GT', 'as_gt') := NULL] 

#' ### Bring gene_name column front
res <- cbind(res[, .(gene_name)], res[, -"gene_name"])

#' ### OUTPUT 
res

#' ### Save result
saveRDS(res, snakemake@output$res_signif_all)

