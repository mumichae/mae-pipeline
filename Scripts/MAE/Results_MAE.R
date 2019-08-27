#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  input:
#'   - mae_res: '`sm expand(parser.getProcResultsDir() + "/mae/samples/{id}_res.Rds", id = parser.getMaeIDs())`'
#'  output:
#'   - res_all: '`sm parser.getProcResultsDir() + "/mae/MAE_results_all.Rds"`' 
#'   - res_signif: '`sm parser.getProcResultsDir() + "/mae/MAE_results.Rds"`' 
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
  library(ggplot2)
  library(tidyr)
})
   
#'
#' ### Read all mae files

res <- lapply(snakemake@input$mae_res, function(m){
  rt <- readRDS(m)
  return(rt)
}) %>% rbindlist()

res <- separate(res, 'sample', into = c('EXOME_ID', 'RNA_ID'), sep = "--", remove = FALSE)
res[, c('GT', 'as_gt') := NULL] 

# Bring gene_name column front
res <- cbind(res[, .(gene_name)], res[, -"gene_name"])

#' Total number of samples
uniqueN(res$sample)

#' Total number of genes
uniqueN(res$gene_name)

#' ### Subset for significant events (pad < .05 & frequency of alternative > .8)
res[, MAE := padj < .2]
res[, MAE_ALT := MAE == TRUE & alt_freq > .8]

#' Number of samples with significant MA for alternative events
uniqueN(res[MAE_ALT == TRUE, sample])

#' Number of samples with significant MA for alternativeevents
uniqueN(res[MAE_ALT == TRUE, gene_name])

#' ### Save results 
saveRDS(res, snakemake@output$res_all)
saveRDS(res[MAE_ALT == TRUE], snakemake@output$res_signif)


#+echo=F
res[, N := .N, by = sample]
res[MAE == TRUE, N_MAE := .N, by = sample]
res[MAE_ALT == TRUE, N_MAE_ALT := .N, by = sample]

rd <- unique(res[,.(sample, N, N_MAE, N_MAE_ALT)])
melt_dt <- melt(rd)
melt_dt[variable == 'N', variable := '+10 counts']
melt_dt[variable == 'N_MAE', variable := 'MAE']
melt_dt[variable == 'N_MAE_ALT', variable := 'MAE for ALT']

#' 
#' ## Cascade plot 
ggplot(melt_dt, aes(variable, value)) + geom_boxplot() +
  scale_y_log10() + theme_bw(base_size = 14) +
  labs(y = 'Heterozygous SNVs per patient', x = '')

