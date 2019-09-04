#'---
#' title: DNA-RNA matching matrix
#' author: vyepez
#' wb:
#'  input:
#'    - mat_qc: '`sm parser.getProcResultsDir()+"/mae/"+config["qc_group"]+"/dna_rna_qc_matrix.Rds"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/MAE/qc_hist.snakemake"))
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/MAE/qc_hist.snakemake"))

suppressPackageStartupMessages({
  library(dplyr)
})

mtxs <- readRDS(snakemake@input$mat_qc)

hist(qc_mat, xlab = '% of overlapping variants from DNA and RNA', main = '')

identity_cutoff <- .85

match_mat <- which(qc_mat > identity_cutoff, arr.ind = TRUE)
match_dt <- data.table(EXOME_ID = row.names(qc_mat)[match_mat[,1]], 
                       RNA_ID_MATCHED = colnames(qc_mat)[match_mat[,2]],
                       ID_value = qc_mat[match_mat])

setorder(match_dt, -ID_value)
DT::datatable(match_dt, filter = 'top')

