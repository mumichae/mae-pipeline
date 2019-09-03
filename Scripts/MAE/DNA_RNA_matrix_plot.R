#'---
#' title: DNA-RNA matching matrix
#' author: vyepez
#' wb:
#'  input:
#'    - mat_qc: '`sm parser.getProcResultsDir() + "/mae/dna_rna_qc_matrix.Rds"`'
#'  output:
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


mat <- readRDS(snakemake@input$mat_qc)
hist(mat, xlab = '% of overlapping variants from DNA and RNA', main = '')

identity_cutoff <- .85

match_mat <- which(mat > identity_cutoff, arr.ind = TRUE)
match_dt <- data.table(EXOME_ID = row.names(mat)[match_mat[,1]], 
                       RNA_ID_MATCHED = colnames(mat)[match_mat[,2]],
                       ID_value = mat[match_mat])

setorder(match_dt, -ID_value)

DT::datatable(match_dt, filter = 'top')
