#'---
#' title: DNA-RNA matching matrix
#' author: vyepez
#' wb:
#'  params:
#'    - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'    - mat_qc: '`sm parser.getProcResultsDir()+"/mae/"+config["mae"]["qcGroup"]+"/dna_rna_qc_matrix.Rds"`'
#' output: 
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+echo=F
saveRDS(snakemake, paste0(snakemake@params$tmpdir, "qc_hist.snakemake"))
# snakemake <- readRDS(".drop/tmp/MAE/qc_hist.snakemake")

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
})


#'
#' ## Plot DNA - RNA matching matrix
qc_mat <- readRDS(snakemake@input$mat_qc)
hist(qc_mat, xlab = '% of overlapping variants from DNA and RNA', main = '')
melt_mat <- melt(qc_mat)

#' Logarithmic scale of the y axis provides a better visualization
ggplot(melt_mat, aes(value)) + geom_histogram(fill = 'cadetblue4', bins = 25) + 
  theme_bw(base_size = 14) + labs(x = '% of matching DNA - RNA variants', y = 'Count') + 
  scale_y_log10()  + xlim(c(NA,1))

#' ## Identify matching samples
identityCutoff <- .85

#' Number of samples that match with another
length(qc_mat[qc_mat > identityCutoff])

#' Median of matching samples value
median(qc_mat[qc_mat > identityCutoff])

#' Median of not matching samples value
median(qc_mat[qc_mat < identityCutoff])

# TODO: subset for qc_group only!!!
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)[, .(DNA_ID, RNA_ID)]
sa[, ANNOTATED_MATCH := TRUE]
colnames(melt_mat)[1:2] <- c('DNA_ID', 'RNA_ID')
sa <- sa[RNA_ID %in% melt_mat$RNA_ID]


#' ### Samples that were annotated to match but do not 
sa2 <- left_join(sa, melt_mat, by = c('DNA_ID', 'RNA_ID')) %>% as.data.table()
DT::datatable(sa2[value < identityCutoff])

#' ### Samples that were not annotated to match but actually do
sa3 <- left_join(melt_mat, sa, by = c('DNA_ID', 'RNA_ID')) %>% as.data.table()
DT::datatable(sa3[is.na(ANNOTATED_MATCH) & value > identityCutoff])
