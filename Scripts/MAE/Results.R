#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  input:
#'   - mae_res: '`sm lambda wildcards: expand(parser.getProcResultsDir() + "/mae/samples/{id}_res.Rds", id = parser.getMaeByGroup({wildcards.dataset}))`'
#'   - gene_name_mapping: '`sm parser.getProcDataDir() + "/mae/gene_name_mapping_{annotation}.tsv"`'
#'  output:
#'   - res_all: '`sm parser.getProcResultsDir() + "/mae/{dataset}/MAE_results_all_{annotation}.Rds"`' 
#'   - res_signif: '`sm parser.getProcResultsDir() + "/mae/{dataset}/MAE_results_{annotation}.Rds"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/mae/{dataset}--{annotation}_results.html"`'
#'  type: noindex
#'---

#+ echo=F
saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/MAE/mae_res_all.Rds") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/MAE/mae_res_all.Rds"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(tidyr)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(dplyr)
})
   
# Read all MAE results files
rmae <- lapply(snakemake@input$mae_res, function(m){
  rt <- readRDS(m)
  return(rt)
}) %>% rbindlist()


# Add gene names
gene_annot_dt <- fread(snakemake@input$gene_name_mapping)

# Subtract the genomic ranges from the annotation and results and overlap them
gene_annot_ranges <- GRanges(seqnames = gene_annot_dt$seqnames, 
                             IRanges(start = gene_annot_dt$start, end = gene_annot_dt$end), 
                             strand = gene_annot_dt$strand)
rmae_ranges <- GRanges(seqnames = rmae$contig, IRanges(start = rmae$position, end = rmae$position), strand = '*')

fo <- findOverlaps(rmae_ranges, gene_annot_ranges)

# Add the gene names
res_annot <- cbind(rmae[from(fo), ],  gene_annot_dt[to(fo), .(gene_name, gene_type)])

# Prioritze protein coding genes
res_annot <- rbind(res_annot[gene_type == 'protein_coding'], res_annot[gene_type != 'protein_coding'])

# Write all the other genes in another column
res_annot[, aux := paste(contig, position, sep = "-")]
rvar <- unique(res_annot[, .(aux, gene_name)])
rvar[, N := 1:.N, by = aux]

r_other <- rvar[N > 1, .(other_names = paste(gene_name, collapse = ',')), by = aux]
res <- left_join(res_annot, r_other, by = 'aux') %>% as.data.table()
res[, c('aux') := NULL]

# Bring gene_name column front
res <- cbind(res[, .(gene_name)], res[, -"gene_name"])

#'
#' Total number of samples
uniqueN(res$MAE_ID)

#' Total number of genes
uniqueN(res$gene_name)

#' ### Subset for significant events (padj < .05 & frequency of alternative > .8)
res[, MAE := padj < .05]
res[, MAE_ALT := MAE == TRUE & altFreq > .8]

#' Number of samples with significant MA for alternative events
uniqueN(res[MAE_ALT == TRUE, MAE_ID])

#' ### Save results 
saveRDS(res, snakemake@output$res_all)
saveRDS(res[MAE_ALT == TRUE], snakemake@output$res_signif)


#+echo=F
res[, N := .N, by = MAE_ID]
res[MAE == TRUE, N_MAE := .N, by = MAE_ID]
res[MAE_ALT == TRUE, N_MAE_ALT := .N, by = MAE_ID]
res[MAE_ALT == TRUE & rare == TRUE, N_MAE_ALT_RARE := .N, by = MAE_ID]

rd <- unique(res[,.(MAE_ID, N, N_MAE, N_MAE_ALT, N_MAE_ALT_RARE)])
melt_dt <- melt(rd)
melt_dt[variable == 'N', variable := '+10 counts']
melt_dt[variable == 'N_MAE', variable := 'MAE']
melt_dt[variable == 'N_MAE_ALT', variable := 'MAE for ALT']
melt_dt[variable == 'N_MAE_ALT_RARE', variable := 'MAE for ALT\n& RARE']

#' 
#' ## Cascade plot 
ggplot(melt_dt, aes(variable, value)) + geom_boxplot() +
  scale_y_log10() + theme_bw(base_size = 14) +
  labs(y = 'Heterozygous SNVs per patient', x = '')


#' 
#' ## Results table
DT::datatable(res[MAE_ALT == TRUE], filter = 'top')
