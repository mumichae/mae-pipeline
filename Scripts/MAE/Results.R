#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  input:
#'   - mae_res: '`sm lambda wildcards: expand(parser.getProcResultsDir() + "/mae/samples/{id}_res.Rds", id = parser.getMaeByGroup({wildcards.dataset}))`'
#'   - gene_name_mapping: '`sm parser.getProcDataDir() + "/mae/gene_name_mapping_{annotation}.Rds"`'
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

# Clean res
rmae <- separate(rmae, 'sample', into = c('EXOME_ID', 'RNA_ID'), sep = "--", remove = FALSE)
rmae[, c('GT', 'RNA_GT') := NULL] 

# Add gene names
gene_annot_dt <- fread(snakemake@input$gene_name_mapping)

# Subtract the genomic ranges from the annotation and results and overlap them
gene_annot_ranges <- GRanges(seqnames = gene_annot_dt$seqnames, 
                             IRanges(start = gene_annot_dt$start, end = gene_annot_dt$end), 
                             strand = gene_annot_dt$strand)
rmae_ranges <- GRanges(seqnames = rmae$chr, IRanges(start = rmae$pos, end = rmae$pos), strand = '*')

fo <- findOverlaps(rmae_ranges, gene_annot_ranges)

# Add the gene names
res_annot <- cbind(rmae[from(fo), ],  gene_annot_dt[to(fo), .(gene_name, gene_type)])

# Prioritze protein coding genes
res_annot <- rbind(res_annot[gene_type == 'protein_coding'], res_annot[gene_type != 'protein_coding'])

# Write all the other genes in another column 
res_annot[, aux := paste(chr, pos, sep = "-")]
res_annot[, N := 1:.N, by = aux]
r_other <- res_annot[N > 1, .(other_names = paste(gene_name, collapse = ',')), by = aux]
res <- left_join(res_annot[N == 1], r_other, by = 'aux') %>% as.data.table()
res[, c('aux', 'N') := NULL]

# Bring gene_name column front
res <- cbind(res[, .(gene_name)], res[, -"gene_name"])

#'
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

