#'---
#' title: Get MAE results
#' author: vyepez
#' wb:
#'  input:
#'   - mae_counts: '`sm parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz" `'
#'   - v29_dt: '`sm parser.getProcDataDir() + "/mae/v29/gene_name_mapping.Rds" `'
#'  output:
#'   - mae_res: '`sm parser.getProcResultsDir() + "/mae/samplesNEW/{vcf}--{rna}_res.Rds"`'
#'  type: script
#'---

print("Started with deseq")
saveRDS(snakemake, paste0(snakemake@config$tmpdir,'/MAE/res_mae.snakemake'))
# snakemake <- readRDS('paste0(snakemake@config$tmpdir,'MAE/res_mae.snakemake'))

suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(magrittr)
    library(tidyr)
    library(GenomicRanges)
    library(SummarizedExperiment)
})


source("src/R/runDESeq2ForMAE.R")

mae_raw <- fread(snakemake@input$mae_counts, fill=TRUE)
mae_raw[, sample := paste(snakemake@wildcards$vcf, snakemake@wildcards$rna, sep = "--")]
# Function from MAE pkg
#print(names(mae_raw))
rmae <- run_deseq_all_mae(mae_raw) ## build test for counting REF and ALT in MAE


v29_dt <- fread(snakemake@input$v29_dt)  # gene_annotation, TODO: fix saving it as tsv!!!

# Subtract the genomic ranges from the annotation and results and overlap them
v29_granges <- GRanges(seqnames = v29_dt$seqnames, IRanges(start = v29_dt$start, end = v29_dt$end), strand = v29_dt$strand)
rmae_ranges <- GRanges(seqnames = rmae$chr, IRanges(start = rmae$pos, end = rmae$pos), strand = '*')

fo <- findOverlaps(rmae_ranges, v29_granges)

# Add the gene names
rmae_annot <- cbind(rmae[from(fo), ],  v29_dt[to(fo), .(gene_name, gene_type)])

# Prioritze protein coding genes
rmae_annot <- rbind(rmae_annot[gene_type == 'protein_coding'], rmae_annot[gene_type != 'protein_coding'])

# Write all the other genes in another column 
rmae_annot[, aux := paste(chr, pos, sep = "-")]
rmae_annot[, N := 1:.N, by = aux]
r_other <- rmae_annot[N > 1, .(other_names = paste(gene_name, collapse = ',')), by = aux]
rmae_unique <- left_join(rmae_annot[N == 1], r_other, by = 'aux') %>% as.data.table()
rmae_unique[, c('aux', 'N') := NULL]

setorderv(rmae_unique, c('chr', 'pos'))

print("Done with deseq")
saveRDS(rmae_unique, snakemake@output$mae_res)



