#'---
#' title: MAE Results table
#' author: vyepez
#' wb:
#'  input:
#'   - mae_res: '`sm expand(parser.getProcResultsDir() + "/mae/samples/{id}_res.Rds", id = parser.getMaeIDs())`'
#'  output:
#'   - res_signif_all: '`sm parser.getProcResultsDir() + "/mae/MAE_results.Rds"`'
#'   - res_signif_rare: '`sm parser.getProcResultsDir() + "/mae/MAE_results_rare.Rds"`'
#' output: 
#'   html_document:
#'    code_folding: show
#'    code_download: TRUE
#'---

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

#' ## Read all mae files
res <- lapply(snakemake@input$mae_res, function(m){
    rt <- readRDS(m)
    # HARDCODED
    #rt <- rt[padj < .05 & alt_freq > .8]
    return(rt)
}) %>% rbindlist()

res <- separate(res, 'sample', into = c('EXOME_ID', 'RNA_ID'), sep = "-", remove = FALSE)
res[, c('GT', 'as_gt') := NULL]


# Remove mismatches
# HARDCODED 
#res <- res[! sample %in% c("EXT_JAP_PT008-103165R", "EXT_JAP_PT875-103207R", "EXT_JAP_PT1146-103229R")]

# Add gene info
setnames(res, "hgncid", "gene_name")
res <- add_all_gene_info(res, dis_genes = F)

# Add sample annotation
sa <- fread(snakemake@config$SAMPLE_ANNOTATION)
res[, gene_name := toupper(gene_name)]
res <- left_join(res, sa[, .(RNA_ID, FIBROBLAST_ID, EXOME_ID, PEDIGREE, KNOWN_MUTATION,
                             CANDIDATE_GENE, BATCH, COMMENT)],
                 by = c("EXOME_ID", "RNA_ID")) %>% as.data.table
setorder(res, padj)

#+ echo=F
res[, aux := paste(chr, pos, REF, ALT, sep = "-")]


#'
# Bring gene_name column front
res = cbind(res[, .(gene_name)], res[, -"gene_name"])
res_rare <- res[MAX_AF < .001 | is.na(MAX_AF)]
# tail(sort(table(res_rare[! aux %in%c('chr6-57513221-C-T','chr6-57513248-C-G','chr6-57513317-T-C','chrY-14107314-T-C','chrY-21154426-G-A','chr10-47133833-A-G','chr2-91809227-A-G','chr9-68433537-T-G','chr6-57512775-T-G'),aux])), 10)

saveRDS(res, snakemake@output$res_signif_all)
saveRDS(res_rare, snakemake@output$res_signif_rare)

#' ### Download results tables
write.table(res_rare, paste0(snakemake@config["webDir"], "/results/MAE_results_rare.tsv"), sep = "\t", quote = F, row.names = F)

#' [Download MAE rare results table](https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/MAE_results_rare.tsv)
DT::datatable(res_rare, caption = "MAE results", style = 'bootstrap', filter = 'top')


#' ## Plots
hist(res$alt_freq, breaks = 20)
hist(res_rare$alt_freq, breaks = 20)

ggplot(res[, .N, by = c('sample', 'BATCH')], aes(BATCH, N)) +
    geom_boxplot()
median(res[, .N, by = sample]$N)

ggplot(res_rare[, .N, by = c('sample', 'BATCH')], aes(BATCH, N)) +
    geom_boxplot()
median(res_rare[, .N, by = sample]$N)

setkey(res, sample, chr, pos, REF, ALT)
setkey(res_rare, sample, chr, pos, REF, ALT)
res[, is_rare := F]
res[res_rare, is_rare := T]

events <- res[, .N, by = c('sample', 'is_rare')]
stat_box_data <- function(y, upper_limit = max(events$N)) {
    data.frame(y = upper_limit, label = paste('median =', round(median(y), 1)))
}
ggplot(events, aes(is_rare, N)) +
    geom_boxplot() +
    stat_summary(fun.data = stat_box_data, geom = "text", vjust = -1) +
    coord_trans(y = 'log10')
