#'---
#'
#' title: Full MAE analysis over all datasets
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - res_all: '`sm expand(parser.getProcResultsDir() + "/mae/{dataset}/MAE_results_all.Rds", dataset=config["mae_groups"])`' 
#'   - res_signif: '`sm expand(parser.getProcResultsDir() + "/mae/{dataset}/MAE_results.Rds", dataset=config["mae_groups"])`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/MAE/overview.snakemake") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/MAE/overview.snakemake")

# groups <- names(snakemake@config$outrider_filtered)
# gene_annotation_names <- names(snakemake@config$GENE_ANNOTATION)
# summaries_titles <- paste(gene_annotation_names, groups)
# summaries <- paste('[', summaries_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
# summaries <- paste(summaries, sep = '\n')
# #' Summaries:  `r summaries`