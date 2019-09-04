#'---
#' title: Full MAE analysis over all datasets
#' author: Michaela Mueller
#' wb:
#'  py:
#'   - |
#'    annotations = list(config["GENE_ANNOTATION"].keys())
#'    datasets = config["mae_groups"]
#'  input:
#'   - html: '`sm expand(config["htmlOutputPath"] + "/mae/{dataset}--{annotation}_results.html", annotation=annotations, dataset=datasets)`'
#' output:
#'  html_document
#'---


saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/MAE/overview.snakemake") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/MAE/overview.snakemake")

groups <- names(snakemake@config$mae_groups)
gene_annotation_names <- names(snakemake@config$GENE_ANNOTATION)

titles <- paste(gene_annotation_names, groups)
summaries <- paste('[', titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$html), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries:  `r summaries`