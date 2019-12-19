#'---
#' title: Full MAE analysis over all datasets
#' author: Michaela Mueller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - html: '`sm expand(config["htmlOutputPath"] + 
#'             "/MAE/{dataset}--{annotation}_results.html",
#'              annotation=list(config["geneAnnotation"].keys()),
#'              dataset=config["mae"]["groups"])`'
#' output:
#'  html_document
#'---


saveRDS(snakemake, file.path(snakemake@params$tmpdir, "overview.snakemake") )
# snakemake <- readRDS(".drop/tmp/MAE/overview.snakemake")

groups <- names(snakemake@config$mae$groups)
gene_annotation_names <- names(snakemake@config$geneAnnotation)

titles <- paste(gene_annotation_names, groups)
summaries <- paste('[', titles ,'](', 
                   gsub(snakemake@config$htmlOutputPath, ".", 
                        snakemake@input$html), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries:  `r summaries`
