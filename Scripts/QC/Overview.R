#'---
#' title: Full VCF-BAM Matching Analysis over All Datasets
#' author: Michaela Mueller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - html: '`sm expand(config["htmlOutputPath"] + "/QC/{dataset}.html",
#'             dataset=config["mae"]["qcGroups"])`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "overview_qc.snakemake") )
# snakemake <- readRDS(".drop/tmp/MAE/overview_qc.snakemake")

groups <- snakemake@config$mae$qcGroups
summaries <- paste('    * [', groups ,'](', 
                   gsub(snakemake@config$htmlOutputPath, ".", 
                        snakemake@input$html), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries:
#' `r summaries`
