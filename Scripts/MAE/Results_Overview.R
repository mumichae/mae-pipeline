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

