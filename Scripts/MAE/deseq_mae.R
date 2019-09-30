#'---
#' title: Get MAE results
#' author: vyepez
#' wb:
#'  input:
#'   - mae_counts: '`sm parser.getProcDataDir() + "/mae/allelic_counts/{vcf}--{rna}.csv.gz" `'
#'  output:
#'   - mae_res: '`sm parser.getProcResultsDir() + "/mae/samples/{vcf}--{rna}_res.Rds"`'
#'  type: script
#'---

print("Started with deseq")
saveRDS(snakemake, paste0(snakemake@config$tmpdir,'/MAE/deseq_mae.snakemake'))
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, '/MAE/deseq_mae.snakemake'))

suppressPackageStartupMessages({
    ## LOAD tMAE package
    devtools::load_all("tMAE")
    library(stringr)
})

# Read mae counts
mae_counts <- fread(snakemake@input$mae_counts, fill=TRUE)
mae_counts[, position := as.numeric(position)]

# Sort by chr
mae_counts$contig <- str_sort(mae_counts$contig, numeric = TRUE)

print("Running DESeq...")
# Function from MAE pkg
rmae <- DESeq4MAE(mae_counts) ## build test for counting REF and ALT in MAE

### Add AF information from gnomAD
print("Adding gnomAD allele frequencies...")
rmae <- add_gnomAD_AF(rmae, gene_assembly = snakemake@config$gene_assembly,
                      max_af_cutoff = snakemake@config$MAX_AF)

saveRDS(rmae, snakemake@output$mae_res)
