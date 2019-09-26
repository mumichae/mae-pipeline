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
    library(dplyr)
    library(data.table)
    library(magrittr)
    library(tidyr)
    library(GenomicRanges)
    library(SummarizedExperiment)
    ## LOAD tMAE package
    devtools::load_all("tMAE")
})

# Read mae counts
mae_counts <- fread(snakemake@input$mae_counts, fill=TRUE)

# Function from MAE pkg
rmae <- DESeq4MAE(mae_counts) ## build test for counting REF and ALT in MAE

print("Done with deseq")

### Add AF information from gnomAD
if(gene_assembly == 'hg19'){
    library(MafDb.gnomAD.r2.1.hs37d5)
    mafdb <- MafDb.gnomAD.r2.1.hs37d5 
} else if(gene_assembly == 'hg38'){
    library(MafDb.gnomAD.r2.1.GRCh38)
    mafdb <-MafDb.gnomAD.r2.1.GRCh38 
}

# Transform into GRanges object
gr <- GRanges(seqnames = rmae$contig, ranges = IRanges(start = rmae$position, end = rmae$position), strand = '*')
# Add score of all, African, American, East Asian and Non-Finnish European
pt <- score(mafdb, gr, pop=c('AF', 'AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe'))
# Compute the MAX_AF
pt$MAX_AF = apply(pt, 1, max, na.rm=TRUE)

rmae <- cbind(rmae, pt) %>% as.data.table()

saveRDS(rmae, snakemake@output$mae_res)
