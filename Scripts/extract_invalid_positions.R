#'---
#' title: Extract invalid positions in vcf
#' author: salazar, figueira
#' wb:
#'  input:
#'   - vcf: '`sm lambda wildcards: parser.getProcDataDir() + "/mae/temp_processed_{vcf}.vcf.gz"`'
#'  output:
#'   - invalidIntervals: '`sm parser.getProcDataDir() + "/mae/{vcf}_invalid.intervals"`'
#'  threads: 1
#'  type: script
#'---

suppressPackageStartupMessages({
    library(vcfR)
    library(data.table)
})

filename <- snakemake@input$vcf
vcf <- read.vcfR(filename)
vcf <- data.table(getFIX(vcf))

chrom_pos <- data.table(vcf[grepl(",",ALT, fixed=TRUE), paste0(CHROM,":",POS)])
fwrite(chrom_pos, file = snakemake@output$invalidIntervals, col.names = FALSE)


#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

#library(vcfR)
#library(data.table)

# Folder creation
#folder1 = "/s/project/mitoMultiOmics/raw_data/helmholtz/"
#folder2 = "/exomicout/paired-endout/processedData"
#file1 = "temp_processed_"
#file2 = ".vcf.gz"

#sample_number <- args[1]

#folder <- paste0(folder1,sample_number,folder2)
#file <- paste0(file1,sample_number,file2)

# Read VCF
#vcf <- read.vcfR(file.path(folder, file))
#vcf <- data.table(getFIX(vcf))

#chrom_pos <- data.table(vcf[grepl(",",ALT, fixed=TRUE), paste0(CHROM,":",POS)])
#fwrite(chrom_pos, file = file.path(folder,"invalid.intervals"), col.names = FALSE)
