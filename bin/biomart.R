#!/usr/bin/env Rscript
library(biomaRt)
# create gene annotation table for structural variant annotation
args <- commandArgs(trailingOnly = TRUE)
genome <- args[1]
# fetch gene table from biomart
if (genome == 'hg38') {
    host = 'https://useast.ensembl.org'
} else if (genome == 'hg19') {
    host = 'https://grch37.ensembl.org'
} else {
    stop("Only hg19 and hg38 are supported")
}
# Disable biomaRt cache
options(timeout = 600)
Sys.setenv("BIOMART_CACHE" = "/tmp/biomart_cache")
dir.create("/tmp/biomart_cache", showWarnings = FALSE, recursive = TRUE)
# get annotations
ensembl <- useMart(host=host, biomart='ENSEMBL_MART_ENSEMBL',dataset = 'hsapiens_gene_ensembl')
ensemblAnnotation <- getBM(attributes=c('ensembl_gene_id',
                                    'hgnc_symbol',
                                    'gene_biotype',
                                    'description',
                                    'chromosome_name',
                                    'start_position','end_position',
                                    'strand'),
                                    mart = ensembl)
# create gene annotation table
write.table(ensemblAnnotation,file=paste(genome,'-genes.txt',sep=''),col.names=T,row.names=F,sep='\t',quote=F)
