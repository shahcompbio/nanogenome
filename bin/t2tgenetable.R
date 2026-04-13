#!/usr/bin/env Rscript
library(GenomicFeatures)
library(BiocT2T)
# create gene annotation table for T2T-CHM13v2.0 structural variant annotation
# install and load the T2T TxDb
BiocT2T::install_early_t2t_txdb()
library(TxDb.Hsapiens.NCBI.CHM13v2)
txdb <- TxDb.Hsapiens.NCBI.CHM13v2

# extract gene information
g <- genes(txdb)

# build annotation table matching biomart output format
t2tAnnotation <- data.frame(
    ensembl_gene_id = g$gene_id,
    hgnc_symbol     = g$gene_id,
    gene_biotype    = NA,
    description     = NA,
    chromosome_name = sub("^chr", "", as.character(seqnames(g))),
    start_position  = start(g),
    end_position    = end(g),
    strand          = ifelse(as.character(strand(g)) == "+", 1, -1),
    stringsAsFactors = FALSE
)

# write gene annotation table
write.table(t2tAnnotation, file = "t2t-genes.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
