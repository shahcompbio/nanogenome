#!/usr/bin/env Rscript
library(karyoploteR)
library(dplyr)
# code for plotting circos from SV and CNA data
args <- commandArgs(trailingOnly = TRUE)
sv.data.path <- args[1]
out.file <- args[2]
genome.version <- args[3]
dd <- read.table(sv.data.path, sep="\t", header=TRUE, stringsAsFactors = FALSE)
# filter to ensure only things on same chromosome are shown
dd.filtered <- dd %>% filter(chrom1 == chrom2)
# Convert to GRanges object
sv.data <- toGRanges(dd.filtered[,c("chrom1", "base1", "base2", "SV_Type")])
# Assign colors based on SV type
sv.colors <- ifelse(dd.filtered$SV_Type == "INS", "#CC6677",
                    ifelse(dd.filtered$SV_Type == "DEL", "#332288", "#CCBB44"))

# plot
svg(file=out.file, height=11, width=8);
# Create karyotype plot
kp <- plotKaryotype(genome=genome.version, plot.type = 2)
kpPlotDensity(kp, data=sv.data)
kpPlotRegions(kp, data=sv.data, col=sv.colors, data.panel = 2)
legend("bottomright",
       legend=c("Insertion", "Deletion", "Other"),
       fill = c("#CC6677", "#332288", "#CCBB44"),
       title.adj = 0.2,
       title = "SV Types"
)
# Close the SVG graphics device
dev.off()
