#!/usr/bin/env Rscript
library(circlize)
library(dplyr)
# code for plotting circos from SV and CNA data
args <- commandArgs(trailingOnly = TRUE)
sv.data.path <- args[1]
hp1.bed.path <- args[2]
hp2.bed.path <- args[3]
circos.path <- args[4]
# functions
#' sv link colors
#'
#' define colors for SV chords in circos plot
#' @param sv_type type of SV ("INS", "DEL", "TRA", "INV", "DUP", "BND")
#' @param pal color palette
#' @param alpha opacity of link
#' @return returns link colors
#' @export
get_link_color <- function(sv_type, pal, alpha=0.5) {
  if (sv_type == "INS") {
    link.color <- pal[1]
  } else if (sv_type == "DEL") {
    link.color <- pal[2]
  } else if (sv_type == "TRA") {
    link.color <- pal[3]
  } else if (sv_type == "INV") {
    link.color <- pal[4]
  } else if (sv_type == "DUP") {
    link.color <- pal[5]
  } else if (sv_type == "BND") {
    link.color <- pal[6]
  } else {
    link.color <- NA  # Handle unknown sv_type
  }
  link.color <- adjustcolor(link.color, alpha.f = alpha)
  return(link.color)
}
# load sv data
sv.data <- read.csv(sv.data.path, sep="\t")
# load in gene info
sv.gene.data <- sv.data[, c("chrom1", "base1", "base1", "gene_name_1")]
new_cols <- c("Chromosome", "chromStart", "chromEnd", "Gene")
colnames(sv.gene.data) <- new_cols
sv.gene.data.1 <- sv.data[, c("chrom2", "base2", "base2", "gene_name_2")]
colnames(sv.gene.data.1) <- new_cols
sv.gene.data <- dplyr::bind_rows(sv.gene.data, sv.gene.data.1)
sv.gene.data <- sv.gene.data %>%
  filter(!(sv.gene.data$Gene == "")) %>%
  distinct(Gene, .keep_all = TRUE)
# load in haplotype-resolved CNA data from wakhan
hp1.bed <- read.table(hp1.bed.path, sep="\t", skip=8, header=TRUE,
                      comment.char="", stringsAsFactors=FALSE)
hp2.bed <- read.table(hp2.bed.path, sep="\t", skip=8, header=TRUE,
                      comment.char="", stringsAsFactors=FALSE)
# restrict to columns of interest
cols <- c("X.chr", "start", "end", "copynumber_state")
hp1.bed <- hp1.bed[, cols]
hp2.bed <- hp2.bed[, cols]
# merge the two dataframes
cna.data <- inner_join(hp1.bed, hp2.bed, by=c("X.chr", "start", "end"), suffix=c("_hp1", "_hp2"))
# Using pmin to clip the maximum value at 4
cna.data$copynumber_state_hp1 <- pmin(cna.data$copynumber_state_hp1, 4)
cna.data$copynumber_state_hp2 <- pmin(cna.data$copynumber_state_hp2, 4)
# take sv data and prep for circlize
# prep bed files for circlize
bed1 <- sv.data[, c("chrom1", "base1", "base1", "SV_Type")]
bed2 <- sv.data[, c("chrom2", "base2", "base2", "SV_Type")]
# give colors based on SV type ... let's do a color-blind friendly palette
pal <- palette.colors(palette = "Okabe-Ito")
link.colors <- c()
for (i in 1:nrow(bed1)){
  sv_type <- bed1$SV_Type[i]
  link.color <- get_link_color(sv_type, pal)
  link.colors <- c(link.colors, link.color)
}
# plot
out.file <- circos.path;
svg(file=out.file, height=8, width=8);

circos.initializeWithIdeogram(plotType = NULL)
# load in gene labels
circos.genomicLabels(sv.gene.data, labels.column = 4, side = "outside", cex=0.5)
circos.genomicIdeogram(species="hg38")
# add in names of chromosomes
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  terms = strsplit(chr, "chr")[[1]]
  chr = terms[2]
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rgb(runif(1), runif(1), runif(1)))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.5, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.1, bg.border = NA)
# load in CNV data
circos.genomicTrackPlotRegion(cna.data,
  ylim = c(0, 4.5),
  numeric.column = c("copynumber_state_hp1",
                     "copynumber_state_hp2"),
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region,
                        value,
                        lwd = 2,
                        type = "segment",
                        col=c("#E41A1C", "#377EB8"),
                        cex=0.1, ...)
  }
)
# use colors as before
circos.genomicLink(bed1, bed2, col = link.colors,
                   border = NA)
# Add a legend for the colors
legend("topright",
       legend = c("Insertion",
                  "Deletion",
                  "Translocation",
                  "Inversion",
                  "Duplication",
                  "Breakends",
                  "Haplotype 1 CN",
                  "Haplotype 2 CN"),
       col = c(pal[1:6], "#E41A1C", "#377EB8"),
       lty = c(1, 1, 1, 1, 1, 1, NA, NA),
       pch = c(NA, NA, NA, NA, NA, NA, 16, 16),
       cex = 0.9,
       lwd = 3.5,
       pt.cex = 1.5,
       title.adj = 0.2,
       y.intersp = 1,
       x.intersp = 0.7)

# Close the SVG graphics device
dev.off()
