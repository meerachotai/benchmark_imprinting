#!/usr/bin/env Rscript

if (is.element('VennDiagram', installed.packages()[,1])==FALSE) { stop("Error: R package VennDiagram could not be loaded") }
suppressPackageStartupMessages(library(VennDiagram))
if (is.element('edgeR', installed.packages()[,1])==FALSE) { stop("Error: R package edgeR could not be loaded") }
suppressPackageStartupMessages(library(edgeR))

parser = ArgumentParser()

parser$add_argument("outdir", nargs=1, help = "outdirectory")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

outdir = opt$outdir

picard_megs = scan(paste0(outdir,"/picard_MEGs.txt"), character())
picard_pegs = scan(paste0(outdir,"/picard_PEGs.txt"), character())

wyder_megs = scan(paste0(outdir,"/wyder_MEGs.txt"), character())
wyder_pegs = scan(paste0(outdir,"/wyder_PEGs.txt"), character())

anderson_megs = scan(paste0(outdir,"/anderson_MEGs.txt"), character())
anderson_pegs = scan(paste0(outdir,"/anderson_PEGs.txt"), character())

png(paste0(outdir,'/MEGs_venn.png'))
draw.triple.venn(length(wyder_megs), length(anderson_megs),length(picard_megs), 
                 length(intersect(wyder_megs,anderson_megs)), length(intersect(anderson_megs,picard_megs)),
                 length(intersect(wyder_megs,picard_megs)),length(intersect(intersect(picard_megs,wyder_megs),anderson_megs)),
                 category = c("Wyder-MEGs", "Anderson-MEGs","Picard-MEGs"), lty = rep(1, 3), col = rep("gray30",3),scaled = TRUE,
                 fill = c("lightblue3","darkorange3", "plum3"), alpha = rep(0.5, 3), 
                 cex=1.5, cat.cex = 1.5, cat.dist = c(0.05,0.030,0.025),cat.pos = c(-15, 15,150)) 
dev.off()
png(paste0(outdir,'/PEGs_venn.png'))
draw.triple.venn(length(wyder_pegs), length(anderson_pegs),length(picard_pegs), 
                 length(intersect(wyder_pegs,anderson_pegs)), length(intersect(anderson_pegs,picard_pegs)),
                 length(intersect(wyder_pegs,picard_pegs)),length(intersect(intersect(picard_pegs,wyder_pegs),anderson_pegs)),
                 category = c("Wyder-PEGs", "Anderson-PEGs","Picard-PEGs"), lty = rep(1, 3), col = rep("gray30",3),scaled = TRUE,
                 fill = c("lightblue3","darkorange3", "plum3"), alpha = rep(0.5, 3), 
                 cex=1.5, cat.cex = 1.5, cat.dist = c(0.05,0.030,0.025),cat.pos = c(-15, 15,150))
dev.off()
# ----------------------------------------------------------

# true_megs = scan(paste0(outdir,"/true_MEGs.txt"), character())
# true_pegs = scan(paste0(outdir,"/true_PEGs.txt"), character())
# all_genes = scan(paste0(outdir,"/all_genes.txt"), character())

# # for MEGs
# true_pos = intersect(true_megs, wyder_megs)
# false_pos = wyder_megs[!wyder_megs %in% true_megs]
# false_neg = true_megs[!true_megs %in% wyder_megs]
# true_neg = intersect(all_genes[!all_genes %in% true_megs], all_genes[!all_genes %in% wyder_megs])

# # for PEGs
# true_pos = intersect(true_pegs, wyder_pegs)
# false_pos = wyder_pegs[!wyder_pegs %in% true_pegs]
# false_neg = true_pegs[!true_pegs %in% wyder_pegs]
# true_neg = intersect(all_genes[!all_genes %in% true_pegs], all_genes[!all_genes %in% wyder_pegs])

