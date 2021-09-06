#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('edgeR', installed.packages()[,1])==FALSE) { stop("Error: R package edgeR could not be loaded") }
suppressPackageStartupMessages(library(edgeR))

parser = ArgumentParser()

parser$add_argument("-c", type = "character", default = "", help = "concatenated counts filename (only one contrast at a time) / prefix if you need to concatenate first")
parser$add_argument("-r", type="double", default = "", help = "number of reciprocal cross pairs")
parser$add_argument("-p", type = "double", default = 0.05, help = "edgeR pvalue/FDR-cutoff, default = 0.05")
parser$add_argument("-l", type="double", default = 0, help = "edgeR log2fc-cutoff (based on ratio), default = 0")
parser$add_argument("-C", default = FALSE, action="store_true", help = "need to concatenate files?")

parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

outprefix = opt$outprefix

logfc = opt$l
fdr_cutoff = opt$p
infile = opt$c
rep=opt$r

need_concat = opt$C

if(need_concat == TRUE) { # concatenate in the order AxB_1_A AxB_1_B AxB_2_A ... BxA_1_A BxA_1_B BxA_2_A ...
  cat("Concatenating files...\n")
  counts = read.table(paste0(infile, "AxB_", 1, ".txt"), sep = "\t")
  colnames(counts) = c("feature", "AxB_1_A", "AxB_1_B")
  if(rep > 1) {
    for (i in 2:rep) {
      dat = read.table(paste0(infile, "AxB_", i, ".txt"), sep = "\t")
      colnames(dat) = c("feature", paste0("AxB_",i, "_A"), paste0("AxB_",i, "_B"))
      counts = merge(counts, dat)
    }
  }
  for (i in 1:rep) {
    dat = read.table(paste0(infile, "BxA_", i, ".txt"), sep = "\t")
    colnames(dat) = c("feature", paste0("BxA_",i, "_A"), paste0("BxA_",i, "_B"))
    counts = merge(counts, dat)
  }
  rownames(counts) = counts$feature
  counts = counts[2:length(counts)]
} else {
  counts = read.table(infile, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
}

counts[is.na(counts)] = 0
counts_summed = counts[rowSums(counts) >= 10, ] # <10 removed

mother = c(rep("A",rep), rep("B", rep))
type = c(rep(c("mother","father"),rep), rep(c("father","mother"), rep)) # for BxA, the order of parents switches
cross = c(rep("1", rep), rep("2", rep))

# design <- data.frame(row.names=colnames(counts_summed), mother=ifelse(is_AxB_Sample, A, B), type=ifelse(grepl("mat_", colnames(counts_summed)),"mother","father"), cross=ifelse(is_AxB_Sample, "1","2"))
design <- data.frame(row.names=colnames(counts_summed), mother, type, cross)

edgeR.design <- model.matrix(~design$cross + design$type)
#print(colnames(edgeR.design))

edgeR <- DGEList(counts=counts_summed, genes=row.names(counts_summed)) 
#print(head(counts_summed, 5))

edgeR <- calcNormFactors(edgeR)
#print(edgeR$samples$norm.factors)

edgeR <- estimateGLMCommonDisp(edgeR, edgeR.design)
edgeR <- estimateGLMTrendedDisp(edgeR, edgeR.design)
edgeR <- estimateGLMTagwiseDisp(edgeR, edgeR.design)

edgeR.fit <- glmFit(edgeR, edgeR.design)
edgeR.tr <- glmTreat(edgeR.fit, coef="design$typemother", lfc=logfc)
biased = decideTestsDGE(edgeR.tr, p=fdr_cutoff, adjust="BH")

cat("\n-----------------------------------\n")
cat("Wyder edgeR summary:\nlogFC:", logfc, ", FDR-cutoff",fdr_cutoff, "\n")
cat("maternally-biased: ",length(biased[biased==1]), "\n") 
cat("paternally-biased: ",length(biased[biased==-1]), "\n")
cat("insignificant: ",length(biased[biased==0]), "\n")
cat("-----------------------------------\n")

# -(logfc) --> PEGs
# +(logfc) --> MEGs
cat("\nWriting Wyder-edgeR MEG and PEG lists...")
write.table(row.names(counts_summed)[decideTestsDGE(edgeR.tr, p=fdr_cutoff, adjust="BH")==-1], file=paste0(outprefix,"_PEGs.txt"), row.names=F, quote=F, col.names=F)
write.table(row.names(counts_summed)[decideTestsDGE(edgeR.tr, p=fdr_cutoff, adjust="BH")==1], file=paste0(outprefix,"_MEGs.txt"), row.names=F, quote=F, col.names=F)
cat(" Done.\n")

all_stats = topTags(edgeR.tr, n=nrow(edgeR.tr))

cat("Writing Wyder-edgeR all-stats file...")
write.table(all_stats, paste0(outprefix,"_stats.txt"), row.names=T, quote=F, col.names=NA)
cat(" Done.\n")

megs_counts = counts_summed[decideTestsDGE(edgeR.tr, p=fdr_cutoff, adjust="BH")==1, ]
pegs_counts = counts_summed[decideTestsDGE(edgeR.tr, p=fdr_cutoff, adjust="BH")==-1, ]

cat("Writing MEG/PEG + counts file...")
write.table(megs_counts,paste0(outprefix,"_MEGs_counts.txt"), row.names=T, quote=F, col.names=NA, sep = "\t")
write.table(pegs_counts,paste0(outprefix,"_PEGs_counts.txt"), row.names=T, quote=F, col.names=NA, sep = "\t")
cat(" Done.\n\n")
