#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('edgeR', installed.packages()[,1])==FALSE) { stop("Error: R package edgeR could not be loaded") }
suppressPackageStartupMessages(library(edgeR))

parser = ArgumentParser()

parser$add_argument("-c", type = "character", default = "", help = "concatenated counts filename (only one contrast at a time) / prefix if you need to concatenate first")
parser$add_argument("-i", type="character", default = "", help = "inprefix for rename file from get_counts_roth.R file")
parser$add_argument("-r", type="double", default = "", help = "number of reciprocal cross pairs")
parser$add_argument("-p", type = "double", default = 0.05, help = "edgeR pvalue/FDR-cutoff, default = 0.05")
parser$add_argument("-l", type="double", default = 0, help = "edgeR log2fc-cutoff (based on ratio), default = 0")
parser$add_argument("-C", default = FALSE, action="store_true", help = "need to concatenate files?")
parser$add_argument("-t", type = "double", default = 0.5, help = "threshold/cutoff for consensus of SNPs, default = 0.5")
parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

outprefix = opt$outprefix

logfc = opt$l
fdr_cutoff = opt$p
infile = opt$c
rep = opt$r
inprefix = opt$i
need_concat = opt$C
cutoff = opt$t

need_concat = opt$C

# infile = "roth_counts.txt"
# need_concat = FALSE
# rep = 3
# logfc = 1
# fdr_cutoff = 0.05
# labels = "rename.txt"
# cutoff = 0.3
# outprefix = "roth"

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


cat("Preparing design...\n")
mother = c(rep("A",rep*2), rep("B", rep*2))
type = c(rep(c("mother","father"),rep), rep(c("father","mother"), rep)) # for BxA, the order of parents switches
cross = c(rep("1", rep*2), rep("2", rep*2))

# design <- data.frame(row.names=colnames(counts_summed), mother=ifelse(is_AxB_Sample, A, B), type=ifelse(grepl("mat_", colnames(counts_summed)),"mother","father"), cross=ifelse(is_AxB_Sample, "1","2"))
design <- data.frame(row.names=colnames(counts_summed), mother, type, cross)

edgeR.design <- model.matrix(~design$cross + design$type)
#print(colnames(edgeR.design))

edgeR <- DGEList(counts=counts_summed, genes=row.names(counts_summed)) 
#print(head(counts_summed, 5))

edgeR <- calcNormFactors(edgeR)
#print(edgeR$samples$norm.factors)

cat("Running edgeR...\n")

edgeR <- estimateGLMCommonDisp(edgeR, edgeR.design)
edgeR <- estimateGLMTrendedDisp(edgeR, edgeR.design)
edgeR <- estimateGLMTagwiseDisp(edgeR, edgeR.design)

edgeR.fit <- glmFit(edgeR, edgeR.design)
edgeR.tr <- glmTreat(edgeR.fit, coef="design$typemother", lfc=logfc)
biased = as.data.frame(decideTestsDGE(edgeR.tr, p=fdr_cutoff, adjust="BH"))

# cat("\n-----------------------------------\n")
# cat("Wyder edgeR summary:\nlogFC:", logfc, ", FDR-cutoff",fdr_cutoff, "\n")
# cat("maternally-biased: ",length(row.names(biased)[biased==1]), "\n") 
# cat("paternally-biased: ",length(row.names(biased)[biased==-1]), "\n")
# cat("insignificant: ",length(row.names(biased)[biased==0]), "\n")
# cat("-----------------------------------\n")

pegs = row.names(counts_summed)[decideTestsDGE(edgeR.tr, p=fdr_cutoff, adjust="BH")==-1]
megs = row.names(counts_summed)[decideTestsDGE(edgeR.tr, p=fdr_cutoff, adjust="BH")==1]

imprinted = ifelse(row.names(counts) %in% megs, "maternal", 
                    ifelse(row.names(counts) %in% pegs,"paternal", "none"))

status = data.frame(cbind(row.names(counts), imprinted))
names(status) = c("labels", "imprint")
write.table(status, paste0(outprefix, "_imprint.txt"), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

# --------------------- PART 2 ----------------------------------------

names = read.table(paste0(inprefix,"_rename.txt"), sep = "\t")
names(names) = c("genes", "labels", "start", "end")

imprint = merge(names, status, by = "labels")

genes = unique(imprint$genes)
matCount = vector("numeric", length = length(genes))
patCount = vector("numeric", length = length(genes))
for (i in 1:length(genes)) {
  selectGene = imprint[imprint$genes == genes[i],] 
  selectMat = selectGene[selectGene$imprint == "maternal",]
  matCount[i] = (nrow(selectMat) / nrow(selectGene))
  selectPat = selectGene[selectGene$imprint == "paternal",]
  patCount[i] = (nrow(selectPat) / nrow(selectGene))
}

megs = genes[matCount >= cutoff]
pegs = genes[patCount >= cutoff]

cat("\n-----------------------------------\n")
cat("Roth edgeR summary:\nlogFC:", logfc, ", FDR-cutoff",fdr_cutoff, ", cutoff", cutoff,"\n")
cat("maternally-biased: ",length(megs), "\n")
cat("paternally-biased: ",length(pegs), "\n")
cat("-----------------------------------\n")

write.table(megs, paste0(outprefix, "_MEGs.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
write.table(pegs, paste0(outprefix, "_PEGs.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

MEGs = imprint[imprint$genes %in% megs, ]
MEGs = subset(MEGs, select = -c(labels, imprint))
PEGs = imprint[imprint$genes %in% pegs, ]
PEGs = subset(PEGs, select = -c(labels, imprint))

write.table(MEGs, paste0(outprefix, "_MEGs_snp_report.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
write.table(PEGs, paste0(outprefix, "_PEGs_snp_report.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

