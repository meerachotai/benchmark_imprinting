#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('edgeR', installed.packages()[,1])==FALSE) { stop("Error: R package edgeR could not be loaded") }
suppressPackageStartupMessages(library(edgeR))

parser = ArgumentParser()

parser$add_argument("-c", type = "character", default = "", help = "concatenated counts filename (only one contrast at a time) / prefix if you need to concatenate first")
parser$add_argument("-i", type="character", default = "", help = "inprefix for rename file from get_counts_roth.R file")
parser$add_argument("-r", type="double", default = "", help = "number of reciprocal cross pairs")
parser$add_argument("-f", type = "double", default = 0.05, help = "edgeR pvalue/FDR-cutoff, default = 0.05")
parser$add_argument("-C", default = FALSE, action="store_true", help = "need to concatenate files?")
parser$add_argument("-m", type = "double", default = 83.3, help = "maternal cutoff (default = 83.3)")
parser$add_argument("-p", type = "double", default = 33.3, help = "paternal cutoff (default = 33.3)")

parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

outprefix = opt$outprefix

fdr_cutoff = opt$f
infile = opt$c
rep = opt$r
inprefix = opt$i
need_concat = opt$C
cutoff = opt$t
PEG_bound = opt$p
MEG_bound = opt$m
need_concat = opt$C

PEG_bound = PEG_bound/100
MEG_bound = MEG_bound/100

# setwd("~/Documents/SJ_Lab/Imprinting/Winter2022/roth_imprinting")
# infile = "roth_counts.txt"
# need_concat = FALSE
# rep = 3
# logfc = 1
# fdr_cutoff = 0.1
# cutoff = 0.3
# outprefix = "roth"
# inprefix = "roth"
# PEG_bound = 0.333
# MEG_bound = 0.833

MEG_bound = MEG_bound/100
PEG_bound = PEG_bound/100

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

names = read.table(paste0(inprefix,"_rename.txt"), sep = "\t")
names(names) = c("genes", "labels", "start", "end")

# ---------------------AxB p-values ---------------------------------------

counts[is.na(counts)] = 0
counts_summed = counts[rowSums(counts) >= 10, ] # <10 removed

counts_AxB = counts_summed[1:(rep*2)]
counts_AxB$mat_mean = rowMeans(counts_AxB[,seq(1,rep*2, by = 2)])
counts_AxB$pat_mean = rowMeans(counts_AxB[,seq(2,rep*2, by = 2)])
counts_AxB$tot_mean = counts_AxB$mat_mean + counts_AxB$pat_mean
# counts_AxB$maternal_preference = counts_AxB$mat_mean / (counts_AxB$mat_mean + counts_AxB$pat_mean)

perGene_AxB = merge(names, counts_AxB, by.x = "labels", by.y = 0)

table_AxB = data.frame()
genes = unique(perGene_AxB$genes)
for (gene in genes) {
  matMean = perGene_AxB[which(gene == perGene_AxB$genes),]$mat_mean
  totMean = perGene_AxB[which(gene == perGene_AxB$genes),]$tot_mean
  # matPref = perGene_AxB[which(gene == perGene_AxB$genes),]$maternal_preference
  
  matWeights = matMean / sum(matMean)
  totWeights = totMean / sum(totMean)
  
  matWeight = weighted.mean(matMean, matWeights)
  totWeight = weighted.mean(totMean, totWeights)
  # matPrefWeight = weighted.mean(matPref, weights)
  
  row = c(gene,matWeight, totWeight)#, matPrefWeight)
  table_AxB = rbind(table_AxB, row)
  names(table_AxB) = c("gene", "mat", "tot")#, "maternal_preference")
}
columns <-c("mat", "tot")#, "maternal_preference")
table_AxB[, columns] <- lapply(columns, function(x) as.numeric(table_AxB[[x]]))

pval = vector()
for (i in 1:nrow(table_AxB)) { 
  y = tryCatch(
    expr = {
      mat = table_AxB$mat[i]
      tot = table_AxB$tot[i]
      obs = c(mat, tot - mat)
      exp = c(2/3, 1/3)
      chisq.test(obs, p = exp)$p.value
    },
    error = function(cond) {
      return(NA)
    }
  )
  pval = c(pval, y)
}
table_AxB$pval = pval
adjp = p.adjust(table_AxB$pval, method = "fdr", n = nrow(table_AxB))
table_AxB$fdr = adjp
# table_AxB[table_AxB$fdr < 0.05,]

# ----------------------repeat for BxA ------------------------------------------------

counts_BxA = counts_summed[((rep*2) + 1):(rep*4)]
counts_BxA$mat_mean = rowMeans(counts_BxA[,seq(2,rep*2, by = 2)])
counts_BxA$pat_mean = rowMeans(counts_BxA[,seq(1,rep*2, by = 2)])
counts_BxA$tot_mean = counts_BxA$mat_mean + counts_BxA$pat_mean
counts_BxA$maternal_preference = counts_BxA$mat_mean / (counts_BxA$mat_mean + counts_BxA$pat_mean)

perGene_BxA = merge(names, counts_BxA, by.x = "labels", by.y = 0)

table_BxA = data.frame()
genes = unique(perGene_BxA$genes)
for (gene in genes) {
  matMean = perGene_BxA[which(gene == perGene_BxA$genes),]$mat_mean
  totMean = perGene_BxA[which(gene == perGene_BxA$genes),]$tot_mean
  
  matWeights = matMean / sum(matMean)
  totWeights = totMean / sum(totMean)
  
  matWeight = weighted.mean(matMean, matWeights)
  totWeight = weighted.mean(totMean, totWeights)
  
  row = c(gene,matWeight, totWeight)#, matPrefWeight)
  table_BxA = rbind(table_BxA, row)
  names(table_BxA) = c("gene", "mat", "tot")
}
columns <-c("mat", "tot")
table_BxA[, columns] <- lapply(columns, function(x) as.numeric(table_BxA[[x]]))

pval = vector()
for (i in 1:nrow(table_BxA)) { 
  y = tryCatch(
    expr = {
      mat = table_BxA$mat[i]
      tot = table_BxA$tot[i]
      obs = c(mat, tot - mat)
      exp = c(2/3, 1/3)
      chisq.test(obs, p = exp)$p.value
    },
    error = function(cond) {
      return(NA)
    }
  )
  pval = c(pval, y)
}
table_BxA$pval = pval
adjp = p.adjust(table_BxA$pval, method = "fdr", n = nrow(table_BxA))
table_BxA$fdr = adjp

# -------------------------- call imprinting ------------------------------

imprinted = merge(table_AxB,table_BxA, by = "gene")
names(imprinted) <- gsub(".x",paste0("_AxB"),names(imprinted),fixed = T)
names(imprinted) <- gsub(".y",paste0("_BxA"),names(imprinted),fixed = T)

imprinted$maternal_preference_AxB = imprinted$mat_AxB / imprinted$tot_AxB
imprinted$maternal_preference_BxA = imprinted$mat_BxA / imprinted$tot_BxA

imprinted$fdr_met = ifelse(imprinted$fdr_AxB < fdr_cutoff & imprinted$fdr_BxA < fdr_cutoff, TRUE, FALSE)

imprinted$status = ifelse(imprinted$maternal_preference_AxB < PEG_bound & imprinted$maternal_preference_BxA < PEG_bound,
                          "PEG", ifelse(imprinted$maternal_preference_AxB > MEG_bound & imprinted$maternal_preference_BxA > MEG_bound,
                                        "MEG", "not_imprinted"))

megs = imprinted[imprinted$status == "MEG" & imprinted$fdr_met == T,]$gene
pegs = imprinted[imprinted$status == "PEG" & imprinted$fdr_met == T,]$gene

cat("\n-----------------------------------\n")
cat("Roth summary:\nmaternal cutoff:", MEG_bound, ", paternal cutoff:", PEG_bound, "\n")
cat("maternally-biased: ",length(megs), "\n")
cat("paternally-biased: ",length(pegs), "\n")
cat("-----------------------------------\n")

write.table(megs, paste0(outprefix, "_MEGs.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
write.table(pegs, paste0(outprefix, "_PEGs.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

write.table(imprinted, paste0(outprefix, "_stats.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# MEGs = imprinted[imprinted$gene %in% megs, ]
# PEGs = imprinted[imprinted$gene %in% pegs, ]
# write.table(MEGs, paste0(outprefix, "_MEGs_snp_report.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
# write.table(PEGs, paste0(outprefix, "_PEGs_snp_report.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

