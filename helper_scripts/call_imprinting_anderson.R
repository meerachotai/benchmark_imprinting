#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('DESeq2', installed.packages()[,1])==FALSE) { stop("Error: R package DESeq2 could not be loaded") }
suppressPackageStartupMessages(library(DESeq2))

parser = ArgumentParser()

parser$add_argument("-c", type = "character", default = "", help = "concatenated counts filename (only one contrast at a time) / prefix if you need to concatenate first")

parser$add_argument("-p", type = "double", default = 0.05, help = "p-value/alpha cutoff for DESeq2")
parser$add_argument("-u", type="double", default = 0.8, help = "RER upper limit, default = 0.8")
parser$add_argument("-l", type = "double", default = 0.5, help = "RER lower limit, default = 0.5")
parser$add_argument("-f", type="double", default = 1, help = "log2fc cutoff, BxA / AxB")
parser$add_argument("-r", type="double", default = "", help = "number of replicates")

parser$add_argument("-A", type = "character", default = "", help = "strain A name")
parser$add_argument("-B", type = "character", default = "", help = "strain B name")

parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

outprefix = opt$outprefix

infile = opt$c

alpha = opt$p
upper_lim = opt$u
lower_lim = opt$l
rep = opt$r
log2fc = opt$f

A = opt$A
B = opt$B

# infile = "counts_"
# A = "cviA"
# B = "cviB"
# 
# alpha = 0.05
# upper_lim = 0.7
# lower_lim = 0.4
# log2fc = 1
# rep = 3 # number of replicates


cat("Concatenating files with prefix \"", infile, "\"...\n")

counts_AxB = read.table(paste0(infile, "AxB_", 1, ".txt"), sep = "\t", row.names = 1)
counts_A = data.frame(row.names(counts_AxB), counts_AxB[,1])
names(counts_A) = c("feature", "AxB_1")
counts_B = data.frame(row.names(counts_AxB), counts_AxB[,2])
names(counts_B) = c("feature", "AxB_1")

for (i in 2:rep) {
  dat = read.table(paste0(infile, "AxB_", i, ".txt"), sep = "\t", row.names = 1)
  dat_A = data.frame(row.names(dat), dat[,1])
  names(dat_A) = c("feature", paste0("AxB_", i))
  counts_A = merge(counts_A, dat_A)
  
  dat_B = data.frame(row.names(dat), dat[,2])
  names(dat_B) = c("feature", paste0("AxB_", i))
  counts_B = merge(counts_B, dat_B)
}

for (i in 1:rep) {
  dat = read.table(paste0(infile, "BxA_", i, ".txt"), sep = "\t", row.names = 1)
  dat_A = data.frame(row.names(dat), dat[,1])
  names(dat_A) = c("feature", paste0("BxA_", i))
  counts_A = merge(counts_A, dat_A)
  
  dat_B = data.frame(row.names(dat), dat[,2])
  names(dat_B) = c("feature", paste0("BxA_", i))
  counts_B = merge(counts_B, dat_B)
}

counts_A$genome = A

counts_B$feature = paste0(counts_B$feature, "_B")
counts_B$genome = B

# CREATE GENE KEY HERE
gene_key = data.frame(counts_A$feature, counts_B$feature)
names(gene_key) = c("A", "B")

counts = data.frame(rbind(counts_A, counts_B)) # stack them

rownames(counts) = counts$feature
counts = counts[2:length(counts)]

counts$feature = row.names(counts)

counts$A_mean <- rowMeans(counts[,1:rep])
counts$B_mean <- rowMeans(counts[,(rep+1):((2*rep))])

counts$maternal_preference <- ifelse(counts$genome == A,counts$A_mean/(counts$A_mean + counts$B_mean),
                                     counts$B_mean/(counts$A_mean + counts$B_mean))

# -------------------------- DESeq2 -----------------------------

cat("Using DESeq2 for calculating maternal/paternal bias...\n")
counts_DE = counts[,1:(rep*2)] # extract only counts
counts_DE <- subset(counts_DE,rowSums(counts_DE) >= 10) 

sample_list = c(rep("AxB",rep),rep("BxA",rep))
sample_info = colnames(counts_DE)
sample_info = as.data.frame(cbind(sample_info,sample_list))
row.names(sample_info) = colnames(counts_DE)

# head(sample_info)

dds <- DESeqDataSetFromMatrix(countData = counts_DE,colData = sample_info, design = ~ sample_list)
dds <- DESeq(dds)

# Wald tests of significance: |log2fc| > lfcthreshold - two-tailed
# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/05b_wald_test_results.html
res <- results(dds, lfcThreshold=log2fc, altHypothesis = "greaterAbs")

res = as.data.frame(res)

# ---------------------- classifying MEGs, PEGs --------------

counts_res <- merge(counts,res,by.x="feature",by.y="row.names",all=F)
counts_res <- subset(counts_res,!is.na(counts_res$padj))

# up --> log2fc < -1, matpreference > 0.9 or matpreference < 0.1
# down --> log2fc > 1, matpreference > 0.9 or matpreference < 0.1
counts_res$category <- ifelse(counts_res$padj < alpha & counts_res$log2FoldChange < -(log2fc) & 
                                (counts_res$maternal_preference > upper_lim | counts_res$maternal_preference < lower_lim), "up",
                              ifelse(counts_res$padj < alpha & counts_res$log2FoldChange > (log2fc) & 
                                       (counts_res$maternal_preference > upper_lim | counts_res$maternal_preference < lower_lim),"down","notDE"))

# MEG --> up, genome = A or down, genome = B
# PEG --> down, genome = A or up, genome = B
counts_res$imprint <- ifelse(counts_res$category =="up" & counts_res$genome == A,"MEG",
                             ifelse(counts_res$category =="down" & counts_res$genome == A,"PEG",
                                    ifelse(counts_res$category =="down" & counts_res$genome == B,"MEG",
                                           ifelse(counts_res$category =="up" & counts_res$genome == B,"PEG","no.imprint"))))

# PEGs are not being caught because of strict lims, repeated with new limits from original - 0.8, 0.5
# --------------------------------------------------------------------
# required format - A | B syntelogs side-by-side
A.col = 1
B.col = 2
A.colname = colnames(gene_key)[A.col] 
B.colname = colnames(gene_key)[B.col]

counts_syntelogs = counts_res

counts_syntelogs$syntelog <- ifelse(counts_syntelogs$feature %in% gene_key[,A.col] | 
                                      counts_syntelogs$feature %in% gene_key[,B.col], "syntelog","non.syntelog")


syntelogs_A <- merge(subset(counts_syntelogs,counts_syntelogs$syntelog == "syntelog" & counts_syntelogs$genome == A),
                     gene_key,by.x="feature",by.y=A.colname,all.x=T)
syntelogs_B <- subset(counts_syntelogs,counts_syntelogs$syntelog == "syntelog" & counts_syntelogs$genome == B)

syntelogs_AB <- merge(syntelogs_A,syntelogs_B,by.x = B.colname,by.y="feature",all=T)

names(syntelogs_AB) <- gsub(".x",paste0(".",A),names(syntelogs_AB),fixed = T)
names(syntelogs_AB) <- gsub(".y",paste0(".",B),names(syntelogs_AB),fixed = T)

imprint_A = grep(paste0("imprint.",A), colnames(syntelogs_AB)) # column number
imprint_B = grep(paste0("imprint.",B), colnames(syntelogs_AB))

syntelogs_AB[,imprint_A][is.na(syntelogs_AB[,imprint_A])] <- "no.imprint" # NAs are not imprinted
syntelogs_AB[,imprint_B][is.na(syntelogs_AB[,imprint_B])] <- "no.imprint"

# ---------------- save imprinting lists --------------------

names(syntelogs_AB)[names(syntelogs_AB) == 'feature'] <- A.colname

imprinted_MEGs = syntelogs_AB[,2][syntelogs_AB[imprint_A] == "MEG" & syntelogs_AB[imprint_B] == "MEG"] # both directions imprinted

imprinted_PEGs = syntelogs_AB[,2][syntelogs_AB[imprint_A] == "PEG" & syntelogs_AB[imprint_B] == "PEG"]

cat("\n-----------------------------------\n")
cat("Anderson/DESeq2 imprinted genes summary:\n")
cat("logFC cutoff: ", log2fc, ", p-value cutoff: ",alpha, "\n")
cat("maternally-biased: ",length(imprinted_MEGs), "\n") 
cat("paternally-biased: ",length(imprinted_PEGs), "\n")
cat("-----------------------------------\n")

cat("\nWriting imprinting lists...\n")
write.table(imprinted_MEGs, paste0(outprefix, "_MEGs.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
write.table(imprinted_PEGs, paste0(outprefix, "_PEGs.txt"), quote = F, row.names = F, col.names = F, sep = "\t")

syntelogs_AB = syntelogs_AB[-1]
# either_imprinted = syntelogs_AB[syntelogs_AB[imprint_A] != "no.imprint" | syntelogs_AB[imprint_B] != "no.imprint",]
write.table(syntelogs_AB, paste0(outprefix, "_stats.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


cat("Done.\n")
