#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('DESeq2', installed.packages()[,1])==FALSE) { stop("Error: R package DESeq2 could not be loaded") }
suppressPackageStartupMessages(library(DESeq2))

parser = ArgumentParser()

parser$add_argument("-c", type = "character", default = "", help = "concatenated counts filename (only one contrast at a time) / prefix if you need to concatenate first")
parser$add_argument("-k", type = "character", default = "", help = "gene key file name (can have all together, need appropriate A_ID, B_ID")

parser$add_argument("-p", type = "double", default = 0.05, help = "p-value/alpha cutoff for DESeq2")
parser$add_argument("-u", type="double", default = 0.8, help = "RER upper limit, default = 0.8")
parser$add_argument("-l", type = "double", default = 0.5, help = "RER lower limit, default = 0.5")
parser$add_argument("-f", type="double", default = 1, help = "log2fc cutoff, BxA / AxB")
parser$add_argument("-r", type="double", default = "", help = "number of replicates")

parser$add_argument("-A", type = "character", default = "", help = "strain A name")
parser$add_argument("-B", type = "character", default = "", help = "strain B name")
parser$add_argument("-a", type = "character", default = "", help = "strain A ID unique")
parser$add_argument("-b", type = "character", default = "", help = "strain B ID unique")

parser$add_argument("-C", default = FALSE, action="store_true", help = "need to concatenate files?")
parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

outprefix = opt$outprefix

key = opt$k
infile = opt$c

alpha = opt$p
upper_lim = opt$u
lower_lim = opt$l
rep = opt$r
log2fc = opt$f

A = opt$A
B = opt$B
A_ID = opt$a
B_ID = opt$b

need_concat = opt$C

# outprefix = "simul"

# alpha = 0.05
# upper_lim = 0.7
# lower_lim = 0.4
# log2fc = 1
# rep = 3 # number of replicates
# 
# key = "gene_key.txt"
# 
# infile = "counts_cviA_cviB"
# need_concat = TRUE
# 
# A = "cviA"
# B = "cviB"
# A_ID = ".1A"
# B_ID = ".1B"

if(need_concat == TRUE) {
  cat("Concatenating files...\n")
  counts = read.table(paste0(infile, "AxB_", 1, ".txt"), sep = "\t")
  colnames(counts) = c("feature", "AxB_1")
  for (i in 2:rep) {
    dat = read.table(paste0(infile, "AxB_", i, ".txt"), sep = "\t")
    colnames(dat) = c("feature", paste0("AxB_",i))
    counts = merge(counts, dat)
  }
  for (i in 1:rep) {
    dat = read.table(paste0(infile, "BxA_", i, ".txt"), sep = "\t")
    colnames(dat) = c("feature", paste0("BxA_",i))
    counts = merge(counts, dat)
  }
  rownames(counts) = counts$feature
  counts = counts[2:length(counts)]
} else {
  counts = read.table(infile, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
}


# -------------------- calculate maternal preference --------------------

cat("Calculating maternal preferences...")
# counts = read.table(infile, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)

counts <- subset(counts,!row.names(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique"))

counts$genome <- ifelse(grepl(A_ID,row.names(counts)),A,ifelse(grepl(B_ID,row.names(counts)),B,"Other")) 

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

# ------------------ finding syntelogs -----------------------

cat("Reading gene key...\n")
counts_syntelogs = counts_res

gene_key <- read.table(key,sep="\t",header=T,stringsAsFactors = F)
gene_key <- subset(gene_key,!duplicated(gene_key[,1])) # remove duplicates

# required format - A | B syntelogs side-by-side
A.col = 1
B.col = 2
A.colname = colnames(gene_key)[A.col] 
B.colname = colnames(gene_key)[B.col]


counts_syntelogs$syntelog <- ifelse(counts_syntelogs$feature %in% gene_key[,A.col] | 
                                      counts_syntelogs$feature %in% gene_key[,B.col], "syntelog","non.syntelog")


syntelogs_A <- merge(subset(counts_syntelogs,counts_syntelogs$syntelog == "syntelog" & counts_syntelogs$genome == A),
                     gene_key,by.x="feature",by.y=A.colname,all.x=T)
syntelogs_B <- subset(counts_syntelogs,counts_syntelogs$syntelog == "syntelog" & counts_syntelogs$genome == B)

syntelogs_AB <- merge(syntelogs_A,syntelogs_B,by.x = B.colname,by.y="feature",all=T)

names(syntelogs_AB) <- gsub(".x",paste0(".",A),names(syntelogs_AB),fixed = T)
names(syntelogs_AB) <- gsub(".y",paste0(".",B),names(syntelogs_AB),fixed = T)

imprint_A = grep(paste0("imprint.",A), colnames(syntelogs_AB))
imprint_B = grep(paste0("imprint.",B), colnames(syntelogs_AB))

syntelogs_AB[,imprint_A][is.na(syntelogs_AB[,imprint_A])] <- "no.imprint"
syntelogs_AB[,imprint_B][is.na(syntelogs_AB[,imprint_B])] <- "no.imprint"

# ---------------- save imprinting lists --------------------

cat("Writing imprinting lists...\n")
names(syntelogs_AB)[names(syntelogs_AB) == 'feature'] <- A.colname

MEGs_A_ID = syntelogs_AB[,2][syntelogs_AB[imprint_A] == "MEG" & syntelogs_AB[imprint_B] == "MEG"]
MEGs_B_ID = syntelogs_AB[,1][syntelogs_AB[imprint_A] == "MEG" & syntelogs_AB[imprint_B] == "MEG"]

imprinted_MEGs = data.frame(MEGs_A_ID, MEGs_B_ID)
colnames(imprinted_MEGs) = c(A, B)
write.table(imprinted_MEGs, paste0(outprefix, "_anderson_MEGs.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

PEGs_A_ID = syntelogs_AB[,2][syntelogs_AB[imprint_A] == "PEG" & syntelogs_AB[imprint_B] == "PEG"]
PEGs_B_ID = syntelogs_AB[,1][syntelogs_AB[imprint_A] == "PEG" & syntelogs_AB[imprint_B] == "PEG"]

imprinted_PEGs = data.frame(PEGs_A_ID, PEGs_B_ID)
colnames(imprinted_PEGs) = c(A, B)
write.table(imprinted_PEGs, paste0(outprefix, "_anderson_PEGs.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

# either_imprinted = syntelogs_AB[syntelogs_AB[imprint_A] != "no.imprint" | syntelogs_AB[imprint_B] != "no.imprint",]
write.table(syntelogs_AB, paste0(outprefix, "_anderson_stats.txt"), quote = F, row.names = F, col.names = T, sep = "\t")

cat("Done.\n")
