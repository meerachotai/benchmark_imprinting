#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('stringr', installed.packages()[,1])==FALSE) { stop("Error: R package stringr could not be loaded") }
suppressPackageStartupMessages(library(stringr))

parser = ArgumentParser()

parser$add_argument("-c", type = "character", default = "", help = "concatenated counts filename (only one contrast at a time) / prefix if you need to concatenate first")
parser$add_argument("-k", type = "character", default = "", help = "gene key file name (can have all together, need appropriate A_ID, B_ID")
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
rep = opt$r

A = opt$A
B = opt$B
A_ID = opt$a
B_ID = opt$b

need_concat = opt$C

# need_concat = T
# infile = "counts_cviA_cviB_"
# key= "gene_key.txt"
# rep = 3
# A_ID = ".1A"
# B_ID = ".1B"
# A = "cviA"
# B = "cviB"

if(need_concat == TRUE) {
  cat("Concatenating files...\n")
  counts = read.table(paste0(infile, "AxB_", 1, ".txt"), sep = "\t")
  colnames(counts) = c("feature", "AxB_1")
  dat = read.table(paste0(infile, "BxA_", 1, ".txt"), sep = "\t")
  colnames(dat) = c("feature", "BxA_1")
  counts = merge(counts, dat)
  if(rep > 1) { # not for Anderson, but for others
    for (i in 2:rep) {
      dat = read.table(paste0(infile, "AxB_", i, ".txt"), sep = "\t")
      colnames(dat) = c("feature", paste0("AxB_",i))
      counts = merge(counts, dat)
      dat = read.table(paste0(infile, "BxA_", i, ".txt"), sep = "\t")
      colnames(dat) = c("feature", paste0("BxA_",i))
      counts = merge(counts, dat)
    }
  }
  rownames(counts) = counts$feature
  counts = counts[2:length(counts)]
} else {
  counts = read.table(infile, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
}

counts <- subset(counts,!row.names(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique"))
counts$feature = row.names(counts)
counts$genome <- ifelse(grepl(A_ID,row.names(counts)),A,ifelse(grepl(B_ID,row.names(counts)),B,"Other")) 

counts_syntelogs = counts
cat("Reading gene key...\n")

gene_key <- read.table(key,sep="\t",header=T,stringsAsFactors = F)
gene_key <- subset(gene_key,!duplicated(gene_key[,1])) # remove duplicates

A.col = 1
B.col = 2
A.colname = colnames(gene_key)[A.col] 
B.colname = colnames(gene_key)[B.col]

cat("Using gene key to calculated AxB_A, AxB_B, BxA_A and BxA_B counts for all replicates\n")
counts_syntelogs$syntelog <- ifelse(counts_syntelogs$feature %in% gene_key[,A.col] | 
                                      counts_syntelogs$feature %in% gene_key[,B.col], "syntelog","non.syntelog")


syntelogs_A <- merge(subset(counts_syntelogs,counts_syntelogs$syntelog == "syntelog" & counts_syntelogs$genome == A),
                     gene_key,by.x="feature",by.y=A.colname,all.x=T)

syntelogs_B <- subset(counts_syntelogs,counts_syntelogs$syntelog == "syntelog" & counts_syntelogs$genome == B)

# aligns syntelogs, along with AxB and BxA counts
syntelogs_AB <- merge(syntelogs_A,syntelogs_B,by.x = B.colname,by.y="feature",all=T)

names(syntelogs_AB) <- gsub(".x","_A",names(syntelogs_AB),fixed = T)
names(syntelogs_AB) <- gsub(".y","_B",names(syntelogs_AB),fixed = T)

# remove unnecessary columns
syntelogs_AB = syntelogs_AB[ ,-which(names(syntelogs_AB) %in% c("genome_A","genome_B", "syntelog_A", "syntelog_B", B.colname))]
names(syntelogs_AB)[names(syntelogs_AB) == 'feature'] <- A.colname

# make first col into rownames

rownames(syntelogs_AB) <- syntelogs_AB[,1]
syntelogs_AB = syntelogs_AB[,-1]

AxB_cols = seq(1,rep*4, by = 2) # 1,3,5...
AxB_A_cols = AxB_cols[1:rep]
AxB_B_cols = AxB_cols[rep+1:(rep)]

BxA_cols = seq(2, rep*4, by = 2) # 2,4,6...
BxA_A_cols = BxA_cols[1:rep]
BxA_B_cols = BxA_cols[rep+1:(rep)]

cat(paste0("Writing counts files in form ", outprefix, "_cross_rep.txt...\n"))
# want AxB_A, AxB_B | BxA_A, BxA_B
for (i in 1:rep) {
  AxB = syntelogs_AB[,c(AxB_A_cols[i],AxB_B_cols[i])]
  row.names(AxB) = str_replace_all(row.names(AxB), A_ID, "")
  write.table(AxB, paste0(outprefix,"_AxB_", i, ".txt"), sep = "\t", col.names = F, row.names = T, quote = F) # merged counts
  write.table(cbind(row.names(AxB), AxB[,1]), paste0(outprefix,"_AxB_A_", i, ".txt"), sep = "\t", col.names = F, row.names = F, quote = F) # only A
  write.table(cbind(row.names(AxB), AxB[,2]), paste0(outprefix,"_AxB_B_", i, ".txt"), sep = "\t", col.names = F, row.names = F, quote = F) # only B
  
  BxA = syntelogs_AB[,c(BxA_A_cols[i],BxA_B_cols[i])]
  row.names(BxA) = str_replace_all(row.names(AxB), A_ID, "")
  write.table(BxA, paste0(outprefix,"_BxA_", i, ".txt"), sep = "\t", col.names = F, row.names = T, quote = F) 
  write.table(cbind(row.names(BxA), BxA[,1]), paste0(outprefix,"_BxA_A_", i, ".txt"), sep = "\t", col.names = F, row.names = F, quote = F)
  write.table(cbind(row.names(BxA), BxA[,2]), paste0(outprefix,"_BxA_B_", i, ".txt"), sep = "\t", col.names = F, row.names = F, quote = F)
}
cat("Done\n")
# write.table(AxB, paste0(outprefix, "_", rep, "_", A, "x", B, "_counts_merged.txt"), sep = "\t", col.names = F, row.names = T, quote = F) 
# write.table(cbind(row.names(AxB), AxB[,1]), paste0(outprefix,"_", rep, "_", A, "x", B, "_", A, "_counts.txt"), sep = "\t", col.names = F, row.names = F, quote = F)
# write.table(cbind(row.names(AxB), AxB[,2]), paste0(outprefix,"_", rep, "_", A, "x", B, "_", B, "_counts.txt"), sep = "\t", col.names = F, row.names = F, quote = F)

# write.table(BxA, paste0(outprefix,"_", rep,  "_", B, "x", A, "_counts_merged.txt"), sep = "\t", col.names = F, row.names = T, quote = F) 
# write.table(cbind(row.names(BxA), BxA[,1]), paste0(outprefix, "_", rep, "_", B, "x", A, "_", A, "_counts.txt"), sep = "\t", col.names = F, row.names = F, quote = F)
# write.table(cbind(row.names(BxA), BxA[,2]), paste0(outprefix, "_", rep, "_", B, "x", A, "_", B, "_counts.txt"), sep = "\t", col.names = F, row.names = F, quote = F)

