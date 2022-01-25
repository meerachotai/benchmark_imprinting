#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))

parser = ArgumentParser()

parser$add_argument("-i", type = "character", default = "", help = "inprefix of *_snp_report files")
parser$add_argument("-r", type="double", default = "", help = "number of reciprocal cross pairs")
parser$add_argument("-A", type = "character", default = "", help = "strain A name")
parser$add_argument("-B", type = "character", default = "", help = "strain B name")

parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

outprefix = opt$outprefix
A = opt$A
B = opt$B
rep = opt$r
inprefix = opt$i

# outprefix = "roth"
# A = "strainA"
# B = "strainB"
# rep = 3
# paired = TRUE

cat("Reading replicate _snp_report.bed files...\n")
i = 1
# infile = paste0(dir,"/rep_", i, "_", i, "_", A, "x", B, "_snp_report.bed") 
infile = paste0(inprefix, "AxB_", i, ".bed")


data = read.table(infile, sep = "\t")
data = subset(data, select = -c(V4))
names(data) = c("genes", "start", "end", paste0('AxB_',i,'_A'), paste0('AxB_',i,'_B'))

genes = unique(data[,1])
snpNum = vector(mode="numeric", length=length(genes))
newNames = vector(mode = "character", length = nrow(data))
for (j in 1:nrow(data)) {
  cur = data[j,1]
  index = which(genes == cur)
  snpNum[index] = snpNum[index] + 1
  newNames[j] = paste0("gene",index,"_snp",snpNum[index])
}
data$Labels = newNames

rename = cbind(data$genes, data$Labels,data$start, data$end)

cat("Writing renamed file...\n")
write.table(rename, paste0(outprefix, "_rename.txt"), sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

allData = data

for (i in 2:rep) {
  # infile = paste0("rep_", i, "_", i, "_", A, "x", B, "_snp_report.bed") 
  infile = paste0(inprefix, "AxB_", i, ".bed")
  data = read.table(infile, sep = "\t")
  data = subset(data, select = -c(V4))
  names(data) = c("genes", "start", "end", paste0('AxB_',i,'_A'), paste0('AxB_',i,'_B'))
  
  allData = merge(allData, data, by = c("start", "end", "genes"), all.x = TRUE, all.y = TRUE)
}

for (i in 1:rep) {
  # infile = paste0("rep_", i, "_", i, "_", B, "x", A, "_snp_report.bed") 
  infile = paste0(inprefix, "BxA_", i, ".bed")
  data = read.table(infile, sep = "\t")
  data = subset(data, select = -c(V4))
  names(data) = c("genes", "start", "end", paste0('BxA_',i,'_A'), paste0('BxA_',i,'_B'))
  
  allData = merge(allData, data, by = c("start", "end", "genes"), all.x = TRUE, all.y = TRUE)
}

cat("Writing concatenated renamed + counts file...\n")
saveData = allData[,!(names(allData) %in% c("start", "end", "genes"))] 
rownames(saveData) = saveData$Labels
saveData = saveData[,!(names(saveData) %in% c("Labels"))] 
write.table(saveData,paste0(outprefix, "_counts.txt"), sep = "\t", quote = FALSE)

# dataAlt = data[,5:6]
# rownames(dataAlt) = data$Labels
# write.table(dataAlt, outfile, sep = "\t", col.names = FALSE, quote = FALSE)
