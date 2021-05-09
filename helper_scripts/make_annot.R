#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))

parser = ArgumentParser()

parser$add_argument("input", nargs=1, help = "input file, first two columns from samtools index")
parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

input = opt$input
outprefix = opt$outprefix

input = read.table(input, header=F, sep="\t", stringsAsFactors = F)

seqid = input[,1]
end = input[,2]
start = rep(1, nrow(input))
feature = rep("exon", nrow(input))
fillers = rep(".", nrow(input))
strand = rep("+", nrow(input))
# parent = paste0(rep("Parent:", nrow(input)),input[,1])
parent = paste0(rep("ID=", nrow(input)),input[,1])

cat("##gff-version 3\n",file=outprefix)
annot = cbind(seqid, fillers, feature, start, end, fillers, strand, fillers, parent) #, input[,2])

write.table(annot, outprefix, append=TRUE, row.names = F, quote = F, col.names = F, sep = "\t")
