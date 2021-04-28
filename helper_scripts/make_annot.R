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

# seqid feature type start end 
# note: sequence numbering starting at 1
# add ##gff-version 3 on top
# seqid = vector()
# start = vector()
# end = vector()
# id = vector()

# for(i in 1:nrow(input)) {
#   # seqid[i] = input[i,1]
#   # end[i] = input[i,2]
#   id[i] = paste0("ID=exon",i) # check with CLP
# }

seqid = input[,1]
end = input[,2]
start = rep(1, nrow(input))
feature = rep("exon", nrow(input))
fillers = rep(".", nrow(input))
strand = rep("+", nrow(input))
parent = paste0(rep("Parent:", nrow(input)),input[,1])

cat("##gff-version 3\n",file=paste0(outprefix, "_annot.gff3"))
annot = cbind(seqid, fillers, feature, start, end, fillers, strand, fillers, parent) #, input[,2])

write.table(annot, paste0(outprefix, "_annot.gff3"), append=TRUE, row.names = F, quote = F, col.names = F, sep = "\t")
