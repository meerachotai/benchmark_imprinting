#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))

parser = ArgumentParser()

parser$add_argument("input", nargs=1, help = "input file, first two columns from samtools index")
parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")
parser$add_argument("-A", default = FALSE, action="store_true", help = "Anderson format")
parser$add_argument("-P", default = FALSE, action="store_true", help = "Picard/Gehring format")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

input = opt$input
outprefix = opt$outprefix

if(opt$A) {
  cat("\nCreating annotation file based upon Anderson requirements\n")
  input = read.table(input, header=F, sep="\t", stringsAsFactors = F)
  
  seqid = input[,1]
  end = input[,2]
  start = rep(1, nrow(input))
  feature = rep("exon", nrow(input))
  fillers = rep(".", nrow(input))
  strand = rep("+", nrow(input))
  # parent = paste0(rep("Parent:", nrow(input)),input[,1])
  id = paste0(rep("ID=", nrow(input)),input[,1])
  
  cat("##gff-version 3\n",file=outprefix)
  annot = cbind(seqid, fillers, feature, start, end, fillers, strand, fillers, id) #, input[,2])
  
  write.table(annot, outprefix, append=TRUE, row.names = F, quote = F, col.names = F, sep = "\t")
  
} else if(opt$P) {
  cat("\nCreating annotation file based upon Picard-Gehring requirements\n")
  input = read.table(input, header=F, sep="\t", stringsAsFactors = F)
  
  seqid = input[,1]
  end = input[,2]
  start = rep(1, nrow(input))
  feature = rep("exon", nrow(input))
  fillers = rep(".", nrow(input))
  strand = rep("+", nrow(input))
  gene_id = paste0(rep("gene_id=", nrow(input)),input[,1], rep(";transcript_id ",nrow(input)),input[,1])
  # separator = strand = rep(";", nrow(input))
  # transcript_id = paste0(rep("transcript_id ", nrow(input)),input[,1])
  
  annot = cbind(seqid, fillers, feature, start, end, fillers, strand, fillers, gene_id) #, separator, transcript_id) #, input[,2])
  write.table(annot, outprefix, append=TRUE, row.names = F, quote = F, col.names = F, sep = "\t")
}
