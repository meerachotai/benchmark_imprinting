#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('Biostrings', installed.packages()[,1])==FALSE) { stop("Error: R package Biostrings could not be loaded") }
suppressPackageStartupMessages(library(Biostrings))

parser = ArgumentParser()

parser$add_argument("-A", type = "character", default = "", help = "FASTA file for ref strain A")

parser$add_argument("-m", type = "double", default = 1, help = "score for match (> 0 for rewarding)")
parser$add_argument("-s", type = "double", default = 1, help = "score for snp (converted to < 0 for penalizing)")
parser$add_argument("-i", type = "double", default = 3, help = "score for indel start (converted < 0 for penalizing)")
parser$add_argument("-e", type = "double", default = 1, help = "score (each) for indel extending (converted < 0 for penalizing)")

parser$add_argument("-p", type = "double", default = 20, help = "total % score needed")
parser$add_argument("-r", type = "double", default = 1, help = "snp:indel ratio")

parser$add_argument("-a", type = "double", default = 2, help = "alpha value input into simuG")
parser$add_argument("-C", type = "double", default = 0.5, help = "constant value input into simuG")

parser$add_argument("-o", type = "character", default = "strainB", help = "outdirectory to save temp files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

seq = opt$A

match = opt$m
snp = opt$s
indel = opt$i
extend = opt$e

score_percent = opt$p
s.i.ratio = opt$r

alpha = opt$a
C = opt$C

outdir = opt$o

# match = 1
# extend = 1
# indel = 2
# snp = 2
# score_percent = 90
# s.i.ratio = 1
# alpha = 2
# C = 0.5
# seq = "cviA_genome.fa"

# DO NOT TOUCH---------
per_chrom_stats = paste0(outdir, "/per_chrom_stats.txt")

seq = readDNAStringSet(seq)

# combining extend and indel for ONE score to get indel counts
# SimuG's power-law-fitted dist P(size) = C * (size ^ -alpha) from 0 <= size <= 50
average_indel_size = 0
for(size in 1:50) { # calculating expected value for indel size
  average_indel_size = average_indel_size + (C * (size^-alpha) * size)
}
average_indel_score = ((round(average_indel_size) - 1) * extend) + indel

chr = names(seq)
chr = sub("\\s.*","",chr) # only till the first space, to match simuG's chr name
seq_stats = as.data.frame(chr)
seq_stats$length = as.numeric(as.character(width(seq)))
colnames(seq_stats) = c("chr", "length")
seq_stats$length = as.numeric(as.character(seq_stats$length))
seq_stats$perfect = seq_stats$length * match

# get total required score
seq_stats$total_score = round((score_percent/100) * (seq_stats$perfect))

# get snp score and counts
seq_stats$snps = round((seq_stats$perfect - seq_stats$total_score) / (match + snp + (snp/s.i.ratio)), digits = 0)
# get indel score (using snp:indel ratio)
seq_stats$indel_score = (seq_stats$snps * snp) / s.i.ratio
# get indel counts
seq_stats$indels = round(seq_stats$indel_score / average_indel_score)

total_indels = sum(seq_stats$indels)
total_snps = sum(seq_stats$snps)

# to be read in using grep on bash, inform simuG command (example below)
cat("\nsnp-count:", total_snps)
cat("\nindel-count:", total_indels, "\n")

write.table(seq_stats, per_chrom_stats, row.names = F, quote = F, col.names = T, sep = "\t")

# perl simuG/simuG.pl -r cviA_genome.fa -indel_count 814 -snp_count 1217 -prefix cviB -seed 5
