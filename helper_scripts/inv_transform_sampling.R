#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))

parser = ArgumentParser()
parser$add_argument("-n", type="double", default = 5000, help = "number of elements needed from sampling")
parser$add_argument("-o", type = "character", default = "out.txt", help = "outfile default = out.txt")
parser$add_argument("-s", type="double", default = 5, help = "set seed for reproducibility")
parser$add_argument("infile", nargs=1, help = "infile to sample from")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

infile = opt$infile
pick_num = opt$n
outfile = opt$o
seed = opt$s

# pick_num = 5000
# infile = "genes.gtf"
# outfile = "genes.gtf"

set.seed(seed)
sample = read.table(infile, sep = "\t", header=F)

# inverse transform sampling:
# uniform dist [0,1]
rand = runif(nrow(sample), 0, 1)
sample$rand = rand

# qnorm is the R function that calculates the inverse cdf of the normal distribution
sample$inv = qnorm(sample$rand)
sample = sample[order(sample$inv, decreasing = T), ]
sample = sample[1:pick_num,1:(ncol(sample)-2)]

write.table(sample, outfile, sep = "\t", quote = FALSE, row.names = F, col.names = F)

