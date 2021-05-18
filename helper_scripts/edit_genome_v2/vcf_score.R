#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('vcfR', installed.packages()[,1])==FALSE) { stop("Error: R package vcfR could not be loaded") }
suppressPackageStartupMessages(library(vcfR))
if (is.element('Biostrings', installed.packages()[,1])==FALSE) { stop("Error: R package Biostrings could not be loaded") }
suppressPackageStartupMessages(library(Biostrings))

parser = ArgumentParser()

parser$add_argument("-G", type = "character", default = "", help = "SimuG generated FASTA file")
parser$add_argument("-S", type = "character", default = "", help = "vcf SNP infile for strain B")
parser$add_argument("-I", type = "character", default = "", help = "vcf indel infile for strain B")

parser$add_argument("-s", type = "double", default = 1, help = "score for snp (converted to < 0 for penalizing)")
parser$add_argument("-i", type = "double", default = 3, help = "score for indel start (converted < 0 for penalizing)")
parser$add_argument("-e", type = "double", default = 1, help = "score (each) for indel extending (converted < 0 for penalizing)")

parser$add_argument("-W", type = "double", default = 1, help = "% range window for indel score")

parser$add_argument("-B", type = "character", default = "strainB", help = "FASTA outfile name for strainB")

parser$add_argument("-o", type = "character", default = "strainB", help = "outdirectory to save/get temp files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

indel_vcf = opt$I
snp_vcf = opt$S

seq = opt$G

snp = opt$s
indel = opt$i
extend = opt$e
window = opt$W

outdir = opt$o
outfile = opt$B

# seq = "cviA_genome.fa"
# indel_vcf = "cviB_simul/cviB.refseq2simseq.INDEL.vcf"
# snps_vcf = "cviB_simul/cviB.refseq2simseq.SNP.vcf"
# window = 3
# extend = 1
# indel = 2
# snp = 2
# outfile = "cviB" # HAS to remain the same every time

# DO NOT TOUCH -----------
rows_to_remove = paste0(outdir,"/chromosomes_done.txt")
per_chrom_stats = paste0(outdir, "/per_chrom_stats.txt")

to_remove = scan(rows_to_remove, character(),quiet = TRUE)
per_chrom = read.table(per_chrom_stats, header=TRUE, sep="\t", row.names = 1)
per_chrom = per_chrom[!(row.names(per_chrom) %in% to_remove),]

cat("\nreading files:", indel_vcf, snp_vcf)
snp_data = read.vcfR(snp_vcf, verbose = F)
snp_data = as.data.frame(vcfR2tidy(snp_data, info_only = TRUE)$fix)

indel_data = read.vcfR(indel_vcf, verbose = F)
indel_data = as.data.frame(vcfR2tidy(indel_data, info_only = TRUE)$fix)

# indel_data$ref = (as.numeric(indel_data$ref_end) - as.numeric(indel_data$ref_start))
# indel_data$sim = (as.numeric(indel_data$sim_end) - as.numeric(indel_data$sim_start))
# indel_data$extend = abs(indel_data$sim - indel_data$ref) - 1
# indel_data$score = (indel_data$extend * extend ) + indel

# do snps first because it's a faster calculation (count # rows)
snps_done = vector(mode = "character")
snps_index = vector(mode = "numeric")
snp_stats = data.frame(snps_done = character(), snps_needed = integer(), snps_actual = integer())

for (i in 1:nrow(per_chrom)) {
  chrom = row.names(per_chrom)[i]
  actual_score = length(which(grepl(chrom, snp_data$CHROM)))
  needed_score = per_chrom$snps[i]
  snp_window = ceiling((window/100) * needed_score)
  if(actual_score > (needed_score - snp_window) & actual_score < needed_score + snp_window) {
    snps_done = c(snps_done, chrom)
    snps_index = c(snps_index, i)
    current = data.frame(chrom, needed_score, actual_score)
    snp_stats = rbind(snp_stats, current)
  }
}

# now only get those that have good SNP AND indel score

indels_done = vector(mode = "character")
final_stats = data.frame(chrom = character(), snps_needed = integer(), snps_actual = integer(), indels_needed = integer(), indels_actual = integer())

for (i in 1:length(snps_done)) {
  chrom = snps_done[i]
  rows = subset(indel_data, grepl(chrom, indel_data$CHROM))
  
  rows$ref = as.numeric(rows$ref_end) - as.numeric(rows$ref_start)
  rows$sim = as.numeric(rows$sim_end) - as.numeric(rows$sim_start)
  rows$extend = abs(rows$sim - rows$ref) - 1
  rows$score = (rows$extend * extend ) + indel
  
  actual_score = sum(rows$score)
  needed_score = per_chrom$indel_score[snps_index[i]]
  
  indel_window = ceiling(window/100 * needed_score) # specific for each chromosome
  
  if(actual_score > (needed_score - indel_window) & actual_score < (needed_score + indel_window)) {
    indels_done = c(indels_done, chrom)
    current = data.frame(chrom, snp_stats[i,2], snp_stats[i,3], needed_score, actual_score)
    final_stats = rbind(final_stats, current)
  }
}
# colnames(final_stats) = c("chrom", "snps_needed","snps_actual", "indels_needed", "indels_actual")

# cat("indels done:", indels_done, "\n")
# the intersect is already there for indels_done, so just append to done list
if(length(indels_done) > 0) {
  write(indels_done, rows_to_remove, append = TRUE)
  write.table(final_stats, paste0(outdir, "/score_stats.txt"), append = TRUE, row.names = F, quote = F, col.names = F, sep = "\t")
}

# get sequence, and append to final fasta file for strainB
seqB = readDNAStringSet(seq)
names(seqB) = sub("\\s.*","",names(seqB)) # remove after first space
writeXStringSet(seqB[indels_done], outfile, append = TRUE, format = "fasta")

# seqA = data.frame(as.character(seqA))
# colnames(seqA) = "seq"
# row.names(seqA) = sub("\\s.*","",row.names(seqA))
# seqA = seqA[row.names(seqA) %in% indels_done]
