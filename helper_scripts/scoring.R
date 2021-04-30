#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('tidyverse', installed.packages()[,1])==FALSE) { stop("Error: R package tidyverse could not be loaded") }
suppressPackageStartupMessages(library(tidyverse))

parser = ArgumentParser()

parser$add_argument("-A", type = "character", default = "", help = "gene_id + sequences for ref strain A")
parser$add_argument("-c", type = "character", default = "", help = "chrom number on table")
parser$add_argument("-S", type = "character", default = "", help = "vcf SNP infile for strain B")
parser$add_argument("-I", type = "character", default = "", help = "vcf indel infile for strain B")

parser$add_argument("-m", type = "double", default = 1, help = "score for match (> 0 for rewarding)")
parser$add_argument("-s", type = "double", default = 1, help = "score for snp (converted to < 0 for penalizing)")
parser$add_argument("-i", type = "double", default = 3, help = "score for indel start (converted < 0 for penalizing)")
parser$add_argument("-e", type = "double", default = 1, help = "score (each) for indel extending (converted < 0 for penalizing)")

parser$add_argument("-p", type = "double", default = 20, help = "total % score needed")
parser$add_argument("-r", type = "double", default = 1, help = "snp:indel ratio")

parser$add_argument("-a", type = "double", default = 2, help = "alpha value input into simuG")
parser$add_argument("-C", type = "double", default = 0.5, help = "constant value input into simuG")

parser$add_argument("-P", default = FALSE, action="store_true", help = "get pre-simuG reqs")
parser$add_argument("-V", default = FALSE, action="store_true", help = "get vcf file scores ")

parser$add_argument("-W", type = "double", default = 1, help = "% range window for indel score")

# options

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

seq_A = opt$A
chrom_num = opt$c
indel_file = opt$I
snp_file = opt$S

match = opt$m
snp = opt$s
indel = opt$i
extend = opt$e

score_percent = opt$p
s.i.ratio = opt$r

alpha = opt$a
C = opt$C

window = opt$W

# chrom_num = 1
# snp_file = "strainB.refseq2simseq.SNP.vcf"
# indel_file = "strainB_trial.refseq2simseq.INDEL.vcf"
# seq_A = "colA_seq.txt"
# match = 1
# snp = 2
# indel = 2
# extend = 2.5
# s.i.ratio = 1
# score_percent = 90
# alpha = 2
# C = 0.5

seq_A = read.table(seq_A, sep = "\t", colClasses = "character", header=F)
colnames(seq_A) = c("chrom", "seq")
chrom = seq_A[chrom_num, 2]

cat("\nworking at:",score_percent, "%")
cat("\nworking on: chromosome #", chrom_num)

length = nchar(chrom)
perfect_score = length * match

# required score = ___ % of the score if everything matched
# alternately = ___ % of the length
score = round((score_percent / 100) * perfect_score)

cat("\nlength of chromosome:", length)
cat("\nrequired-score:", score)

# ------------------- find snp and indel count ---------------------

# snp_count = round((perfect_score - score) / (match + snp + (1/s.i.ratio)), digits = 0)
snp_count = round((perfect_score - score) / (match + snp + (snp/s.i.ratio)), digits = 0)
need_indel_score = (snp_count * snp) / s.i.ratio

# to check, should get score
# got_score = ((length - snp_count) * match) - (snp_count * snp) - (snp_count * snp * (1/s.i.ratio))

# power-law-fitted dist P(size) = C * (size ^ -alpha) from 0 <= size <= 50

#  using (default) alpha = 2, C = 0.5
# TODO: edit to allow other alpha and constant values as input

average_indel_size = 0
for(size in 1:50) {
  average_indel_size = average_indel_size + (C * (size^-alpha) * size)
}

average_score = ((round(average_indel_size) - 1) * extend) + indel
indel_count = round(need_indel_score / average_score)

if(opt$P) {  
  cat("\nsnp-count:", snp_count)
  cat("\nsnp-score-needed:", snp_count * snp)
  cat("\nindel-score-needed:", need_indel_score)
  cat("\nusing average indel size", average_indel_size)
  cat("\nindel-count:", indel_count, "\n")
}

# use this snp_count for simuG, use rejection sampling to get the given indel score
# perl simuG/simuG.pl -r simul_trials/strainA.fa -indel_count 300 -snp_count 50 -prefix simul_trials/strainB -s 5

# ------------------- find indel score -------------------

# scoring for rejection sampling 
if(opt$V) {
  
  # if((file.info(indel_file)$size != 0) & (file.info(snp_file)$size != 0)) {
  #   indel_file = read.table(indel_file, sep = "\t")
  #   snp_file = read.table(snp_file, sep = "\t")
  # } else {
  #   stop("Error, either snp/indel file is empty, try-again!")
  # }
  
  open_file = tryCatch( 
    {  
      cat("\nreading files:", indel_file, snp_file)
      indel_file = read.table(indel_file, sep = "\t")
      snp_file = read.table(snp_file, sep = "\t")
    },
    error=function(error_message) {
      cat("\n")
      # message(error_message)
      stop("\nError: either snp/indel file is empty, try-again!")
    }
  )
  
  indel_info = indel_file[,8]
  indel_info = as.data.frame(indel_info)
  indel_data = indel_info %>% separate(indel_info, sep = ";", into = c(NA, NA, "r.start", "r.end", NA, "s.start", "s.end"))
  
  for(i in 1:ncol(indel_data)) {
    dat = indel_data %>% separate(colnames(indel_data)[i], sep = "=", into = c(NA, colnames(indel_data)[i]))
    indel_data[i] = dat[i]
  }
  
  indel_data = data.frame(apply(indel_data, 2, function(x) as.numeric(as.character(x))))
  
  indel_data$ref = (indel_data$r.end - indel_data$r.start)
  indel_data$simul = (indel_data$s.end - indel_data$s.start)
  indel_data$extend = abs(indel_data$simul - indel_data$ref) - 1
  indel_data$score = (indel_data$extend * extend ) + indel
  
  indel_score = sum(indel_data$score)
  cat("\nindel-score-needed:", need_indel_score)
  cat("\nindel-score:", indel_score)

  # window = round((window/100) * need_indel_score) # floor->round->nothing so that it isn't 0 for small values
  # cat("\nwindow:", window)
  # upper = need_indel_score + window
  # lower = need_indel_score - window
  # cat("\nindel score must be >=", lower, "and =<", upper)
  
  # if(indel_score < lower | indel_score > upper) { cat("\ntry-again!") }
  
  snp_count = nrow(snp_file)
  snp_score = snp_count * snp
  cat("\nsnp-score:", snp_score)
  
  actual_score = ((length - snp_count) * match) - snp_score - indel_score
  cat("\nactual-score:", actual_score)
  
  window = floor((window/100) * score)
  cat("\nwindow:", window)
  upper = score + window
  lower = score - window
  cat("\nactual score must be >=", lower, "and =<", upper)

  if(actual_score < lower | actual_score > upper) { cat("\ntry-again!") }
}
