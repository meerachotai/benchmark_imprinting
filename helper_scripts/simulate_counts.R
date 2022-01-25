#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))

parser = ArgumentParser()

parser$add_argument("--seed", type = "integer", help = "choose a seed for reproducibility")
parser$add_argument("--disp", type="character", default = "med", help = "negative binomial dist. size (inverse to dispersion), default = 100")
parser$add_argument("--n", type = "integer", default = 15000, help = "number of non-biased 'genes' to generate, default = 15000")
parser$add_argument("--MEGbias", type = "double", default = 80, help = "% maternal bias for simulated 'MEGs' default = 80")
parser$add_argument("--PEGbias", type = "double", default = 50, help = "% maternal bias for simulated 'PEGs' default = 50")
parser$add_argument("--nMEG", type = "integer", default = 200, help = "number of maternally-biased 'genes' to generate, default = 200")
parser$add_argument("--nPEG", type = "integer", default = 200, help = "number of paternally-biased 'genes' to generate, default = 200")
parser$add_argument("--rep", type = "integer", default = 1, help = "number of replicates, default = 1")
parser$add_argument("--seqDepth", type = "integer", default = 1, help = "sequencing depth, default = 1")

parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

if(!is.null(opt$seed)) {
  seed = opt$seed
  set.seed(seed)
} else {seed = "NA"}

###################
#### functions ####
###################
#           5%          50%          95% 
# 0.0004353327 0.0089020807 0.1101403767 

size_func <- function(disp_type, total_exp) {
  if(disp_type == "low") {
    size = 0.11 * total_exp
  }
  if(disp_type == "med") {
    size = 0.008 * total_exp
  }
  if(disp_type == "high") {
    size = 0.0004 * total_exp
  }
  if(size == 0) {
    size = 0.01
  }
  return(size)
}

generate_totexp <- function(mu_mat, mu_pat, disp_type, n) {
  tot_exp = mu_mat + mu_pat
  size = size_func(disp_type, tot_exp)
  rnbinom(n = n, size = size, mu = tot_exp)
}

nbin_sum_r <- function(mat1_ratio, mat2_ratio, pat_ratio, size, total_exp, n) {
  tot = mat1_ratio + mat2_ratio + pat_ratio
  mu_mat1 = (total_exp * mat1_ratio/tot)
  mu_mat2 = (total_exp * mat2_ratio/tot)
  mu_pat = (total_exp * pat_ratio/tot)
  
  mat1 = rnbinom(n = n, mu = mu_mat1, size = size)
  mat2 = rnbinom(n = n, mu = mu_mat2, size = size) 
  pat = rnbinom(n = n, mu = mu_pat, size = size)
  
  mat = data.frame(mat1, mat2)
  mat = rowSums(mat)
  
  data.frame(mat,pat)
}

nbin_sum_p <- function(percent_mat, size, total_exp, n) {
  percent_mat = percent_mat / 100
  percent_mat = percent_mat / 2
  
  mu_mat = percent_mat * total_exp
  mu_pat = total_exp - (mu_mat * 2)
  
  mat1 = rnbinom(n = n, mu = mu_mat, size = size)
  mat2 = rnbinom(n = n, mu = mu_mat, size = size)
  pat = rnbinom(n = n, mu = mu_pat, size = size)
  
  mat = data.frame(mat1, mat2)
  mat = rowSums(mat)
  
  data.frame(mat,pat)
}

nbin_loop_r <- function(mat1_ratio, mat2_ratio, pat_ratio, disp_type, vec_totexp, n) {
  AxB = data.frame(matrix(NA, nrow = n, ncol = 2))
  BxA = data.frame(matrix(NA, nrow = n, ncol = 2))
  colnames(AxB) = c("mat", "pat")
  colnames(BxA) = colnames(AxB)
  
  for (i in 1:n) {
    size = size_func(disp_type, vec_totexp[i])
    # n = 2, one for each cross direction
    counts = nbin_sum_r(mat1_ratio, mat2_ratio, pat_ratio,size, vec_totexp[i], 2)
    AxB[i,] = counts[1,]
    BxA[i,] = counts[2,]
  }
  data.frame(AxB, BxA)
}

nbin_loop_p = function(percent_mat, disp_type, vec_totexp, n) {
  AxB = data.frame(matrix(NA, nrow = n, ncol = 2))
  BxA = data.frame(matrix(NA, nrow = n, ncol = 2))
  colnames(AxB) = c("mat", "pat")
  colnames(BxA) = colnames(AxB)
  
  for (i in 1:n) {
    size = size_func(disp_type, vec_totexp[i])
    # n = 2, one for each cross direction
    counts = nbin_sum_p(percent_mat,size, vec_totexp[i], 2)
    AxB[i,] = counts[1,]
    BxA[i,] = counts[2,]
  }
  data.frame(AxB, BxA)
}

###############
### options ###
###############

outprefix = opt$outprefix
n = opt$n
megs_bias = opt$MEGbias
pegs_bias = opt$PEGbias
n_megs = opt$nMEG
n_pegs = opt$nPEG
disp_type = opt$disp
rep = opt$rep
seqDepth = opt$seqDepth

# size.m        mu.m      size.p        mu.p      size.t        mu.t       ratio 
# 0.1805655 239.3698305   0.1821203 119.7124271   0.1799987 336.5474092   2.1142447

cat("Sequencing Depth:", seqDepth, "\n")
mu_mat = 239.37 * seqDepth
mu_pat = 119.71 * seqDepth

cat("\nRunning script type 3: total expression generated  \non a negative binomial dist. for each gene\n")

cat("\n--------------------------\n")
cat("seed:", seed, "\n")
cat("dispersion:",disp_type,"\n")
cat("total number of 'genes':", n+n_megs+n_pegs, "\n") 
cat("bias: %maternal bias for maternally-biased 'genes' -",megs_bias,"%\n")
cat("bias: %maternal bias for paternally-biased 'genes' -",pegs_bias,"%\n")
cat("--------------------------\n")

#################
### main code ###
#################

# randomized total expression for each gene in a loop ------------

vec_totexp = generate_totexp(mu_mat, mu_pat, disp_type, n)
vec_totexp_m = generate_totexp(mu_mat, mu_pat, disp_type, n_megs)
vec_totexp_p = generate_totexp(mu_mat, mu_pat, disp_type, n_pegs)

for (i in 1:rep) {
  # counts for unbiased genes first
  unbiased = nbin_loop_r(1,1,1, disp_type, vec_totexp, n)
  # counts for MEGs and PEGs
  true_MEGs = nbin_loop_p(megs_bias, disp_type,vec_totexp_m, n_megs)
  true_PEGs = nbin_loop_p(pegs_bias, disp_type, vec_totexp_p, n_pegs)
  
  counts = rbind(unbiased, true_MEGs, true_PEGs)
  # separate both directions' counts
  counts_AxB = data.frame(counts[,1], counts[,2])
  counts_BxA = data.frame(counts[,3], counts[,4])
  colnames(counts_AxB) = c("mat", "pat")
  colnames(counts_BxA) = colnames(counts_AxB)
  counts_BxA = counts_BxA[,c("pat", "mat")] # so that A is always first column
  
  if(i == 1) { 
    all_counts_AxB = counts_AxB 
    all_counts_BxA = counts_BxA  
  }
  else { 
    # format of all_counts:
    # A-rep1 B-rep1 A-rep2 B-rep2 ...
    all_counts_AxB = cbind(all_counts_AxB, counts_AxB)
    all_counts_BxA = cbind(all_counts_BxA, counts_BxA)
  }
}


cat("\nWriting count tables for both directions...")
write.table(all_counts_AxB,paste0(outprefix,"_AxB.txt"), row.names=T, quote=F, col.names=F, sep="\t")
write.table(all_counts_BxA,paste0(outprefix,"_BxA.txt"), row.names=T, quote=F, col.names=F, sep = "\t")
cat(" Done.\n\n")

cat("\nWriting MEG and PEG 'true' lists from simulated counts file...")
megs = row.names(counts_AxB[(n+1):(n+1+n_megs),])
pegs = row.names(counts_AxB[(n+n_megs+1):(n+n_megs+n_pegs),])
write.table(megs, paste0(outprefix,"_megs.txt"), row.names=F, quote=F, col.names=F)
write.table(pegs, paste0(outprefix,"_pegs.txt"), row.names=F, quote=F, col.names=F)
cat(" Done.\n")
