#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))
if (is.element('tidyverse', installed.packages()[,1])==FALSE) { stop("Error: R package tidyverse could not be loaded") }
suppressPackageStartupMessages(library(tidyverse))

parser = ArgumentParser()
parser$add_argument("-a", type = "character", default = "", help = "AxB counts file (simulated)")
parser$add_argument("-b", type = "character", default = "", help = "BxA counts file (simulated)")
parser$add_argument("-A", type = "character", default = "", help = "gene_id + sequences for A")
parser$add_argument("-B", type = "character", default = "", help = "gene_id + sequences for B")
parser$add_argument("-p", type = "character", default = 'I', help = "ASCII code for PHRED+33 score, default = I")
parser$add_argument("-r", type="double", default = 50, help = "read length, default = 50")
parser$add_argument("-s", type="double", default = 5, help = "set seed for reproducibility")
parser$add_argument("-R", type = "integer", default = 1, help = "number of replicates")
parser$add_argument("outprefix", nargs=1, help = "prefix to use for all output files")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

counts_AxB = opt$a
counts_BxA = opt$b
seq_A = opt$A
seq_B = opt$B

read_length = opt$r
seed = opt$s
phred = opt$p
outprefix = opt$outprefix

replicates = opt$R

# library(tidyverse)
# read_length = 50
# seq_A = "colA_seq.txt"
# seed = 5
# counts_AxB = "count_simul_AxB.txt"
# counts_BxA = "count_simul_BxA.txt"
# phred = 'I'

set.seed(seed)

counts_AxB = read.table(counts_AxB, sep = "\t", header=F, row.names = 1)
counts_BxA = read.table(counts_BxA, sep = "\t", header=F, row.names = 1)

seq_A = read.table(seq_A, sep = "\t", colClasses = "character", header=F)
seq_B = read.table(seq_B, sep = "\t", colClasses = "character", header=F)

phred_vec = paste(rep(phred, read_length), collapse = "")

if(nrow(counts_AxB) != nrow(seq_A) & nrow(counts_BxA) != nrow(seq_A)) { # made this mistake, need fail-check
  stop("Error: number of gene_ids for A do not match number of counts simulated, check again")
}

if(nrow(counts_AxB) != nrow(seq_B) & nrow(counts_BxA) != nrow(seq_B)) { # made this mistake, need fail-check
  stop("Error: number of gene_ids for B do not match number of counts simulated, check again")
}

# 1. randomly pick a start index
# 2. find the next -read_length- nucleotides in the sequence
# 3. do this -counts- times (based on the counts table)

reads_simul = function(counts, seq_table, out, read_length, phred_vec) {
  # simulate reads
  all_reads = vector()
  all_labels = vector()
  for(i in 1:length(counts)) {
    seq = str_split(seq_table[i,2], pattern = "")[[1]] # breaks the seq into vector if bases
    if(length(seq) > read_length) { # if the read has more than -read_length- bases
      index_range = length(seq) - read_length # so that the end is not too short
      reads = vector()
      labels = rep(gsub(" ", "", seq_table[i,1]),counts[i]) # adds gene label (no spaces)
      for(j in 1:counts[i]) {
        # TODO: improve this random sampling, for now it's uniform dist, rounded down
        index = floor(runif(1,1,index_range))
        read = paste(seq[index:(index+(read_length-1))], collapse ="") # updated it so that it does not give (read-length + 1) bp
        reads = c(reads, read)
      }
      all_reads = c(all_reads, reads)
      all_labels = c(all_labels, labels)
    }
  }
  
  # save file
  sink(paste0(out, ".fq"))
  for(i in 1:length(all_reads)) {
    cat(paste0("@",str_split(all_labels[i], '>')[[1]][2]), "\n")
    cat(paste(all_reads[i], "+", phred_vec, sep = "\n"), "\n")
  }
  sink()
  
  # all_labels
}

for (rep in 1:replicates) {
  # changed this to accept replicates
  A_col = ((rep - 1) * 2) + 1
  B_col = ((rep - 1) * 2) + 2
  
  counts_AxB_A = counts_AxB[,A_col]
  counts_BxA_A = counts_BxA[,A_col]
  counts_AxB_B = counts_AxB[,B_col]
  counts_BxA_B = counts_BxA[,B_col]
  
  # updated these file names to fit the map function in simulate_reads.sh
  reads_simul(counts_AxB_A, seq_A, paste0(outprefix,"_AxB_", rep, "_A"), read_length, phred_vec)
  reads_simul(counts_BxA_A, seq_A, paste0(outprefix,"_BxA_", rep, "_A"), read_length, phred_vec)
  reads_simul(counts_AxB_B, seq_B, paste0(outprefix,"_AxB_", rep,"_B"), read_length, phred_vec)
  reads_simul(counts_BxA_B, seq_B, paste0(outprefix,"_BxA_", rep, "_B"), read_length, phred_vec)
}

# seq_A$seq = substring(unlist(strsplit(seq_A[,1], " "))[seq(1,length(unlist(strsplit(seq_A[,1], " "))), 3)],2)
counts_id_A = cbind(seq_A[,1], counts_AxB[,seq(1, (replicates *2), 2)], counts_BxA[,seq(1, (replicates *2), 2)])
colnames(counts_id_A) = c("seq_id",rep("AxB_A",replicates), rep("BxA_A", replicates))
write.table(counts_id_A, paste0(outprefix, "_counts+id_A.txt"), col.names = T, row.names = F, sep = "\t", quote = F)

# seq_B$seq = substring(unlist(strsplit(seq_B[,1], " "))[seq(1,length(unlist(strsplit(seq_B[,1], " "))), 3)],2)
counts_id_B = cbind(seq_B[,1], counts_AxB[,seq(2, (replicates*2), 2)], counts_BxA[,seq(2, (replicates*2), 2)])
colnames(counts_id_B) = c("seq_id",rep("AxB_B",replicates), rep("BxA_B", replicates))
write.table(counts_id_B, paste0(outprefix, "_counts+id_B.txt"), col.names = T, row.names = F, sep = "\t", quote = F)

# ----------------- to double-check --------------------

# all_labels_AxB_A = reads_simul(counts_AxB_A, seq_A, paste0(outprefix,"_AxB_A"), read_length, phred_vec)
# 
# counts_check_AxB_A = vector()
# for(i in 1:nrow(seq_A)) { counts_check_AxB_A[i] = length(grep(seq_A[i,1], all_labels_AxB_A)) }
# inspect = as.data.frame(cbind(counts_AxB_A, counts_check_AxB_A))
# inspect[inspect[,1] != inspect[,2], ]

