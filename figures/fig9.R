library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(stringr)

getTPM = function(dir, seqDepth, rep) {
  counts = list.files(path = dir,pattern = paste0( "_", seqDepth, "_1_counts.txt"))
  counts = read.table(paste0(dir,"/",counts), sep = "\t")
  tpm <- data.frame(matrix(ncol = rep, nrow = nrow(counts)))
  for (j in 1:rep) {
    file = list.files(path = dir, pattern = paste0("_", seqDepth, "_",j,"_counts.txt"))
    counts = read.table(paste0(dir,"/",file), sep = "\t")
    length = list.files(path = dir, pattern = paste0("_",seqDepth,".txt"))
    length = read.table(paste0(dir, "/",length), sep = "\t")
    counts = merge(counts, length, by = "V1")
    names(counts) = c("genes", "A", "B", "length")
    counts$sum = counts$A + counts$B
    
    counts$rpk = counts$sum / (counts$length/1000)
    total = sum(counts$rpk)
    per_mil = total/1000000
    tpm[j] = counts$rpk / per_mil
    names(tpm)[j] = paste0("AxB_", seqDepth, "_", j)
  }
  tpm$avg = rowMeans(tpm[1:rep])
  median(tpm$avg)
}

# functions from benchmarking graphs
find_caption = function(caption, extract) {
  caption = stringr::str_extract_all(caption, "\\S+")
  as.numeric(caption[[1]][which(grepl(paste0(extract,":"), caption[[1]])) + 1])
}

roc_stats <- function(param_table, n, n_megs, n_pegs) {
  param_table$fp_mat = param_table$mat - param_table$tp_mat
  param_table$fp_pat = param_table$pat - param_table$tp_pat
  
  param_table$tp_mat_rate = (param_table$tp_mat / n_megs)
  param_table$tp_pat_rate = (param_table$tp_pat / n_pegs)
  
  param_table$tn_mat_rate = 1 - (param_table$fp_mat / (n - n_megs))
  param_table$tn_pat_rate = 1 - (param_table$fp_pat / (n - n_pegs))
  
  param_table$fp_mat_rate = (param_table$fp_mat / (n - n_megs))
  param_table$fp_pat_rate = (param_table$fp_pat / (n - n_pegs))
  
  param_table
}

table_maker <- function(parent, params, method, parameter, dir, rep) {
  
  all_combined = readLines(params)
  caption = all_combined[1]
  param_table = read.table(textConnection(all_combined[-1]), header = T, sep = "\t")
  
  n = find_caption(caption, "n")
  n_megs = find_caption(caption, "nmegs")
  n_pegs = find_caption(caption, "npegs")
  
  param_table[is.na(param_table)] = 0
  names(param_table)[1] = parameter
  
  param_table = roc_stats(param_table, n, n_megs, n_pegs)
  param_table$method =  rep(method, nrow(param_table))
  
  tpm = vector()
  for (i in param_table$`sequencing depth`) {
    x = getTPM(dir, i, rep)
    tpm = c(tpm, x)
  }
  
  param_table$tpm = tpm
  param_table
}



# DATA ANALYSIS -------------------------------------

setwd("~/Documents/SJ_Lab/Imprinting/Winter2022/benchmark_files/feb20/seq_depth")
parameter = "sequencing depth (TPM)"
graph_header = "sequencing depth"
file_label = "seq_depth"
rep = 3
dir = "~/Documents/SJ_Lab/Imprinting/Winter2022/benchmark_files/feb20/seqDepth_counts"
table <- data.frame()

method="Wyder"
param = list.files(pattern = paste0(method, "_", file_label), ignore.case = TRUE)
param_table = table_maker(parent, param[1], method, parameter, dir, rep)
param_table = param_table[1:(nrow(param_table))-1,]
table = rbind(table, param_table)

method="Picard"
param = list.files(pattern = paste0(method,"_", file_label), ignore.case = TRUE)
param_table = table_maker(parent, param[i], method, parameter, dir, rep)
param_table = param_table[1:(nrow(param_table))-1,]
table = rbind(table, param_table)

method="Anderson"
param = list.files(pattern = paste0(method,"_", file_label), ignore.case = TRUE)
param_table = table_maker(parent, param[i], method, parameter, dir, rep)
param_table = param_table[1:(nrow(param_table))-1,]
table = rbind(table, param_table)

method="Anderson_Picard_Combined"
param = list.files(pattern = paste0(method,"_"), ignore.case = TRUE)
param_table = table_maker(parent, param[i], method, parameter, dir, rep)
param_table = param_table[1:(nrow(param_table))-1,]
table = rbind(table, param_table)

method="Anderson_Wyder_Combined"
param = list.files(pattern = paste0(method,"_"), ignore.case = TRUE)
param_table = table_maker(parent, param[i], method, parameter, dir, rep)
param_table = param_table[1:(nrow(param_table))-1,]
table = rbind(table, param_table)

method="Roth_Wyder_Combined"
param = list.files(pattern = paste0(method,"_"), ignore.case = TRUE)
param_table = table_maker(parent, param[i], method, parameter, dir, rep)
param_table = param_table[1:(nrow(param_table))-1,]
table = rbind(table, param_table)

method="Roth"
param = list.files(pattern = paste0(method,"_", file_label), ignore.case = TRUE)
param_table = table_maker(parent, param[i], method, parameter, dir, rep)
param_table = param_table[1:(nrow(param_table))-1,]
table = rbind(table, param_table)

# mother / MEGs
tp = sym("tp_mat_rate")
fp = sym("fp_mat_rate")
stitle = "showing true-positive/false-positive rates for MEGs"

linetype = c("true-pos" = "solid", "false-pos" = "dashed")
graph1 = ggplot(table) +
  geom_line(aes(x = tpm, y = !!tp, color = method, linetype = "true-pos"), show.legend = TRUE) +
  geom_line(aes(x = tpm, y = !!fp, color = method, linetype = "false-pos"), show.legend = TRUE) +  
  # scale_color_manual(name = "method", values = color) +
  scale_color_brewer(palette = "Paired") +
  scale_linetype_manual(name = "rates", values = linetype) +
  ylab("rate") + xlab(parameter) +
  ylim(0,1) + 
  labs(subtitle = stitle) +
  ggtitle(paste("varying", graph_header))

# father / PEGs
tp = sym("tp_pat_rate")
fp = sym("fp_pat_rate")
stitle = "showing true-positive/false-positive rates for PEGs"

linetype = c("true-pos" = "solid", "false-pos" = "dashed")
graph2 = ggplot(table) +
  geom_line(aes(x = tpm, y = !!tp, color = method, linetype = "true-pos"), show.legend = TRUE) +
  geom_line(aes(x = tpm, y = !!fp, color = method, linetype = "false-pos"), show.legend = TRUE) +  
  # scale_color_manual(name = "method", values = color) +
  scale_color_brewer(palette = "Paired") +
  scale_linetype_manual(name = "rates", values = linetype) +
  ylab("rate") + xlab(parameter) +
  ylim(0,1) + 
  labs(subtitle = stitle) +
  ggtitle(paste("varying", graph_header))

wrap_plots(graph1, graph2, ncol = 2)
ggsave(paste0(indir, "/",file_name, "_graph.png"), width = 13, height = 6, dpi = 150, units = "in", device='png')

write.table(table, "seq_depth.txt", quote = F, sep = "\t", row.names = F)
