library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)

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

table_maker <- function(parent, params, method, parameter) {
  
  all_combined = readLines(params)
  caption = all_combined[1]
  param_table = read.table(textConnection(all_combined[-1]), header = T, sep = "\t")
  
  n = find_caption(caption, "n")
  n_megs = find_caption(caption, "nmegs")
  n_pegs = find_caption(caption, "npegs")
  # print(c(n_megs, n_pegs))
  
  param_table[is.na(param_table)] = 0
  names(param_table)[1] = parameter
  
  param_table = roc_stats(param_table, n, n_megs, n_pegs)
  param_table$method =  rep(method, nrow(param_table))
  
  param_table
}

graph_maker <- function(parameter, parent, graph_header, true_neg, see_table, file_label) {
  # print(parent)
  method="Wyder"
  param = list.files(pattern = paste0(method, "_", file_label), ignore.case = TRUE)
  
  table <- data.frame()
  for(i in 1:length(param)) {
    param_table = table_maker(parent, param[i], method, parameter)
    table = rbind(table, param_table)
  }
  
  method="Picard"
  
  param = list.files(pattern = paste0(method,"_", file_label), ignore.case = TRUE)
  
  for(i in 1:length(param)) {
    param_table = table_maker(parent, param[i], method, parameter)
    table = rbind(table, param_table)
  }
  
  method="Anderson"
  
  param = list.files(pattern = paste0(method,"_", file_label), ignore.case = TRUE)
  
  for(i in 1:length(param)) {
    param_table = table_maker(parent, param[i], method, parameter)
    table = rbind(table, param_table)
  }
  
  method="Anderson_Picard_Combined"
  
  param = list.files(pattern = paste0(method,"_"), ignore.case = TRUE)
  
  for(i in 1:length(param)) {
    param_table = table_maker(parent, param[i], method, parameter)
    table = rbind(table, param_table)
  }
  
  method="Anderson_Wyder_Combined"
  
  param = list.files(pattern = paste0(method,"_"), ignore.case = TRUE)
  
  for(i in 1:length(param)) {
    param_table = table_maker(parent, param[i], method, parameter)
    table = rbind(table, param_table)
  }
  
  method="Roth_Wyder_Combined"
  
  param = list.files(pattern = paste0(method,"_"), ignore.case = TRUE)
  
  for(i in 1:length(param)) {
    param_table = table_maker(parent, param[i], method, parameter)
    table = rbind(table, param_table)
  }
  
  method="Roth"
  
  param = list.files(pattern = paste0(method,"_", file_label), ignore.case = TRUE)
  
  for(i in 1:length(param)) {
    param_table = table_maker(parent, param[i], method, parameter)
    table = rbind(table, param_table)
  }
  
  if(parent == "mother") {
    tp = sym("tp_mat_rate")
    fp = sym("fp_mat_rate")
    tn = sym("tn_mat_rate")
    if(true_neg == TRUE) {
      stitle = "showing true-positive/true-negative rates for MEGs"
    } else {
      stitle = "showing true-positive/false-positive rates for MEGs"
    }
  } else {
    tp = sym("tp_pat_rate")
    fp = sym("fp_pat_rate")
    tn = sym("tn_pat_rate")
    if(true_neg == TRUE) {
      stitle = "showing true-positive/true-negative rates for PEGs"
    } else {
      stitle = "showing true-positive/false-positive rates for PEGs"
    }
  }
  
  table[is.na(table)] = 0
  
  if(see_table == TRUE) {
    print(table)
  }
  
  if(parameter == "dispersion") {
    table$dispersion = factor(table$dispersion, levels = c("low", "med", "high"))
    if(true_neg == TRUE) {
      linetype = c("true-pos" = "solid", "true-neg" = "dashed")
      graph = ggplot(table) +
        geom_line(aes(x = dispersion, y = !!tp, color = method, linetype = "true-pos", group=method), show.legend = TRUE) +
        geom_line(aes(x = dispersion, y = !!tn, color = method, linetype = "true-neg", group=method), show.legend = TRUE) +
        # scale_color_manual(name = "method", values = color) +
        scale_color_brewer(palette = "Paired") +
        scale_linetype_manual(name = "rates", values = linetype) +
        ylab("rate") + xlab(parameter) +
        ylim(0,1) +
        labs(subtitle = stitle) +
        ggtitle(paste("varying", graph_header))
    } else {
      linetype = c("true-pos" = "solid", "false-pos" = "dashed")
      graph = ggplot(table) +
        geom_line(aes(x = dispersion, y = !!tp, color = method, linetype = "true-pos", group=method), show.legend = TRUE) +
        geom_line(aes(x = dispersion, y = !!fp, color = method, linetype = "false-pos", group=method), show.legend = TRUE) +
        # scale_color_manual(name = "method", values = color) +
        scale_color_brewer(palette = "Paired") +
        scale_linetype_manual(name = "rates", values = linetype) +
        ylab("rate") + xlab(parameter) +
        ylim(0,1) +
        labs(subtitle = stitle) +
        ggtitle(paste("varying", graph_header))
    }
  } else {
    if(true_neg == TRUE) {
      linetype = c("true-pos" = "solid", "true-neg" = "dashed")
      graph = ggplot(table) +
        geom_line(aes(x = !!sym(parameter), y = !!tp, color = method, linetype = "true-pos"), show.legend = TRUE) +
        geom_line(aes(x = !!sym(parameter), y = !!tn, color = method, linetype = "true-neg"), show.legend = TRUE) +
        scale_color_brewer(palette = "Paired") +
        scale_linetype_manual(name = "rate", values = linetype) +
        ylab("rate") + xlab(parameter) +
        ylim(0,1) +
        labs(subtitle = stitle) +
        ggtitle(paste("varying", graph_header))
    } else {
      linetype = c("true-pos" = "solid", "false-pos" = "dashed")
      graph = ggplot(table) +
        geom_line(aes(x = !!sym(parameter), y = !!tp, color = method, linetype = "true-pos"), show.legend = TRUE) +
        geom_line(aes(x = !!sym(parameter), y = !!fp, color = method, linetype = "false-pos"), show.legend = TRUE) +  
        scale_color_brewer(palette = "Paired") +
        scale_linetype_manual(name = "rates", values = linetype) +
        ylab("rate") + xlab(parameter) +
        ylim(0,1) +
        labs(subtitle = stitle) +
        ggtitle(paste("varying", graph_header))
    }
  }
  graph
}

true_neg = FALSE
see_table = TRUE
indir = "benchmark_files"

# ------------------------ MATERNAL BIAS -----------------------------------
setwd(paste0(indir,"/patbias"))
file_label = "patbias"
parameter = "%maternal bias"
parent = "father"
graph_header = "%maternal bias"
graph2 = graph_maker(parameter, parent, graph_header, true_neg, see_table, file_label)

setwd(paste0(indir,"/matbias"))
file_label = "matbias"
parameter = "%maternal bias"
parent = "mother"
graph_header = "%maternal bias"
graph1 = graph_maker(parameter, parent, graph_header, true_neg, see_table, file_label)

wrap_plots(graph1, graph2, ncol = 2)
ggsave(paste0(file_label, "_graph.png"), width = 13, height = 6, dpi = 150, units = "in", device='png')

# ------------------------- SIMILARITY SCORES -----------------------------
setwd(paste0(indir,"/sim_score"))
file_label = "sim_score"
parameter = "%similarity"
parent = "mother"
graph_header = "%similarity between strains"
graph1 = graph_maker(parameter, parent, graph_header, true_neg, see_table, file_label)

setwd(paste0(indir,"/sim_score"))
file_label = "sim_score"
parameter = "%similarity"
parent = "father"
graph_header = "%similarity between strains"
graph2 = graph_maker(parameter, parent, graph_header, true_neg, see_table, file_label)

wrap_plots(graph1, graph2, ncol = 2)
ggsave(paste0(file_label, "_graph.png"), width = 13, height = 6, dpi = 150, units = "in", device='png')

# -------------------------- P-VALUE ------------------------------------
setwd(paste0(indir,"/alpha"))
file_label = "alpha"
parameter = "p-value cutoffs"
parent = "mother"
graph_header = "p-value cutoffs"
graph1 = graph_maker(parameter, parent, graph_header, true_neg, see_table, file_label)

setwd(paste0(indir,"/alpha"))
file_label = "alpha"
parameter = "p-value cutoffs"
parent = "father"
graph_header = "p-value cutoffs"
graph2 = graph_maker(parameter, parent, graph_header, true_neg, see_table, file_label)

wrap_plots(graph1, graph2, ncol = 2)
ggsave(paste0(file_label, "_graph.png"), width = 13, height = 6, dpi = 150, units = "in", device='png')

# ---------------------------- DISPERSION ------------------------------------------

setwd(paste0(indir,"/disp"))
file_label = "disp"
parameter = "dispersion"
parent = "mother"
graph_header = "dispersion"
graph1 = graph_maker(parameter, parent, graph_header, true_neg, see_table, file_label)

setwd(paste0(indir,"/disp"))
file_label = "disp"
parameter = "dispersion"
parent = "father"
graph_header = "dispersion"
graph2 = graph_maker(parameter, parent, graph_header, true_neg, see_table, file_label)

wrap_plots(graph1, graph2, ncol = 2)
ggsave(paste0(file_label, "_graph.png"), width = 13, height = 6, dpi = 150, units = "in", device='png')
