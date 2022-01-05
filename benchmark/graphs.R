library(ggplot2)

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

graph_maker <- function(parameter, parent, graph_header, true_neg, see_table) {
  # print(parent)
  method="Wyder"
  param = list.files(pattern = paste0(method,"*"), ignore.case = TRUE)
  
  table <- data.frame()
  for(i in 1:length(param)) {
    param_table = table_maker(parent, param[i], method, parameter)
    table = rbind(table, param_table)
  }
  
  method="Picard"
  
  param = list.files(pattern = paste0(method,"*"), ignore.case = TRUE)
  
  for(i in 1:length(param)) {
    param_table = table_maker(parent, param[i], method, parameter)
    table = rbind(table, param_table)
  }
  
  method="Anderson"
  
  param = list.files(pattern = paste0(method,"*"), ignore.case = TRUE)
  
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
  
  linetype = c("Anderson" = "dotted", "Picard" = "dashed", "Wyder" = "solid")
  table[is.na(table)] = 0
  
  if(see_table == TRUE) {
    print(table)
  }
  
  if(parameter == "dispersion") {
    table$dispersion = factor(table$dispersion, levels = c("low", "med", "high"))
    if(true_neg == TRUE) {
      color = c("true-pos" = "dodgerblue", "true-neg" = "darkorchid3")
      graph = ggplot(table) +
        geom_line(aes(x = dispersion, y = !!tp, color = "true-pos", linetype = method, group=method), show.legend = TRUE) +
        geom_line(aes(x = dispersion, y = !!tn, color = "true-neg", linetype = method, group=method), show.legend = TRUE) +
        scale_color_manual(name = "rates", values = color) +
        scale_linetype_manual(name = "method", values = linetype) +
        ylab("rate") + xlab(parameter) +
        ylim(0,1) +
        labs(subtitle = stitle) +
        ggtitle(paste("varying", graph_header))
    } else {
        color = c("true-pos" = "dodgerblue", "false-pos" = "indianred")
        graph = ggplot(table) +
          geom_line(aes(x = dispersion, y = !!tp, color = "true-pos", linetype = method, group=method), show.legend = TRUE) +
          geom_line(aes(x = dispersion, y = !!fp, color = "false-pos", linetype = method, group=method), show.legend = TRUE) +
          scale_color_manual(name = "rates", values = color) +
          scale_linetype_manual(name = "method", values = linetype) +
          ylab("rate") + xlab(parameter) +
          ylim(0,1) +
          labs(subtitle = stitle) +
          ggtitle(paste("varying", graph_header))
      }
  } else {
      if(true_neg == TRUE) {
        color = c("true-pos" = "dodgerblue", "true-neg" = "darkorchid3")
        graph = ggplot(table) +
          geom_line(aes(x = !!sym(parameter), y = !!tp, color = "true-pos", linetype = method), show.legend = TRUE) +
          geom_line(aes(x = !!sym(parameter), y = !!tn, color = "true-neg", linetype = method), show.legend = TRUE) +
          scale_color_manual(name = "rates", values = color) +
          scale_linetype_manual(name = "method", values = linetype) +
          ylab("rate") + xlab(parameter) +
          ylim(0,1) +
          labs(subtitle = stitle) +
          ggtitle(paste("varying", graph_header))
      } else {
        color = c("true-pos" = "dodgerblue", "false-pos" = "indianred")
        graph = ggplot(table) +
          geom_line(aes(x = !!sym(parameter), y = !!tp, color = "true-pos", linetype = method), show.legend = TRUE) +
          geom_line(aes(x = !!sym(parameter), y = !!fp, color = "false-pos", linetype = method), show.legend = TRUE) +  
          scale_color_manual(name = "rates", values = color) +
          scale_linetype_manual(name = "method", values = linetype) +
          ylab("rate") + xlab(parameter) +
          ylim(0,1) +
          labs(subtitle = stitle) +
          ggtitle(paste("varying", graph_header))
      }
  }
  graph
}

true_neg = TRUE
see_table = TRUE

setwd("~/Documents/SJ_Lab/Imprinting/benchmark_files/patbias")
parameter = "%paternal bias"
parent = "father"
graph_header = "%paternal bias"
graph_maker(parameter, parent, graph_header, true_neg, see_table)

setwd("~/Documents/SJ_Lab/Imprinting/benchmark_files/matbias")
parameter = "%maternal bias"
parent = "mother"
graph_header = "%maternal bias"
graph_maker(parameter, parent, graph_header, true_neg, see_table)

setwd("~/Documents/SJ_Lab/Imprinting/benchmark_files/sim_score")
parameter = "%similarity"
parent = "mother"
graph_header = "%similarity between strains"
graph_maker(parameter, parent, graph_header, true_neg, see_table)

setwd("~/Documents/SJ_Lab/Imprinting/benchmark_files/sim_score")
parameter = "%similarity"
parent = "father"
graph_header = "%similarity between strains"
graph_maker(parameter, parent, graph_header, true_neg, see_table)

setwd("~/Documents/SJ_Lab/Imprinting/benchmark_files/alpha")
parameter = "alpha"
parent = "mother"
graph_header = "alpha/p-value cutoffs"
graph_maker(parameter, parent, graph_header, true_neg, see_table)

setwd("~/Documents/SJ_Lab/Imprinting/benchmark_files/alpha")
parameter = "alpha"
parent = "father"
graph_header = "alpha/p-value cutoffs"
graph_maker(parameter, parent, graph_header, true_neg, see_table)

setwd("~/Documents/SJ_Lab/Imprinting/benchmark_files/disp")
parameter = "dispersion"
parent = "mother"
graph_header = "dispersion"
graph_maker(parameter, parent, graph_header, true_neg, see_table)

setwd("~/Documents/SJ_Lab/Imprinting/benchmark_files/disp")
parameter = "dispersion"
parent = "father"
graph_header = "dispersion"
graph_maker(parameter, parent, graph_header, true_neg, see_table)
