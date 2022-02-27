# note: matbias = 90, patbias = 20, disp = med, similarity = 90
setwd("~/Documents/SJ_Lab/Imprinting/Winter2022/benchmark_sim_90")

anderson= read.table("anderson_stats.txt", sep = "\t", header = T, row.names = 1)

true_megs = scan("true_MEGs.txt", what = "character")
true_pegs = scan("true_PEGs.txt", what = "character")

# true_megs = true_megs[-1]
true = c(true_megs, true_pegs)
true = unique(true)


library(stringr)
row.names(anderson) = str_remove(row.names(anderson), "strainA")

# both directional imprinted
imprinted = anderson[(anderson$imprint.strainA == "MEG" & anderson$imprint.strainB == "MEG") | 
                       (anderson$imprint.strainA == "PEG" & anderson$imprint.strainB == "PEG"),]

not_imprint = anderson[((anderson$imprint.strainA == "no.imprint" | 
                           anderson$imprint.strainB == "no.imprint") & row.names(anderson) %in% true),]


print(nrow(not_imprint))
log2fc = not_imprint[not_imprint$log2fc.strainA == FALSE | not_imprint$log2fc.strainB == FALSE,]
print(nrow(log2fc))
mat = not_imprint[which(not_imprint$matPref_met.strainA == FALSE | not_imprint$matPref_met.strainB == FALSE),]
print(nrow(mat))
fdr = not_imprint[which(not_imprint$fdr_met.strainA == FALSE | not_imprint$fdr_met.strainB == FALSE),]
print(nrow(fdr))

fdr_genes = row.names(fdr)
mat_genes = row.names(mat)
imprint_genes = row.names(imprinted)

library(VennDiagram)
library(RColorBrewer)
library(ggplot2)

myCol <- brewer.pal(4, "Pastel2")

venn.plot = venn.diagram(x = list(true, imprint_genes,fdr_genes, mat_genes), 
                         category.names = c("True", "Anderson","not met FDR", "not met\nm:p ratio"), fill = myCol,
                         filename = NULL, cat.fontfamily = "sans", cat.cex = 1.5, cex = 1.75,  print.mode="percent",
                         fontfamily = "sans")

ggsave("compare_calls.png", venn.plot, device = "png", width = 6.5)
# png("compare_calls.png");
# grid.draw(venn.plot);
# dev.off();
