library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(stringr)

indir = "benchmark_sim_90"

anderson = read.table(paste0(indir, "/anderson_stats.txt"), sep = "\t", header = T, row.names = 1)
row.names(anderson) = str_remove(row.names(anderson), "strainA")

true_megs = scan(paste0(indir,"/true_MEGs.txt"), what = "character")
true_pegs = scan(paste0(indir,"/true_PEGs.txt"), what = "character")

is_true_meg = row.names(anderson) %in% true_megs
is_true_peg = row.names(anderson) %in% true_pegs

anderson_megs = row.names(anderson[(anderson$imprint.strainA == "MEG" & anderson$imprint.strainB == "MEG"),])
anderson_pegs = row.names(anderson[(anderson$imprint.strainA == "PEG" & anderson$imprint.strainB == "PEG"),])


# find false negatives
fn_megs = anderson[(anderson$imprint.strainA == "no.imprint" | 
                      anderson$imprint.strainB == "no.imprint") & is_true_meg,]

fn_pegs = anderson[(anderson$imprint.strainA == "no.imprint" | 
                      anderson$imprint.strainB == "no.imprint") & is_true_peg,]

# MEGs --------------------------------------

print(nrow(fn_megs))
log2fc_megs = fn_megs[fn_megs$log2fc.strainA == FALSE | fn_megs$log2fc.strainB == FALSE,]
print(nrow(log2fc_megs))
mat_megs = fn_megs[which(fn_megs$matPref_met.strainA == FALSE | fn_megs$matPref_met.strainB == FALSE),]
print(nrow(mat_megs))
fdr_megs = fn_megs[which(fn_megs$fdr_met.strainA == FALSE | fn_megs$fdr_met.strainB == FALSE),]
print(nrow(fdr_megs))

fdr_megs = row.names(fdr_megs)
mat_megs = row.names(mat_megs)

myCol <- brewer.pal(4, "Pastel2")

venn.plot = venn.diagram(x = list(true_megs, anderson_megs, fdr_megs, mat_megs), 
                         category.names = c("True MEGs","Anderson MEGs","not met FDR", "not met\nm:p ratio"), fill = myCol,
                         filename = NULL, cat.fontfamily = "sans", cat.cex = 1, cex = 1.5,#print.mode="percent",
                         fontfamily = "sans")

ggsave("fig3_MEGs.png", venn.plot, device = "png", width = 7.5)

# PEGs --------------------------------------

print(nrow(fn_pegs))
log2fc_pegs = fn_pegs[fn_pegs$log2fc.strainA == FALSE | fn_pegs$log2fc.strainB == FALSE,]
print(nrow(log2fc_pegs))
mat_pegs = fn_pegs[which(fn_pegs$matPref_met.strainA == FALSE | fn_pegs$matPref_met.strainB == FALSE),]
print(nrow(mat_pegs))
fdr_pegs = fn_pegs[which(fn_pegs$fdr_met.strainA == FALSE | fn_pegs$fdr_met.strainB == FALSE),]
print(nrow(fdr_pegs))

fdr_pegs = row.names(fdr_pegs)
mat_pegs = row.names(mat_pegs)

myCol <- brewer.pal(4, "Pastel2")

venn.plot = venn.diagram(x = list(true_pegs, anderson_pegs, fdr_pegs, mat_pegs), 
                         category.names = c("True PEGs","Anderson PEGs","not met FDR", "not met\nm:p ratio"), fill = myCol,
                         filename = NULL, cat.fontfamily = "sans", cat.cex = 1, cex = 1.5, #print.mode="percent",
                         fontfamily = "sans")

ggsave("fig3_PEGs.png", venn.plot, device = "png", width = 7.5)
