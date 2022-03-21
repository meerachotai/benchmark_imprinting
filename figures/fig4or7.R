library(stringr)
library(ggplot2)
library(patchwork)

indir = "benchmark_sim_90/counts" # for fig. 7, change to benchmark_sim_95/counts/
infile = paste0(indir,"/anderson_AxB_1.txt")

anderson = read.table(infile, header = F, sep = "\t", stringsAsFactors = F, row.names = 1)
row.names(anderson) = str_remove(row.names(anderson), "strainA")

infile = paste0(indir,"/rep_1_1_strainAxstrainB_counts_merged.txt")
picard = read.table(infile, header = F, sep = "\t", stringsAsFactors = F, row.names = 1)

infile = paste0(indir,"/simul_counts+id_A.txt")
trueA = read.table(infile, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
trueA = trueA[,1, drop = FALSE] # keep just the first replicate

infile = paste0(indir,"/simul_counts+id_B.txt")
trueB = read.table(infile, header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
trueB = trueB[,1, drop = FALSE]

true = merge(trueA, trueB, by = 0) # AxB_A and AxB_B - first replicate
row.names(true) = true$Row.names
true = subset(true, select = -c(Row.names))
names(true) = c("mat.true", "pat.true")

all = merge(anderson, picard, by = 0, all = TRUE)
row.names(all) = all$Row.names
all = subset(all, select = -c(Row.names))

names(all) <- gsub(".x",paste0(".anderson"),names(all),fixed = T)
names(all) <- gsub(".y",paste0(".picard"),names(all),fixed = T)

names(all) <- gsub("V2",paste0("mat"),names(all),fixed = T)
names(all) <- gsub("V3",paste0("pat"),names(all),fixed = T)

all = merge(all, true, by = 0, all = TRUE)
all = subset(all, select = -c(Row.names))

# ---------------------- MATERNAL -----------------------------
pval = vector()
for (i in 1:nrow(all)) {
  y = tryCatch(
    expr = { chisq.test(c(all$mat.anderson[i], all$mat.true[i]), p = c(0.5, 0.5))$p.value },
    error = function(cond) { return(NA) }
  )
  pval = c(pval, y)
}
adjp = p.adjust(pval, method = "fdr", n = length(pval))
print(length(adjp[adjp < 0.05]))

sig = ifelse(adjp < 0.05 & !is.na(adjp), "maroon", "black")

one = ggplot(all, aes(x = mat.true, y = mat.anderson)) +
  geom_point(color = sig) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  xlab("True") + ylab("Anderson") +
  ggtitle("True vs. Anderson maternal counts")

pval = vector()
for (i in 1:nrow(all)) {
  y = tryCatch(
    expr = { chisq.test(c(all$mat.picard[i], all$mat.true[i]), p = c(0.5, 0.5))$p.value },
    error = function(cond) { return(NA) }
  )
  pval = c(pval, y)
}
adjp = p.adjust(pval, method = "fdr", n = length(pval))
print(length(adjp[adjp < 0.05]))

sig = ifelse(adjp < 0.05 & !is.na(adjp), "maroon", "black")

two = ggplot(all, aes(x = mat.true, y = mat.picard)) +
  geom_point(color = sig) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  xlab("True") + ylab("Picard") +
  ggtitle("True vs. Picard maternal counts")

wrap_plots(one, two, nrow = 1, ncol = 2)


# ---------------- PATERNAL -----------------------------
pval = vector()
for (i in 1:nrow(all)) {
  y = tryCatch(
    expr = { chisq.test(c(all$pat.anderson[i], all$pat.true[i]), p = c(0.5, 0.5))$p.value },
    error = function(cond) { return(NA) }
  )
  pval = c(pval, y)
}
adjp = p.adjust(pval, method = "fdr", n = length(pval))
print(length(adjp[adjp < 0.05])) # significantly different

sig = ifelse(adjp < 0.05 & !is.na(adjp), "maroon", "black")

three = ggplot(all, aes(x = pat.true, y = pat.anderson)) +
  geom_point(color = sig) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  xlab("True") + ylab("Anderson") +
  ggtitle("True vs. Anderson paternal counts")

pval = vector()
for (i in 1:nrow(all)) {
  y = tryCatch(
    expr = { chisq.test(c(all$pat.picard[i], all$pat.true[i]), p = c(0.5, 0.5))$p.value },
    error = function(cond) { return(NA) }
  )
  pval = c(pval, y)
}
adjp = p.adjust(pval, method = "fdr", n = length(pval))
print(length(adjp[adjp < 0.05]))

sig = ifelse(adjp < 0.05 & !is.na(adjp), "maroon", "black")

four = ggplot(all, aes(x = pat.true, y = pat.picard)) +
  geom_point(color = sig) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  xlab("True") + ylab("Picard") +
  ggtitle("True vs. Picard paternal counts")

wrap_plots(three, four, nrow = 1, ncol = 2)

