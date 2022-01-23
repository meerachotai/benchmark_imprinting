outcounts = "roth_counts.txt"
labels = "rename.txt"
A = "strainA"
B = "strainB"
rep = 3
# paired = TRUE

i = 1
infile = paste0("rep_", i, "_", i, "_", A, "x", B, "_snp_report.bed") 

data = read.table(infile, sep = "\t")
data = subset(data, select = -c(V4))
names(data) = c("genes", "start", "end", paste0('AxB_',i,'_A'), paste0('AxB_',i,'_B'))

genes = unique(data[,1])
snpNum = vector(mode="numeric", length=length(genes))
newNames = vector(mode = "character", length = nrow(data))
for (j in 1:nrow(data)) {
  cur = data[j,1]
  index = which(genes == cur)
  snpNum[index] = snpNum[index] + 1
  newNames[j] = paste0("gene",index,"_snp",snpNum[index])
}
data$Labels = newNames

rename = cbind(data$genes, data$Labels,data$start, data$end)

write.table(rename, labels, sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)


allData = data

for (i in 2:rep) {
  infile = paste0("rep_", i, "_", i, "_", A, "x", B, "_snp_report.bed") 
  data = read.table(infile, sep = "\t")
  data = subset(data, select = -c(V4))
  names(data) = c("genes", "start", "end", paste0('AxB_',i,'_A'), paste0('AxB_',i,'_B'))
  
  allData = merge(allData, data, by = c("start", "end", "genes"), all.x = TRUE, all.y = TRUE)
}

for (i in 1:rep) {
  infile = paste0("rep_", i, "_", i, "_", B, "x", A, "_snp_report.bed") 
  data = read.table(infile, sep = "\t")
  data = subset(data, select = -c(V4))
  names(data) = c("genes", "start", "end", paste0('BxA_',i,'_A'), paste0('BxA_',i,'_B'))
  
  allData = merge(allData, data, by = c("start", "end", "genes"), all.x = TRUE, all.y = TRUE)
}


saveData = allData[,!(names(allData) %in% c("start", "end", "genes"))] 
rownames(saveData) = saveData$Labels
saveData = saveData[,!(names(saveData) %in% c("Labels"))] 
write.table(saveData,outcounts, sep = "\t", quote = FALSE)

# dataAlt = data[,5:6]
# rownames(dataAlt) = data$Labels
# write.table(dataAlt, outfile, sep = "\t", col.names = FALSE, quote = FALSE)
