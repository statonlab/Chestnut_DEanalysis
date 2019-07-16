setwd("~/Desktop/Jiali/UTK/chestnut/Chestnut_DEanalysis")

library(DESeq2)

# Read count files
mydata <- read.table("data/gene_counts.txt", check.names = F, stringsAsFactors = F, header = TRUE, row.names = 1)
names(mydata)
sample_types = c("American","American","American","Chinese","Chinese","Chinese","Chinese","Chinese","American","American")
condition = c("canker","healthy","healthy","canker","healthy","healthy","healthy","healthy","healthy","healthy")
MydataDesign = data.frame(row.names = colnames(mydata), species = sample_types, conditions = condition)
MydataDesign
# DESeq
dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign, design = ~species+conditions)
dim(dds)
dds = DESeq(dds)
resultsNames(dds)
res <- results(dds)
res
