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
# gene count log transformation
rld <- rlog(dds, blind = FALSE, fitType = "mean")
# canker vs healthy DEGs
resultsNames(dds)
res <- results(dds, contrast = c("conditions","canker","healthy"),lfcThreshold = 1, alpha = 0.05)
summary(res)
# 181 genes upregulated and 19 downregulated in canker

res_Sig <- data.frame(res[which(res$padj<0.05),])
write.csv(res_Sig,"canker_vs_healthy_all.csv")

# heatmap of the top 50 DEGs comparing canker and healthy
rl.heatmap <- assay(rld[rownames(res[order(res$padj), ])[1:50], ])
library(pheatmap)
pheatmap(rl.heatmap, annotation_col=MydataDesign[,c(1,2)], filename = "./heatmap.png", width = 6, height = 8)

# Interaction of species and conditions
MydataDesign$inter <- paste(MydataDesign$species, MydataDesign$conditions, sep="_")
ddsMF <- DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign, design = ~inter)
ddsMF = DESeq(ddsMF)
# Chinese canker vs healthy
res_Cm <- results(ddsMF, contrast = c("inter", "Chinese_canker","Chinese_healthy"), lfcThreshold = 1, alpha = 0.05)
summary(res_Cm)
res_Cm_Sig <- data.frame(res_Cm[which(res_Cm$padj<0.05),])
write.csv(res_Cm_Sig,"Chinese_canker_vs_healthy.csv")
# 86 upregulated and 6 downregulated genes in Chinese canker compared to Chinese healthy

# American canker vs healthy
res_Cd <- results(ddsMF, contrast = c("inter","American_canker", "American_healthy"), lfcThreshold = 1, alpha = 0.05)
summary(res_Cd)
res_Cd_Sig <- data.frame(res_Cd[which(res_Cd$padj<0.05),])
write.csv(res_Cd_Sig,"American_canker_vs_healthy.csv")
# 226 upregulated and 10 downregulated in American canker compared to American healthy

# vienn diagram of overlap DEGs
Cd_up
Cd_down
Cm_up
Cm_down
All_up
All_down