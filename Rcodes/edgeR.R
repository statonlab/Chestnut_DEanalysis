setwd("~/Desktop/Jiali/UTK/chestnut/Chestnut_DEanalysis/")
library(plyr)
library(edgeR)
# Read count files
mydata <- read.table("data/gene_counts.txt", check.names = F, stringsAsFactors = F, header = TRUE, row.names = 1)
names(mydata)

# keep v4.3 genes
Counts <- mydata[rownames(mydata) %in% gff$gene,]
sample_types = c("American","American","American","Chinese","Chinese","Chinese","Chinese","Chinese","American","American")
condition = c("canker","healthy","healthy","canker","healthy","healthy","healthy","healthy","healthy","healthy")

dgList <- DGEList(counts=Counts, genes=rownames(Counts))
countsPerMillion <- cpm(dgList)
summary(countsPerMillion)

# Filtering
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
summary(cpm(dgList))

# normalization
dgList <- calcNormFactors(dgList, method="TMM")
designMat <- model.matrix(~condition)
designMat

dgList <- estimateCommonDisp(dgList, design=designMat)
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef = ncol(fit$design))
normalized_counts <- data.frame(dgList$pseudo.counts)
colnames(normalized_counts) <- paste0(sample_types,"_",condition)
normalized_counts$gene <- rownames(normalized_counts)

summary(decideTests(lrt))
deGenes <- decideTestsDGE(lrt,adjust.method="BH", p=0.05)
deGenes <- rownames(lrt)[as.logical(deGenes)]
DEG_table <- topTags(lrt, n = Inf, p = 0.05)$table

# Find the DEGs in QTLs
gene_table <- read.csv("blight_QTL_genes_v2.csv")
QTL_DEGs <- gene_table[gene_table$gene %in% DEG_table$genes,]
QTL_DEGs <- merge(QTL_DEGs, DEG_table[,c(1,2,5,6)], by.x = "gene", by.y="genes", all.x=T)
gene_functions <- read.csv("chestnut_gene_functions.csv", header = T)
names(gene_functions)[2] <- "gene"
QTL_DEGs_annot <- join(QTL_DEGs, gene_functions[,-3], by="gene", match="first")

QTL_DEGs_count <- merge(QTL_DEGs[,-c(5:7)], normalized_counts, by = "gene", all.x=T)

write.csv(QTL_DEGs_annot, "blight_QTL_DEGs.csv", row.names = F)
write.csv(QTL_DEGs_count, "blight_QTL_DEGs_normalized_counts.csv", row.names = F)

# Find the DEGs in blight QTLs - 2021 March
gene_table <- read.csv("new_blight_QTL_genes2021.csv")
QTL_DEGs <- gene_table[gene_table$gene %in% DEG_table$genes,]
QTL_DEGs <- merge(QTL_DEGs, DEG_table[,c(1,2,5,6)], by.x = "gene", by.y="genes", all.x=T)
gene_functions <- read.csv("chestnut_gene_functions.csv", header = T)
names(gene_functions)[2] <- "gene"
QTL_DEGs_annot <- join(QTL_DEGs, gene_functions[,-3], by="gene", match="first")

QTL_DEGs_count <- merge(QTL_DEGs[,-c(5:7)], normalized_counts, by = "gene", all.x=T)

QTL_DEGs_annot$American_canker <- QTL_DEGs_count$American_canker
QTL_DEGs_annot$Chinese_canker <- QTL_DEGs_count$Chinese_canker
QTL_DEGs_annot$mean_A_healthy <- rowMeans(QTL_DEGs_count[,c(6,7,13,14)])
QTL_DEGs_annot$mean_C_healthy <- rowMeans(QTL_DEGs_count[,c(9:12)])

write.csv(QTL_DEGs_annot, "new_blight_QTL_DEGs2021.csv", row.names = F)

