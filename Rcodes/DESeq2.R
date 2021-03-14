setwd("~/Desktop/Jiali/UTK/chestnut/Chestnut_DEanalysis")

library(DESeq2)
library(Rmisc)
library(ggplot2)

# Read count files
mydata <- read.table("data/gene_counts.txt", check.names = F, stringsAsFactors = F, header = TRUE, row.names = 1)
names(mydata)
sample_types = c("American","American","American","Chinese","Chinese","Chinese","Chinese","Chinese","American","American")
condition = c("canker","healthy","healthy","canker","healthy","healthy","healthy","healthy","healthy","healthy")
MydataDesign = data.frame(row.names = colnames(mydata), species = sample_types, conditions = condition)
MydataDesign
# Read functional annotation file
ips <- read.table("data/Cm_aa_v4.1_pfam.tsv", check.names = F,stringsAsFactors = F, header = F, sep = '\t',fill = TRUE)
KOALA <- read.table("data/chestnut_annotations.txt", check.names = F, stringsAsFactors = F, header = T, sep = '\t')
blast <- read.table("data/chestnut2peach.blast.txt", check.names = F, stringsAsFactors = F, header = F)
peach_annot <- read.csv("../../Apricot/Ppersica_298_v2.1.annotation_info.txt", header = TRUE, comment.char = '', sep='\t')
# remove duplicate genes/isoforms
ips$V1 <- gsub('.t\\d','',ips$V1)
ips_dedup <- ips[-which(duplicated(ips$V1)), c(1,13,14)]
KOALA$gene_callers_id <- paste0("Cm_",KOALA$gene_callers_id)
KOALA$gene_callers_id <- gsub('.t\\d','',KOALA$gene_callers_id)
KOALA_dedup <- KOALA[-which(duplicated(KOALA$gene_callers_id)),c(1,3,4)]
blast$V1 <- gsub(".t\\d","",blast$V1) # remove chestnut isoforms
blast$V2 <- gsub('0\\.\\d*','0',blast$V2) # remove peach isoforms
blast_dedup <- blast[!duplicated(blast[,1:2]),][,c(1,2,11)]
names(blast_dedup) <- c("Cm_ID","Ppe_ID","Evalue")
blast_dedup_annot <- merge(blast_dedup, peach_annot[,c(2,11,12,13)], by.x = "Ppe_ID", by.y = "locusName", all.x=T)
blast_dedup_annot <- blast_dedup_annot[!duplicated(blast_dedup_annot),]

# DESeq
dds = DESeqDataSetFromMatrix(countData = mydata, colData = MydataDesign, design = ~species+conditions)
dim(dds)
dds = DESeq(dds)
normalized_counts <- counts(dds, normalized = T)
colnames(normalized_counts) <- paste0(sample_types,"_",condition)
# gene count log transformation
rld <- rlog(dds, blind = FALSE, fitType = "mean")
# canker vs healthy DEGs
resultsNames(dds)
res <- results(dds, contrast = c("conditions","canker","healthy"),lfcThreshold = 1, alpha = 0.05)
summary(res)
# 181 genes upregulated and 19 downregulated in canker

res_Sig <- data.frame(res[which(res$padj<0.05),])
res_Sig_anno <- merge(res_Sig, ips_dedup, by.x = "row.names", by.y = 'V1', all.x = T)
res_Sig_anno <- merge(res_Sig_anno, KOALA_dedup, by.x ="Row.names", by.y = 'gene_callers_id', all.x = T)
names(res_Sig_anno)[8:11] <- c("pfam","GO","KEGG","KEGG.function")
res_Sig_anno <- merge(res_Sig_anno, blast_dedup_annot, by.x="Row.names",by.y = "Cm_ID", all.x=T)
write.csv(res_Sig_anno,"canker_vs_healthy_all.csv")

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
res_Cm_Sig_anno <- merge(res_Cm_Sig, ips_dedup, by.x = "row.names", by.y = 'V1', all.x = T)
res_Cm_Sig_anno <- merge(res_Cm_Sig_anno, KOALA_dedup, by.x ="Row.names", by.y = 'gene_callers_id', all.x = T)
names(res_Cm_Sig_anno)[8:11] <- c("pfam","GO","KEGG","KEGG.function")
res_Cm_Sig_anno <- merge(res_Cm_Sig_anno, blast_dedup_annot, by.x="Row.names",by.y = "Cm_ID", all.x=T)
write.csv(res_Cm_Sig_anno,"Chinese_canker_vs_healthy.csv")
# 86 upregulated and 6 downregulated genes in Chinese canker compared to Chinese healthy

# American canker vs healthy
res_Cd <- results(ddsMF, contrast = c("inter","American_canker", "American_healthy"), lfcThreshold = 1, alpha = 0.05)
summary(res_Cd)
res_Cd_Sig <- data.frame(res_Cd[which(res_Cd$padj<0.05),])
res_Cd_Sig_anno <- merge(res_Cd_Sig, ips_dedup, by.x = "row.names", by.y = 'V1', all.x = T)
res_Cd_Sig_anno <- merge(res_Cd_Sig_anno, KOALA_dedup, by.x ="Row.names", by.y = 'gene_callers_id', all.x = T)
names(res_Cd_Sig_anno)[8:11] <- c("pfam","GO","KEGG","KEGG.function")
res_Cd_Sig_anno <- merge(res_Cd_Sig_anno, blast_dedup_annot, by.x="Row.names",by.y = "Cm_ID", all.x=T)
write.csv(res_Cd_Sig_anno,"American_canker_vs_healthy.csv")
# 226 upregulated and 10 downregulated in American canker compared to American healthy

# vienn diagram of overlap DEGs
library(VennDiagram)
Cd_up <- rownames(res_Cd_Sig[which(res_Cd_Sig$log2FoldChange>0),])
Cd_down <- rownames(res_Cd_Sig[which(res_Cd_Sig$log2FoldChange<0),])
Cm_up <- rownames(res_Cm_Sig[which(res_Cm_Sig$log2FoldChange>0),])
Cm_down <- rownames(res_Cm_Sig[which(res_Cm_Sig$log2FoldChange<0),])
All_up <- rownames(res_Sig[which(res_Sig$log2FoldChange>0),])
All_down <- rownames(res_Sig[which(res_Sig$log2FoldChange<0),])

# upregulation genes
area1=length(All_up)
area2=length(Cd_up)
area3=length(Cm_up)

#---pairs
n12=length(intersect(All_up,Cd_up))
n13=length(intersect(All_up,Cm_up))
n23=length(intersect(Cd_up,Cm_up))

#---trios
n123=length(Reduce(intersect,list(All_up,Cd_up,Cm_up)))

grid.newpage()
draw.triple.venn(area1, area2, area3, n12, n13, n23, n123,
                 scaled=FALSE, 
                 fill=c("red", "blue","yellow"), 
                 c("ALL", "American","Chinese"))

# downregulation genes
area1=length(All_down)
area2=length(Cd_down)
area3=length(Cm_down)

#---pairs
n12=length(intersect(All_down,Cd_down))
n13=length(intersect(All_down,Cm_down))
n23=length(intersect(Cd_down,Cm_down))

#---trios
n123=length(Reduce(intersect,list(All_down,Cd_down,Cm_down)))

grid.newpage()
draw.triple.venn(area1, area2, area3, n12, n13, n23, n123,
                 scaled=FALSE, 
                 fill=c("red", "blue","yellow"), 
                 c("ALL", "American","Chinese"))

## 20190930 ------ plot gene expression og given genes
# Functions for plotting
x = "Cm_g21596"
drawPlot <- function(x) {
  data <- plotCounts(dds, gene=x,intgroup=c("species","conditions"), returnData=TRUE)
  data[data == 0.5] <- 0
  datasum <- summarySE(data, measurevar = "count", groupvars = c("species","conditions"))
  datasum$conditions <- factor(datasum$conditions, levels = c("healthy", "canker"))
  plot <-
    ggplot(datasum, aes(x=species, y=count, fill=conditions)) +
    geom_bar(stat = "identity",position=position_dodge()) + geom_errorbar(aes(ymin=count-se,ymax=count+se), width=0.1, position=position_dodge(.9))+
    labs(title = x) + ylab("Normalized expression") +
    scale_fill_manual(values=c('gray20','lightgray')) +
    theme_classic(base_size = 8)
  #theme(text = element_text(size=14),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(plot)
  #  return()
}

# Cm_g21345, Cm_g19316 - ppa003014m orthologs
# Cm_g21596
drawPlot("Cm_g21345")
ggsave("Cm_g21345.png")
drawPlot("Cm_g19316")
ggsave("Cm_g19316.png")
drawPlot("Cm_g21596","gene")
ggsave("Cm_g21596.png")

# Get expression of genes from table1 Int. J. Mol. Sci. 2019, 20(16), 3999
## read chestnut blast to arabidopsis file
blast_arab <- read.table("data/chestnut2arab.blast.txt", header = F)
names(blast_arab) <- c("query", "ref", "identity","length", "mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
## read gene table
inositol <- read.csv('data/inositol genes.csv', header = T)
genes_inositol <- inositol[grep("At",inositol$Gene.ID), c(1,2)]
genes_inositol$Gene.ID <- gsub("t","T",genes_inositol$Gene.ID)
names(genes_inositol)[2] <- "ref"
### remove .1, .2....in the arabidopsis IDs, and get the chestnut ID
blast_arab$ref <- gsub("\\.\\d","",blast_arab$ref)
chestnut_inositol <- join(genes_inositol, blast_arab[,c(1,2)], by = "ref", match = "all")
### compute expression data and add into the datafram
chestnut_inositol$query <- gsub("\\.t\\d","",chestnut_inositol$query)
chestnut_inositol <- unique(chestnut_inositol)
chestnut_inositol <- na.omit(chestnut_inositol) # remove rows with NA

### create an empty data frame
inositol_df <- data.frame(matrix(ncol = 7, nrow = 0))
names(inositol_df) <- c("name","arab","chestnut","A_canker","A_healthy","C_canker","C_healthy")
for (i in 1:length(chestnut_inositol$query)) {
  x = chestnut_inositol$query[i]
  data <- plotCounts(dds, gene=x,intgroup=c("species","conditions"), returnData=TRUE)
  datasum <- summarySE(data, measurevar = "count", groupvars = c("species","conditions"))
  expression <- as.data.frame(c(chestnut_inositol[i,],datasum$count))
  names(expression) <- c("name","arab","chestnut","A_canker","A_healthy","C_canker","C_healthy")
  inositol_df <- rbind(inositol_df, expression)
}

# plot gene expression for each gene
library(gridExtra)
p <- list()
for(i in 1:length(chestnut_inositol$query)){
  x = chestnut_inositol$query[i]
  gene_name <- chestnut_inositol$Name[i]
  p[[i]] <- drawPlot(x,gene_name)
}
do.call(grid.arrange,p)


## plot candidate genes from the paper
candidiate_genes <- read.csv("data/Candidate genes.csv", header = T)
candidiate_genes$Gene.s.. <- gsub("\\s+","",candidiate_genes$Gene.s..)
candidiate_genes <- candidiate_genes[grep("Cm",candidiate_genes$Gene.s..),]
## qPcE.1 genes
E1 <- candidiate_genes[grep("1",candidiate_genes$LGE.QTL.),]
## qPcE.2 genes
E2 <- candidiate_genes[grep("2",candidiate_genes$LGE.QTL.),]
## aPcE.3 genes
E3 <- candidiate_genes[grep("3",candidiate_genes$LGE.QTL.),]
## non-QTL genes
others <- candidiate_genes[which(candidiate_genes$LGE.QTL. == ""),]
p <- list()
for(i in 1:length(others$Gene.s..)){
  x = others$Gene.s..[i]
  p[[i]] <- drawPlot(x)
}
do.call(grid.arrange,p)

drawPlot("Cm_g10607")

write.csv(blast_dedup_annot, "chestnut_gene_functions.csv", row.names = F)
