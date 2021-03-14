setwd("~/Desktop/Jiali/UTK/Apricot/")
library(DESeq2)
library(ggplot2)
library(plyr)
library(Rmisc)
library(gridExtra)

# Get expression of genes from table1 Int. J. Mol. Sci. 2019, 20(16), 3999
## read gene table
inositol <- read.csv('~/Desktop/Jiali/UTK/chestnut/Chestnut_DEanalysis/data/inositol genes.csv', header = T)
genes_inositol <- inositol[grep("At",inositol$Gene.ID), c(1,2)]
genes_inositol$Gene.ID <- gsub("t","T",genes_inositol$Gene.ID)

# Read peach info
Peach_genes <- read.csv("Ppersica_298_v2.1.annotation_info.txt", header = T, comment.char = '', sep='\t')
Peach_genes$Gene.ID <- gsub("\\.\\d","",Peach_genes$Best.hit.arabi.name)
head(Peach_genes)

# get corresponding peach genes
Peach_inositol <- join(genes_inositol, Peach_genes[,c(2,14)], by = "Gene.ID", match = "all")
Peach_inositol <- unique(Peach_inositol)
names(chestnut_inositol)[2] <- "Gene.ID"
combined_inositol <- join(chestnut_inositol,Peach_inositol[,2:3], by="Gene.ID",match="all")
write.csv(combined_inositol,"At_Cm_Pm orthologs.csv")

## apricot gene plotting
DrawPlot <- function(x, name) {
  data <- plotCounts(dds, gene=x,intgroup=c("Stage","genotype"), returnData=TRUE)
  datasum <- summarySE(data, measurevar = "count", groupvars = c("Stage","genotype"))
  datasum$genotype <- factor(datasum$genotype, levels = c("A2137","A1956","A660", "A1267"))
  plot <-
    ggplot(datasum, aes(x=Stage, y=count, color=genotype, group=genotype)) +
    geom_point(size=2) + geom_line(size=1) + geom_errorbar(aes(ymin=count-se,ymax=count+se), width=0.1)+
    theme_classic(base_size = 10) + ylab("Normalized expression level") +
    labs(title = paste0(gsub('\\.v2.1','',x),"-",name)) + scale_color_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF"))
  print(plot)
  #  return()
}


Peach_inositol_filtered <- na.omit(Peach_inositol)
p <- list()
for(i in 1:length(Peach_inositol_filtered$locusName)){
  x = paste0(Peach_inositol_filtered$locusName[i],".v2.1")
  name = Peach_inositol_filtered$Name[i]
  p[[i]] <- DrawPlot(x, name)
}
do.call(grid.arrange,p)

## peach gene plotting
DrawPlot <- function(x, name) {
  data <- plotCounts(dds, gene=x,intgroup=c("Timepoint","Genotype"), returnData=TRUE)
  datasum <- summarySE(data, measurevar = "count", groupvars = c("Timepoint","Genotype"))
  datasum$Timepoint <- gsub("preChill","0", datasum$Timepoint)
  datasum$Genotype <- factor(datasum$Genotype, levels = c("A209", "A340", "A318", "A323"))
  plot <-
    ggplot(datasum, aes(x=Timepoint, y=count, color=Genotype, group=Genotype)) +
    geom_point(size=2) + geom_line(size=1) + geom_errorbar(aes(ymin=count-se,ymax=count+se), width=0.1)+#scale_y_continuous(limits = c(0,15000))+
    scale_color_manual(values=c("#FF9999","#FF3333","#99CCFF","#0033FF")) +
    scale_x_discrete(limits=c("0","100","600","preBloom-early","1000","preBloom-late"))+
    labs(title = paste0(gsub('\\.v2.1','',x),"-",name)) + ylab("Normalized expression level") +
    theme_classic(base_size = 10)

  print(plot)
  #  return()
}
Peach_inositol_filtered <- na.omit(Peach_inositol)
p <- list()
for(i in 1:length(Peach_inositol_filtered$locusName)){
  x = paste0(Peach_inositol_filtered$locusName[i],".v2.1")
  name = Peach_inositol_filtered$Name[i]
  p[[i]] <- DrawPlot(x,name)
}
do.call(grid.arrange,p)

## Add coordinate to the genes
setwd("~/Desktop/Jiali/UTK/chestnut/Chestnut_DEanalysis/")
# read the gene list
inositol_genes <- read.csv("At_Cm_Pm orthologs.csv", header = T)
# read gff files
gff <- read.table("../blast_oak/Castanea_mollissima_scaffolds_v3.2_ALLmRNA.gff", header = F)
names(gff) <- c("contig", "pipeline","type","start","end","empty1","strand","empty2","gene_name")
head(gff$gene_name)
gff$gene_name <- gsub("ID=.*Parent=","",gff$gene_name)
names(gff)[9] <- "Chestnut"
names(gff)
names(inositol_genes)[3] <- "Chestnut"
inositol_genes <- join(inositol_genes, gff[,c(1,4,5,7,9)],by="Chestnut",match = "first")
write.csv(inositol_genes, "inositol_genes_coordinate.csv")
