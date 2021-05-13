setwd("~/Desktop/Jiali/UTK/chestnut/Chestnut_DEanalysis/")
library(WGCNA)

## read gene normalized counts
options(stringsAsFactors = F)
normalized_data <- read.csv("gene_normalized_counts.csv", header = T)
rownames(normalized_data) <- rownames(Counts)
names(normalized_data)
## read gene function
gene_function <- read.csv("chestnut_gene_functions.csv", header = T)
gene_function <- gene_function[!duplicated(gene_function$Cm_ID),]

## only keep the genes from Bert
SNP_genes <- read.csv("data/sorted_SNPs_gene_expression_sigtests_cdn.csv", header = T)
filtered_Data <- normalized_data[rownames(normalized_data) %in% SNP_genes$gene,]
# run co-express
datExpr0 = as.data.frame(t(normalized_data))
datExpr0 <- as.data.frame(lapply(datExpr0,as.numeric)) #datExpr0 has to be numeric.
rownames(datExpr0) <- colnames(normalized_data)

# check missing value and outliers
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
# Returns True means all genes pass the cuts, if not, the following scripts can remove the outliers.
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Auto one step construction
allowWGCNAThreads() # Work in Rstudio

# Choose a set of soft-thresholding powers
powers = c(c(1:16))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# One-step network construction and module detection
net = blockwiseModules(datExpr0, power = 8,
                       TOMType = "unsigned", minModuleSize = 1800,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       #saveTOMs = TRUE,
                       #saveTOMFileBase = "shortTOM",
                       verbose = 3)
table(net$colors)

# Load genotype data
datTraits0 = read.csv("data/dataTrait.csv", header = T,row.names = 1,stringsAsFactors = T)
rownames(datTraits0)  = rownames(datExpr0)
rownames(datExpr0)


# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
moduleColors = labels2colors(net$colors)
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEs1 = net$MEs
MEs1 = MEs1[,order(names(MEs1))]
moduleTraitCor = cor(MEs1, datTraits0, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(11,11)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
pdf( "./MEbar.pdf", width = 8, height = 4 )
par(mfrow = c(1,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits0),
               yLabels = names(MEs1),
               ySymbols = names(MEs1),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-factors relationships"))
ME_num <- data.frame(table(net$colors))
ME_num$Var1 <- paste0("ME",ME_num$Var1)
rownames(ME_num) <- ME_num$Var1
ME_num <- ME_num[rev(names(MEs1)),]
barplot(ME_num$Freq, main="Number of genes in each module",
        xlab="Number of genes", las=2,names.arg = ME_num$Var1,
        horiz = T)
dev.off()


## export to cytoscape
TOM = TOMsimilarityFromExpr(datExpr0, power = 8);
moduleLabels = net$colors
# Select modules
for(i in 0:4){
modules = i;
# Select module probes
probes = names(datExpr0)
inModule = is.finite(match(moduleLabels, modules));
modProbes = probes[inModule];
# find interact genes
genes_inModule <- SNP_genes[SNP_genes$gene %in% modProbes,]
ME_genes <- data.frame(Cm_ID = modProbes)
ME_genes_annot <- merge(ME_genes, gene_function, by = "Cm_ID", all.x=T)
write.csv(genes_inModule,paste0("coexpress/SNP genes in"," ME", i,".csv"), row.names = F)
write.csv(ME_genes_annot,paste0("coexpress/All genes in ", "ME",i,".csv"), row.names = F)
}

