setwd("~/Desktop/Jiali/UTK/chestnut/Chestnut_DEanalysis/")
gff <- read.table("data/Castanea_mollissima_scaffolds_v4.3_HQ_gene.gff", header = F, sep="\t")
names(gff) <- c("contig","method","type","start","stop","empty1","strand","empty2","gene")
gff$gene <- gsub("ID=","", gff$gene)
contig_pos <- read.csv("data/Castanea_mollissima_scaffolds_v4.3_contigPositions.csv", header = T)
contig_pos$cc_contig <- paste0("lcl_",contig_pos$cc_contig)
contig_pos <- contig_pos[order(contig_pos$cc_LG),]

QTL_interval <- read.csv("data/sorted_blight_interval.csv", header = T)
gene_table <- data.frame()
for (i in 1:13) {
  start_contig <- gff[gff$contig %in% QTL_interval$start_contig[i],] 
  start_gene <- start_contig[which(start_contig$start > QTL_interval$start[i]), c(1,9)]
  a <- which(contig_pos$cc_contig == QTL_interval$start_contig[i]) + 1
  b <- which(contig_pos$cc_contig == QTL_interval$end_contig[i]) - 1
  middle_contigs <- contig_pos$cc_contig[a:b]
  middle_genes <- gff[gff$contig %in% middle_contigs, c(1,9)]
  end_contig <- gff[gff$contig %in% QTL_interval$end_contig[i],] 
  end_gene <- end_contig[which(end_contig$stop < QTL_interval$end[i]), c(1,9)]
  
  QTL_genes <- rbind(start_gene, middle_genes, end_gene)
  QTL_genes$LG <- rep(QTL_interval$cc_LG[i], length(QTL_genes$gene))
  gene_table <- rbind(gene_table, QTL_genes)
}

gene_table <- merge(gene_table, contig_pos[,c(1,3)], by.x = "contig", by.y="cc_contig",all.x=T)
gene_table <- gene_table[,c("LG","cc_pos","contig","gene")]

write.csv(gene_table, "new_blight_QTL_genes2021.csv", row.names = F)

# all genes in a SNP contig
SNP_table <- read.csv("data/significant_SNPs2021.csv", header = T)
SNP_genes <- merge(SNP_table, gff[c(1,9)], by.x="genome_contig", by.y="contig", all.x = T)
SNP_genes$gene <- gsub(";","",SNP_genes$gene)
# gene fold change from edgeR
Expression_all <- topTags(lrt, n = Inf, p = 1)$table
names(Expression_all)[1]<-"gene"
SNP_genes <- join(SNP_genes, Expression_all[,c(1,2,6)], by ="gene", match="first")
SNP_genes_annot <- join(SNP_genes, gene_functions[,-3], by="gene", match="first")
SNP_genes_annot_count <- merge(SNP_genes_annot, normalized_counts, by = "gene", all.x=T)

write.csv(SNP_genes_annot_count, "significant_SNPs_gene_expression.csv", row.names = F)
