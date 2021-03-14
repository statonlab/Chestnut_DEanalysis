setwd("~/Desktop/Jiali/UTK/chestnut/Chestnut_DEanalysis/")
gff <- read.table("data/Castanea_mollissima_scaffolds_v4.3_HQ_gene.gff", header = F, sep="\t")
names(gff) <- c("contig","method","type","start","stop","empty1","strand","empty2","gene")
gff$gene <- gsub("ID=","", gff$gene)
contig_pos <- read.csv("data/Castanea_mollissima_scaffolds_v4.3_contigPositions.csv", header = T)
contig_pos$cc_contig <- paste0("lcl_",contig_pos$cc_contig)
contig_pos <- contig_pos[order(contig_pos$cc_LG),]

QTL_interval <- read.csv("blight_QTL_interval_edit.csv", header = T)
gene_table <- data.frame()
for (i in 1:9) {
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

result_table <- merge(gene_table, normalized_counts, by.x = "gene", by.y = "row.names", all.x=T)
write.csv(result_table, "blight_QTL_genes.csv", row.names = F)

all_genes <- normalized_counts[rownames(normalized_counts) %in% gff$gene,]
write.csv(all_genes, "gene_normalized_counts.csv", row.names = F)

# DE genes
DEG_American <- read.csv("American_canker_vs_healthy.csv", header = T, row.names = 1)
DEG_QTL <- result_table[result_table$gene %in% DEG_American$Row.names,]
write.csv(DEG_QTL, "QTL_DEG_expression_canker_vs_healthy_american.csv")

DEG_Chinese <- read.csv("Chinese_canker_vs_healthy.csv", header = T, row.names = 1)
DEG_QTL <- result_table[result_table$gene %in% DEG_Chinese$Row.names,]
write.csv(DEG_QTL, "QTL_DEG_expression_canker_vs_healthy_chinese.csv")
