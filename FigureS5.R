###FigureS5A
FeatureDimPlot(PEI, features = "MAFB", reduction = "umap")

###FigureS5B
regulon <- "MAFB.gmt"

PEI <- scgmt(PEI, signatures = regulon, method = "AUCell")
FeatureDimPlot(PEI, features = "MAFB_regulon", reduction = "umap")

###FigureS5C
library(clusterProfiler)
library(org.Hs.eg.db)
test <- read.gmt("/data/Documents/MZM/WJY project/single cell ADT/data1/pyscenic/regulons.gmt")
test <- subset(test, term == "MAFB(573g)")
genes <- bitr(test$gene, fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = "org.Hs.eg.db")

go <- enrichGO(gene = genes$ENTREZID,
                 OrgDb = org.Hs.eg.db,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   minGSSize = 10,
                   maxGSSize = 500,
                 readable = T,
                 ont = "BP")
test <- go@result
test <- subset(test, subset = Description %in% c("myeloid cell differentiation",
                                                 "mononuclear cell differentiation",
                                                 "macrophage differentiation",
                                                 "macrophage activation",
                                                 "wound healing",
                                                 "regulation of inflammatory response",
                                                 "interleukin-10 production",
                                                 "interleukin-4 production",
                                                 "interleukin-13 production",
                                                 "macrophage cytokine production"))
test$Description <- factor(test$Description, levels = rev(c("myeloid cell differentiation",
                                                            "mononuclear cell differentiation",
                                                            "macrophage activation",
                                                            "macrophage differentiation",
                                                            "macrophage cytokine production",
                                                            "regulation of inflammatory response",
                                                            "wound healing",
                                                            "interleukin-10 production",
                                                            "interleukin-4 production",
                                                            "interleukin-13 production")))
ggplot(data = test, aes(x = Count, y = Description, fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "YlOrRd", direction = 1)+
  geom_bar(stat = "identity", width = 0.8, alpha = 0.7)+
  labs(x = "Number of Gene", y = "Pathway", title = "GO_BP enrichment")+
  geom_text(aes(x = 0.03, label = Description), hjust = 0)+
  theme_classic()+
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        axis.text.y = element_blank())

###FigureS5H-I
CellDimPlot(seobt, group.by = c("seurat_clusters", "sample"), show_stat = F, reduction = "umap", legend.position = "none", label = T, label_insitu = T)

###FigureS5J
CellStatPlot(seobt, stat.by = "seurat_clusters", group.by = "sample")

###FigureS5K
genes <- c("Cd14", "Egr1", "Mrc1", "Trem2")

library(scgmt)
Mar <- "M1_M2_macrophage-mouse.gmt"
seobt <- scgmt(seobt, signatures = Mar, method = "AUCell")

seobt$seurat_clusters <- factor(seobt$seurat_clusters, levels = c("1", "6", "3", "4", "5", "0", "8", "7", "2"))
GroupHeatmap(
  srt = seobt,
  features = genes,
  group.by = "seurat_clusters",
  heatmap_palette = "YlOrRd",
  cell_annotation = c("M2"),
  cell_annotation_palette = c("Paired"),
  #show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)

###FigureS5L
library(ggradar2)

TCGA1 <- subset(TCGA, subset = group == "Primitive")
TCGA2 <- subset(TCGA, subset = group == "Monocytic")

TARGET1 <- subset(TARGET, subset = group == "Primitive")
TARGET2 <- subset(TARGET, subset = group == "Monocytic")

Beat1 <- subset(Beat, subset = group == "Primitive")
Beat2 <- subset(Beat, subset = group == "Monocytic")

TCGA1 <- TCGA1[,c(1,2,4,3, 5:12)]
TCGA2 <- TCGA2[,c(1,2,4,3, 5:12)]
TARGET1 <- TARGET1[,c(1,2,4,3, 5:12)]
TARGET2 <- TARGET2[,c(1,2,4,3, 5:12)]
Beat1 <- Beat1[,c(1,2,4,3, 5:12)]
Beat2 <- Beat2[,c(1,2,4,3, 5:12)]



cohort <- list(TCGA1, TCGA2, TARGET1, TARGET2, Beat1, Beat2)
names(cohort) <- c("TCGA1", "TCGA2", "TARGET1", "TARGET2", "Beat1", "Beat2")
list <- list()

for ( i in seq_along(cohort)) {
  df <- cohort[[i]]
  correlation_results <- sapply(df[4:11], function(x) cor(df$MAFB, x))
  correlation_df <- data.frame(Column = names(correlation_results), Correlation = correlation_results, group = names(table(df$group)))                              
  list[[i]] <- correlation_df
  
}

names(list) <- names(cohort)

TCGA <- rbind(list$TCGA1, list$TCGA2)
TARGET <- rbind(list$TARGET1, list$TARGET2)
Beat <- rbind(list$Beat1, list$Beat2)
TCGA <- spread(TCGA, key = "Column", value = "Correlation")
TARGET <- spread(TARGET, key = "Column", value = "Correlation")
Beat <- spread(Beat, key = "Column", value = "Correlation")

ggradar2(Beat,polygonfill = FALSE,
         gridline.label = c(-1, -0.5, 0, 0.5, 1),
         group.colours = c("#dc625f", "#64acd5"))

###FigureS5M
counts <- read.table("GSE205719_MFABsh_TPM.tsv", header = T, sep = "\t", check.names = F)

counts <- remove_duplicate_genes(counts, column_of_symbol = "GeneID", method = "mean")
counts <- counts[complete.cases(counts),]

count <- counts[,c(7:12)]
markers <- c("MAFB","BLVRB","CD163", "CD36", "MRC1","CSF1R", "S100A9",
             "IL10")
counts1 <- count[markers,]
counts1 <- as.matrix(counts1)

counts1 <- t(scale(t(counts1)))

library(RColorBrewer)
col <- colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(100)
##提前对矩阵scale与直接在pheatmap中使用scale效果相类似
#pheatmap(counts1, color = col, cluster_rows = F, cluster_cols = F, scale = "row")
pheatmap(counts1, cellwidth = 12, cellheight = 12,
         cluster_cols = F,cluster_rows = F, color = rev(col))

