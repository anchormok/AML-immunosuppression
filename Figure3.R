###Figure3K
library(pheatmap)
library(ComplexHeatmap)

##tumor
counts <- read.table("tpm.txt", header = T, row.names = 1, sep = "\t", check.names = F)
markers <- c("BLVRB","CD163", "CD36", "MRC1", "CCL5", "CCL2", "CSF1R", "S100A8", "S100A9",
             "IL10")
markers <- c("IFNG", "CSF2", "TCGF", "IL4", "IL6", "IL8", "CXCL8", "IL10", "TNF", "TGFB1", "TGFB2",
             "IL1B", "IL12A", "IL17A", "IL5", "IL23A", "CCL3", "CCL4", "IL2RA")
counts1 <- counts[markers,]
counts1 <- as.matrix(counts1)

counts1 <- t(scale(t(counts1)))

library(RColorBrewer)
col <- colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(100)
pheatmap(counts1, cellwidth = 12, cellheight = 12,
         cluster_cols = F,cluster_rows = F, color = rev(col))

##CAR-T
df_adjusted <- count2tpm(counts, idType = "SYMBOL")
df_adjusted <- read.table("tpm.txt", header = T, check.names = F)
df_adjusted <- remove_duplicate_genes(df_adjusted, column_of_symbol = "GeneName", method = "mean")

final <- c("PRF1", "TNFRSF9","CXCL10", "ICOSLG", "IFNGR2", "IL1A", "IL1B", "NFKBIA", "RELA")
df_adjusted_sub <- df_adjusted[final,]
df_adjusted_sub <- t(scale(t(df_adjusted_sub)))
col <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
pheatmap(df_adjusted_sub, 
         cellwidth = 12, cellheight = 12,
         cluster_cols = F,cluster_rows = F, color = col)


###Figure3L
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(msigdbr)

geneList <- resdata$log2FoldChange
names(geneList) <- resdata$Gene
geneList <- sort(geneList, decreasing = T)
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
hall <- msigdbr(species = "Homo sapiens", category = "H")

gsea_result <- GSEA(geneList = geneList,
                    minGSSize = 1, 
                    maxGSSize = 1000,
                    pvalueCutoff = 1,
                    TERM2GENE =  reactome[, c("gs_name", "gene_symbol")])

library(GseaVis)
gseaNb(object = gsea_result, geneSetID = c("REACTOME_CD163_MEDIATING_AN_ANTI_INFLAMMATORY_RESPONSE", 
                                           "REACTOME_INTERLEUKIN_10_SIGNALING",
                                           "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING"),
curveCol = cols4all::c4a("hcl.set2", 3), addPval = T)

gseaNb(object = gsea_result, geneSetID = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                                           "HALLMARK_INFLAMMATORY_RESPONSE"),
curveCol = cols4all::c4a("hcl.set2", 3), addPval = T)
