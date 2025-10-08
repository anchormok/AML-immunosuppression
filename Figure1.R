
load("/mnt/mydisk/Documents/MZM/WJY project/single cell ADT/data1/seobtcell1220.rda")

##Figure1B
CellDimPlot(seob1, group.by = "celltypes", show_stat = F, reduction = "wnn.map", legend.position = "none", label = T, label_insitu = T)

##Figure1C
ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

meta <- seob1@meta.data
summary <- table(meta[,c("celltype", "group")])

roe <- as.data.frame(ROIE(summary))

convert_to_label <- function(x) {
  if (x > 1) {
    return("+++")
  } else {
    return("-")
  }
}

label_matrix <- apply(roe, c(1, 2), convert_to_label)


Heatmap(roe,
        name = "Value",
        col = rev(colorRampPalette(brewer.pal(n = 9, name = "RdYlGn"))(10)),
        row_names_side = "left",
        column_names_side = "bottom",
        show_row_names = TRUE,
        show_column_names = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(label_matrix[i, j], x, y, gp = gpar(fontsize = 10))
        })

###Figure1D and Figure1E
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME")
PD1 <- subset(reactome, gs_description == "PD-1 signaling")
library(scgmt)
PD1 <- "PD1signaling.gmt"
Tex <- "Tex.gmt"
Teff <- "Teff.gmt"

seob1 <- scgmt(seob1, signatures = Tex, method = "singscore")
seob1 <- scgmt(seob1, signatures = Teff, method = "singscore")
seob1 <- scgmt(seob1, signatures = PD1, method = "singscore")

FeatureDimPlot(seob1, features = c("PD.1.signaling", "IKZF2", "BHLHE40"), reduction = "wnn.umap")
FeatureDimPlot(seob1, features = "ADT-PD-1", reduction = "wnn.umap", assay = "ADT")


data <- seob1@meta.data
average_scores <- data %>%
  group_by(celltype) %>%
  summarise(teff = mean(Teff, na.rm = TRUE),
            tex = mean(Tex, na.rm = TRUE))
cols <- c("NK like T cells" = "#a3c9db",
          "MAIT" = "#2271aa",
          "CD8+ exhausted T cells" = "#aed487",
          "CD8+ PD-1+ T cells" = "#399938",
          "CD8+ memory T cells" = "#f5b970",
          "Gamma Delta T cells" = "#ec7a19",
          "TGFbeta-responsive T cells" = "#f49696")

ggplot(average_scores, aes(x = teff, y = tex, color = celltype))+
  geom_point()+
  scale_color_manual(values = cols)+
  geom_text_repel(label = average_scores$celltype3)+
  theme_bw()+
  theme(legend.position="none")+
  xlab("Effector score")+
  ylab("Exhaustion score")


###Figure1F
####TIDE Prim vs Mono
TCGA <- read.csv("/mnt/mydisk/Documents/MZM/AML/TIDE/TCGA_TIDE_results.csv")
TARGET <- read.csv("/mnt/mydisk/Documents/MZM/AML/TIDE/TARGET_TIDE_result.csv")
Beat <- read.csv("/mnt/mydisk/Documents/MZM/AML/TIDE/Beat_TIDE_result.csv")


TCGA_clinic <- read.csv("TCGA/TCGA_Clinical_anno.csv", header = T)
TCGA_clinic <- subset(TCGA_clinic, subset = FAB %in% c("M0", "M1", "M2", "M5"))
TCGA_clinic$group <- ifelse(TCGA_clinic$FAB %in% c("M0", "M1", "M2"), "Primitive", "Monocytic")
TCGA_clinic <- TCGA_clinic[,c(1,3,18,61,84)]
colnames(TCGA_clinic)[1] <- "Patient"
TCGA_clinic$cohort <- "TCGA"

TCGA$Patient <- str_sub(TCGA$Patient, start = 1, end = 12)

TCGA <- left_join(TCGA, TCGA_clinic, by = "Patient")
TCGA <- dplyr::filter(TCGA, !is.na(group))




TARGET_clinic <- read.csv("TARGET/Clinical_TARGET_AML.csv", header = T)


TARGET_clinic <- subset(TARGET_clinic, subset = FAB %in% c("M0", "M1", "M2", "M5"))
TARGET_clinic$group <- ifelse(TARGET_clinic$FAB == "M5", "Monocytic", "Primitive")
colnames(TARGET_clinic)[1] <- "Patient"
TARGET_clinic$cohort <- "TARGET"

TARGET$Patient <- reshape2::colsplit(TARGET$Patient, "[_]", names = c("c1","c2"))$c1
TARGET <- left_join(TARGET, TARGET_clinic, by = "Patient")
TARGET <- dplyr::filter(TARGET, !is.na(group))



clinical <- read.table("/mnt/mydisk/Documents/MZM/AML/Decon/beataml_wv1to4_clinical.txt", sep = "\t", check.names = F, header = T)
clinical <- clinical %>% 
  dplyr::arrange(desc(dbgap_rnaseq_sample)) %>% 
  dplyr::slice(1:523)
clinical <- clinical[, c(3, 78)]

clinical <- clinical %>% 
  dplyr::arrange(desc(fabBlastMorphology)) %>% 
  dplyr::slice(6:254) 
colnames(clinical)[1] <- "Patient"
colnames(clinical)[2] <- "FAB"
clinical$FAB <- str_sub(clinical$FAB, start = 1, end = 2)
clinical <- subset(clinical, subset = FAB %in% c("M0", "M1", "M2", "M5"))
clinical$group <- ifelse(clinical$FAB == "M5", "Monocytic", "Primitive")

Beat$Patient <- reshape2::colsplit(Beat$Patient, "[_]", names = c("c1","c2"))$c1
Beat <- left_join(Beat, clinical, by = "Patient")
Beat <- dplyr::filter(Beat, !is.na(group))
Beat$cohort <- "Beat"

TCGA <- TCGA[,-c(16,17)]
TARGET <- TARGET[,-c(16,17)]


data <- rbind(TCGA, TARGET, Beat)
ggviolin(data, x = "group", y = "Dysfunction", fill = "group",
         palette = rev(c("#f08a61", "#65ba9d")),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = list(c("Primitive", "Monocytic")), label = "p.signif", size = 10,
                     bracket.size = 0.5, tip.length = 0.02, method = "t.test")+
  facet_grid(~cohort)



###Figure1G
TCGA_ciber <- read.table("TCGA/AML ciber.txt", header = T, sep = "\t", check.names = F)
TARGET_ciber <- read.csv("TARGET/TARGET_AML_CIBERSORT.csv", header = T, check.names = F)

Beat_ciber <- deconvo_tme(eset=final,method="cibersort",
                          arrays=TRUE,perm=100)

colnames(TCGA_ciber)[1] <- "patient"
TCGA <- left_join(TCGA, TCGA_ciber, by = "patient")
TCGA <- dplyr::filter(TCGA, !is.na(Monocyte))
TCGA$cohort <- "TCGA"
rownames(TCGA) <- TCGA$patient
TCGA <- TCGA[,-1]


colnames(TARGET_ciber)[1] <- "patient"
TARGET <- left_join(TARGET, TARGET_ciber, by = "patient")
TARGET <- dplyr::filter(TARGET, !is.na(Monocyte))
TARGET$cohort <- "TARGET"
rownames(TARGET) <- TARGET$patient
TARGET <- TARGET[,-1]


colnames(Beat_ciber)[1] <- "patient"
Beat <- cbind(Beat_ciber, Beat)
Beat$cohort <- "BEAT"


data <- rbind(TCGA, TARGET, Beat)
data$group <- factor(data$group, levels = c("Primitive", "Monocytic"))
ggboxplot(data, x = "cohort", y = "CD8+ T cells",
          color = "group")+
  scale_color_manual(values = c("#b0d3e7","#e7afae"))+
  stat_compare_means(aes(group = group))+
  ylab("Proportion")+
  theme()

##Figure1H
###GSEA
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(msigdbr)

geneList <- resdata$log2FoldChange
names(geneList) <- resdata$Gene
geneList <- sort(geneList, decreasing = T)

hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
gsea_result <- GSEA(geneList = geneList,
                    minGSSize = 1, 
                    maxGSSize = 1000,
                    TERM2GENE = reactome[, c("gs_name", "gene_symbol")],
                    pvalueCutoff = 1)
gseaNb(object = gsea_result, geneSetID = c("REACTOME_CD163_MEDIATING_AN_ANTI_INFLAMMATORY_RESPONSE", 
                                           "REACTOME_INTERLEUKIN_10_SIGNALING",
                                           "REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING"),
curveCol = cols4all::c4a("hcl.set2", 3), addPval = T)

