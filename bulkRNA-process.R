##pubilc cohorts
##TCGA
setwd("Decon/TCGA/")
load("/mnt/mydisk/Documents/MZM/BLVRB/scRNA1008/music1219/TCGA-LAML_mRNA.Rdata")
rowdata <- rowData(data)
coldata <- colData(data)
mRNA_data <- data[rowdata$gene_type == "protein_coding",]
counts <- assay(mRNA_data, "tpm_unstrand")

symbol <- rowData(mRNA_data)$gene_name
count <- cbind(data.frame(symbol), as.data.frame(counts))
table(duplicated(count$symbol))

suppressPackageStartupMessages(library(tidyverse))

final <- count %>% 
  as_tibble() %>% 
  mutate(meanrow = rowMeans(.[,-1]), .before = 2) %>% 
  dplyr::arrange(desc(meanrow)) %>% 
  distinct(symbol, .keep_all = T) %>% 
  dplyr::select(-meanrow) %>% 
  column_to_rownames(var = "symbol") %>% 
  as.data.frame()

TCGA_ciber <- read.table("TCGA/AML ciber.txt", header = T, sep = "\t")
TCGA_clinic <- read.csv("TCGA/TCGA_Clinical_anno.csv", header = T)

TCGA_clinic <- subset(TCGA_clinic, subset = FAB %in% c("M0", "M1", "M2", "M5"))
TCGA_clinic$group <- ifelse(TCGA_clinic$FAB %in% c("M0", "M1", "M2"), "Primitive", "Monocytic")
TCGA_clinic <- TCGA_clinic[,c(1,3,18,61,84)]
rownames(TCGA_clinic) <- TCGA_clinic$patient
TCGA_clinic <- TCGA_clinic[,-1]
TCGA_clinic$cohort <- "TCGA"

colnames(final) <- str_sub(colnames(final), start = 1, end = 12)
final <- as.data.frame(t(final))
final$patient <- rownames(final)
TCGA_clinic$patient <- rownames(TCGA_clinic)
TCGA <- left_join(final, TCGA_clinic, by = "patient")
TCGA <- TCGA[,c("patient", "group", "BLVRB", "FAB")]
TCGA <- dplyr::filter(TCGA, !is.na(group))
TCGA$cohort <- "TCGA"
rownames(TCGA) <- TCGA$patient
TCGA <- TCGA[,-1]


####
colnames(TCGA)[1] <- "patient"
TCGA <- left_join(TCGA, TCGA_clinic, by = "patient")

data <- data[order(data$FAB),]

rownames(data) <- data$patient
data <- data[,2:23]

data <- as.matrix(data)
data <- t(data)


##TARGET
TARGET <- read.csv("/mnt/mydisk/Documents/MZM/AML/Decon/TARGET/TARGET_AML_TPM.csv", check.names = F)
TARGET_ciber <- read.csv("TARGET/TARGET_AML_CIBERSORT.csv", header = T)
TARGET_clinic <- read.csv("TARGET/Clinical_TARGET_AML.csv", header = T)

colnames(TARGET)[1] <- "gene"
rownames(TARGET) <- TARGET$gene
TARGET <- TARGET[,c(-1)]
TARGET <- as.data.frame(t(TARGET))
TARGET <- TARGET[grepl("03A|09A", rownames(TARGET)), ]

rownames(TARGET) <- str_sub(rownames(TARGET), start = 1, end = 16)
TARGET$patient <- rownames(TARGET)

TARGET_clinic <- subset(TARGET_clinic, subset = FAB %in% c("M0", "M1", "M2", "M5"))
TARGET_clinic$group <- ifelse(TARGET_clinic$FAB == "M5", "Monocytic", "Primitive")
colnames(TARGET_clinic)[1] <- "patient"
TARGET_clinic <- TARGET_clinic[,c(1,4,5)]

TARGET <- left_join(TARGET, TARGET_clinic, by = "patient")
TARGET <- TARGET[,c("patient", "BLVRB", "group", "FAB")]
TARGET <- dplyr::filter(TARGET, !is.na(group))
TARGET$cohort <- "TARGET"
rownames(TARGET) <- TARGET$patient
TARGET <- TARGET[,-1]

#test
data1 <- left_join(TARGET, TARGET_clinic, by = "sample_name")
data1 <- subset(data1, subset = gender %in% c("male", "female"))
data1 <- data1[order(data1$FAB),]
rownames(data1) <- data1$sample_name
TARGET_clinic <- data1[,24:26]
data1 <- data1[,2:23]

data1 <- as.matrix(data1)
data1 <- t(data1)

TARGET_clinic$major_type <- ifelse(TARGET_clinic$FAB %in% c("M0", "M1", "M2", "M3"), "Prim", "Mono")


TARGET_clinic$cohort <- "TARGET"



pheatmap(data1, show_colnames = F, annotation_col = TARGET_clinic2, column_split = TARGET_clinic$major_type,
         cluster_cols = F)



pheatmap(data, show_colnames = F, annotation_col = TCGA_clinic, column_split = TCGA_clinic$major_type,
         cluster_cols = F)


###Beat-AML
load("/mnt/mydisk/Documents/MZM/AML/Decon/Beat-LAML_mRNA.Rdata")
rowdata <- rowData(data)
coldata <- colData(data)
mRNA_data <- data[rowdata$gene_type == "protein_coding",]
counts <- assay(mRNA_data, "tpm_unstrand")

symbol <- rowData(mRNA_data)$gene_name
count <- cbind(data.frame(symbol), as.data.frame(counts))
table(duplicated(count$symbol))

suppressPackageStartupMessages(library(tidyverse))

final <- count %>% 
  as_tibble() %>% 
  mutate(meanrow = rowMeans(.[,-1]), .before = 2) %>% 
  dplyr::arrange(desc(meanrow)) %>% 
  distinct(symbol, .keep_all = T) %>% 
  dplyr::select(-meanrow) %>% 
  column_to_rownames(var = "symbol") %>% 
  as.data.frame()

#colnames(final) <- str_sub(colnames(final), 1, 12)
clinical <- read.table("/mnt/mydisk/Documents/MZM/AML/Decon/beataml_wv1to4_clinical.txt", sep = "\t", check.names = F, header = T)
clinical <- clinical %>% 
  dplyr::arrange(desc(dbgap_rnaseq_sample)) %>% 
  dplyr::slice(1:523)
clinical <- clinical[, c(3, 78)]

clinical <- clinical %>% 
  dplyr::arrange(desc(fabBlastMorphology)) %>% 
  dplyr::slice(6:254) 
colnames(clinical)[1] <- "patient"
colnames(clinical)[2] <- "FAB"
clinical$FAB <- str_sub(clinical$FAB, start = 1, end = 2)

final <- as.data.frame(t(final))
final$patient <- rownames(final)
final <- left_join(final, clinical, by = "patient")
final <- final[, c((ncol(final)-1):ncol(final), 1:(ncol(final)-2))]
final <- final[complete.cases(final$FAB), ]
rownames(final) <- final$patient

clinical <- final[,c(1,2)]

final <- subset(final, subset = FAB %in% c("M0", "M1", "M2", "M5"))
final$group <- ifelse(final$FAB == "M5", "Monocytic", "Primitive")
final <- final[,c("BLVRB", "group", "FAB")]
final$cohort <- "Beat"
Beat <- final



#####ZJU-cohort
counts <- read.csv("count.csv", header = T, sep = "\t")
library(DESeq2)
library(IOBR)
colnames(counts)[1] <- "gene"
counts <- remove_duplicate_genes(counts, column_of_symbol = "gene", method = "mean")

countR <- counts[,c(2, 3, 5, 7, 9, 11, 13:19)]

condition <- factor(c("pri", "pri", "pri", "pri", "mono", "pri",
                      "mono", "pri", "mono", "pri", "mono", "mono", "pri"))
coldata <- data.frame(row.names = colnames(countR), condition)

dds <- DESeqDataSetFromMatrix(countData = countR, colData = coldata, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "mono", "pri"))

resOrdered <- res[order(res$padj), ]

resdata <- merge(as.data.frame(resOrdered), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

names(resdata)[1] <- "Gene"





