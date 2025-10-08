###Figure4A
library(data.table)
library(pbapply)
library(plyr)
library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)

rasMat <- read.csv("./aucell.csv",row.names = 1, check.names = F)# data.frame

colnames(rasMat) <- sub("(+)", "", colnames(rasMat), fixed = T)

cellinfo <- seob2@meta.data
celltypes <- names(table(cellinfo$celltype))

ctMat <- lapply(celltypes, function(i) {
  as.numeric(cellinfo$celltype == i)
})
ctMat <- do.call(cbind, ctMat)
colnames(ctMat) <- celltypes
rownames(ctMat) <- rownames(cellinfo)

rssMat <- pblapply(colnames(rasMat), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
  })
})
rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(rasMat)
colnames(rssMat) <- colnames(ctMat)

source("utils/plotRegulonRank.R")

PlotRegulonRank(rssMat, "BLVRB Mono")

###Figure4B
FeatureDimPlot(seob2, features = "MAFB", reduction = "wnn.umap")

###Figure4C
regulon <- "MAFB.gmt"

seob2 <- scgmt(seob2, signatures = regulon, method = "AUCell")
FeatureDimPlot(seob2, features = "MAFB_regulon", reduction = "wnn.umap")

###Figure4D
FeatureStatPlot(seob, stat.by = MAFB_regulon, group.by = "celltype", sort=T)

###Figure4F
counts <- counts["MAFB",]
counts <- as.data.frame(t(counts))
counts$group <- c(rep("Scramble",3), rep("shBLVRB-1", 3), rep("shBLVRB-2",3))
counts$group1 <- c(rep("Scramble",3), rep("shBLVRB", 6))

counts$group <- as.factor(counts$group)
result <- t.test(counts$MAFB ~ group, data = counts)
library(ggplot2)
library(ggpubr)

ggplot(data=counts,aes(x=group,y=log2(MAFB+1),
                   colour = group))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8, 
    scale = 'width',
    #linewidth = 1, 
    trim = F)+ 
  geom_boxplot(mapping=aes(x=group,y=log2(MAFB+1),
                           colour=group,fill=group), 
               alpha = 0.5)+ 
  geom_jitter(mapping=aes(x=group,y=log2(MAFB+1),colour = group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Scramble","shBLVRB-1","shBLVRB-2"), 
                    values =c("#A8ACB9","#b86020","#249261"))+
  scale_color_manual(limits=c("Scramble","shBLVRB-1","shBLVRB-2"), 
                     values=c("#A8ACB9","#b86020","#249261"))+ 
  geom_signif(mapping=aes(x=group,y=log2(MAFB+1)), 
              comparisons = list(c("Scramble","shBLVRB-1"), 
                                 c("Scramble","shBLVRB-2")),
              map_signif_level=F, 
              tip_length=c(0,0,0,0,0,0),
              y_position = c(3,3.2), 
              size=1, 
              textsize = 4, 
              test = "t.test", 
              color = "black")+ 
  theme_bw()+ 
  guides(fill = guide_legend(title = "Group"),  
         color = guide_legend(title = "Group"))+  
  labs(x="Group",y="log2(TPM+1)")+ 
  scale_y_continuous(limits = c(4,5), breaks = c(4,5,0.2))



###Figure4G
test <- read.gmt("/data/Documents/MZM/WJY project/single cell ADT/data1/pyscenic/regulons.gmt")
test <- subset(test, term == "MAFB(573g)")
gsea_result <- GSEA(geneList = geneList,
                    minGSSize = 1, 
                    maxGSSize = 1000,
                    pvalueCutoff = 1,
                    TERM2GENE = test)
library(GseaVis)
gseaNb(object = gsea_result, geneSetID = "MAFB(573g)",
       curveCol = cols4all::c4a("hcl.set2", 3), addPval = T)