###Figure2A
CellDimPlot(seob2, group.by = c("seurat_clusters", "sample", "celltype"), show_stat = F, reduction = "wnn.map", legend.position = "none", label = T, label_insitu = T)

###Figure2C
FeatureDimPlot(seob2, features = "BLVRB", reduction = "wnn.umap")

##Figure2D
countR <- count2tpm(countR, idType = "SYMBOL")
countR <- as.data.frame(t(countR))
countR <- countR[, "BLVRB", drop=F]
countR$patient <- rownames(countR)
countR$group <- ifelse(colsplit(countR$patient, "[_]", names = c("c1", "c2", "c3"))$c3 == "M", "Monocytic", "Primitive")


shapiro.test(countR$BLVRB)
bartlett.test(countR$BLVRB ~ group, data = countR)
countR$group <- factor(countR$group, levels = c("Primitive", "Monocytic"))


ggplot(data=countR,aes(x=group,y=log2(BLVRB+1),
                       colour = group))+ 
  geom_violin(#color = 'grey',
    alpha = 0.8, 
    scale = 'width',
    #linewidth = 1,
    trim = F)+ 
  geom_boxplot(mapping=aes(x=group,y=log2(BLVRB+1),
                           colour=group,fill=group), 
               alpha = 0.5,
               size=1.5,
               width = 0.3)+ 
  geom_jitter(mapping=aes(x=group,y=log2(BLVRB+1),colour = group), 
              alpha = 0.3,size=3)+
  scale_fill_manual(limits=c("Primitive","Monocytic"), 
                    values =c("#69b7e5","#e76060"))+
  scale_color_manual(limits=c("Primitive","Monocytic"), 
                     values=c("#69b7e5","#e76060"))+ 
  geom_signif(mapping=aes(x=group,y=log2(BLVRB+1)), 
              comparisons = list(c("Primitive","Monocytic")),
              map_signif_level=F, 
              tip_length=c(0,0,0,0,0,0),
              y_position = c(9,10), 
              size=1, 
              textsize = 4, 
              test = "t.test", 
              color = "black")+ 
  theme_bw()+ 
  guides(fill = guide_legend(title = "Group"),  
         color = guide_legend(title = "Group"))+  
  labs(x="Group",y="Log2(TPM+1)") 

###Figure2E
data <- rbind(TCGA, TARGET, Beat)
data$test <- log2(data$BLVRB+1)
data$group <- factor(data$group, levels = c("Primitive", "Monocytic"))
ggboxplot(data, x = "cohort", y = "test",
          color = "group")+
  scale_color_manual(values = c("#66c2a5","#fc8d62"))+
  stat_compare_means(aes(group = group))+
  ylab("log2(TPM+1)")+
  theme()

###Figure2F
data <- rbind(TCGA, TARGET, Beat)
data$patient <- rownames(data)
data1 <- rbind(TCGA_ciber, TARGET_ciber, Beat_ciber)

data2 <- left_join(data, data1, by = "patient")
data2$test <- log2(data$BLVRB+1)
data2$Monocyte <- data2$Monocyte*100

data2$cohort <- factor(data2$cohort, levels = c("Beat", "TCGA", "TARGET"))
ggplot(data2, aes(x = Monocyte, y = test))+
  geom_point(aes(color = cohort)) +
  geom_rug(aes(color =cohort)) +
  geom_smooth(aes(color = cohort), method = lm, 
              se = FALSE, fullrange = TRUE)+
  facet_grid(~cohort)+
  scale_color_manual(values = c("#00AFBB", "#FC4E07", "#E7B800"))+
  ggpubr::stat_cor(aes(color = cohort), label.x = 0.3, label.y = 9.5)+
  ylab("log2(TPM+1)")+
  theme_test()


###Figure2G
ciber <- deconvo_tme(eset = countR, method = "cibersort", arrays = FALSE, perm = 10)
BL <- as.data.frame(t(countR))
BL <- BL[,'BLVRB', drop=F]
BL$ID <- rownames(BL)
ciber <- left_join(ciber, BL, by = "ID")
ciber$group <- ifelse(colsplit(ciber$ID, "[_]", names = c("c1", "c2", "c3"))$c3 == "P", "Primitive", "Monocytic")

library(ggforce)
library(ggpubr)
ggplot(ciber, aes(x = Macrophages_M2_CIBERSORT, y = log2(BLVRB+1), color = group))+
  geom_point()+
  geom_smooth(aes(color = group),
              method = "lm", se=F, 
              formula = y ~ x)+
  stat_cor(data = ciber, method = "pearson")+
  theme_bw()

ggplot(ciber, aes(x = Monocytes_CIBERSORT, y = log2(BLVRB+1), color = group))+
  geom_point()+
  geom_smooth(aes(color = group),
              method = "lm", se=F, 
              formula = y ~ x)+
  stat_cor(data = ciber, method = "pearson")+
  theme_bw()

ggplot(ciber, aes(x = T_cells_CD8_CIBERSORT, y = log2(BLVRB+1), color = group))+
  geom_point()+
  geom_smooth(aes(color = group),
              method = "lm", se=F, 
              formula = y ~ x)+
  stat_cor(data = ciber, method = "pearson")+
  theme_bw()




###Figure2H
countR$group1 <- ifelse(countR$BLVRB < median(countR$BLVRB), "Low", "High")
colnames(countR)[2] <- "patient"
TIDE <- read.csv("TIDE/TIDE_results.csv")
TIDE <- TIDE %>% mutate(patient = sapply(strsplit(as.character(Patient), "_"), function(x) paste(x[1:3], collapse = "_")))
countR <- left_join(countR, TIDE, by = "patient")
countR$TCD <- ifelse(countR$Dysfunction <= median(countR$Dysfunction), "Low", "High")
load("/data/Documents/MZM/WJY project/0811counts/cibersort/ciber_BLVRB.rda")
ciber$groupMo <- ifelse(ciber$Monocytes_CIBERSORT < median(ciber$Monocytes_CIBERSORT), "Low", "High")
ciber$groupM2 <- ifelse(ciber$Macrophages_M2_CIBERSORT < median(ciber$Macrophages_M2_CIBERSORT), "Low", "High")
ciber$groupT <- ifelse(ciber$T_cells_CD8_CIBERSORT < median(ciber$T_cells_CD8_CIBERSORT), "Low", "High")
ciber <- ciber[,c(1,29:31)]
colnames(ciber)[1] <- "patient"
countR <- left_join(countR, ciber, by = "patient")
countR1 <- countR[,c(2,3,4,20,21,22,23)]

countR3 <- countR1 %>% make_long(patient, group, group1, groupMo,TCD)
countR3$node <- dplyr::recode(countR3$node,
                              "ZYL_Tumor_P" = "Patient G",
                              "ZLG_Tumor_P" = "Patient I",
                              "YJD_Tumor_P" = "Patient K",
                              "WK_Tumor_P" = "Patient L",
                              "SSH_Tumor_P" = "Patient J",
                              "NML_Tumor_M" = "Patient B",
                              "LZC_Tumor_M" = "Patient A",
                              "LH_Tumor_P" = "Patient M",
                              "JJH_Tumor_P" = "Patient H",
                              "GXW_Tumor_M" = "Patient E",
                              "GLJ_Tumor_M" = "Patient D",
                              "FW_Tumor_M" = "Patient C",
                              "DT_Tumor_P" = "Patient F",
                              "Primitive" = "Primitive",
                              "Monocytic" = "Monocytic",
                              "Low" = "Low",
                              "High" = "High") 

col<- rev(c4a("rainbow",21))
ggplot(countR3, aes(x = x, 
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node),
               label = node)) +  
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size=3,color=1,fill="white")+  
  scale_fill_manual(values = col) +  
  theme_sankey(base_size = 16)


###Figure2I
countR <- as.data.frame(t(countR))
countR$group <- ifelse(countR$BLVRB > median(countR$BLVRB), "high", "low")
rownames(countR) <- paste(rownames(countR), sep = "_", countR$group)
countR <- countR[,-35047]
countR <- as.data.frame(t(countR))
expr <- t(apply(countR, 1, function(x){x-(mean(x))}))
write.table(expr, file = "TIDE/TIDE.txt", sep = "\t", quote = F)
result <- read.csv("TIDE/TIDE_results.csv")
result$group <- reshape2::colsplit(result$Patient, "[_]", c("c1", "c2", "c3", "c4"))$c4

ZJU <- result
ZJU$cohort <- "ZJU_cohort"

list <- list(ZJU, tcga, target, beat)
final <- do.call(rbind, list)

ggplot(final,aes(x=cohort,y=Dysfunction,fill=group))+  
  geom_violin(aes(fill=group),color="white",cex=1,alpha=0.5, trim = F)+
  geom_boxplot(width=0.1,position=position_dodge(0.9),color="white")+ 
  theme_classic()+    
  scale_fill_manual(values = c("#f08a61", "#65ba9d"))+  
  labs(x=NULL,title =NULL)+ 
  theme(axis.title=element_text(size=12),axis.text.y=element_text(size=10,color='black'),
        axis.text.x=element_text(size=12,color='black'))+
  stat_compare_means(label="p.signif",label.y =.5)+
  scale_y_continuous(limits = c(-0.6, 0.6))+
  ylab("T cell Dysfunction score from TIDE")


###Figure2J-K
seob2 <- subset(seob2, subset = celltype %nin% c("DCs", "B", "MEP"))
seob2 <- RunSlingshot(seob2, group.by = "celltype", reduction = "wnn.umap", start = "HSC-like", align_start = T)
seob2 <- RunDynamicFeatures(srt = seob2, lineages = "Lineage1", n_candidates = 200)

results <- cytotrace2(seob2, species = "human",
                      is_seurat = T, slot_type = "counts",
                      full_model = FALSE, ncores = 10)


seob2$cytotrace <- results$CytoTRACE2_Score

IL10 <- "GOBP_INTERLEUKIN_10_PRODUCTION.v2024.1.Hs.gmt"
IL4 <- "GOBP_INTERLEUKIN_4_PRODUCTION.v2024.1.Hs.gmt"
IL13 <- "GOBP_INTERLEUKIN_13_PRODUCTION.v2024.1.Hs.gmt"
Ros <- "GOBP_REACTIVE_OXYGEN_SPECIES_BIOSYNTHETIC_PROCESS.v2025.1.Hs.gmt"

seob2 <- scgmt(seob2, signatures = IL10, method = "AUCell")
seob2 <- scgmt(seob2, signatures = IL4, method = "AUCell")
seob2 <- scgmt(seob2, signatures = IL13, method = "AUCell")
seob2 <- scgmt(seob2, signatures = Ros, method = "AUCell")

ht <- DynamicHeatmap(
  srt = seob2, lineages = "Lineage1",
  pseudotime_palette = "Spectral",
  heatmap_palette = "viridis", cell_annotation = "celltype",
  separate_annotation = c("cytotrace", "IL10_production", "IL4_production", "IL13_production", "ROS_biosynthetic"),
  separate_annotation_palcolor = c4a("superfishel_stone", 6),
  features = c("BLVRB"),
  show_row_names = T
)
print(ht$plot)

CellDimPlot(seob2, group.by = "celltype", reduction = "wnn.umap", lineages = "Lineage1", lineages_span = 0.1)

###Figure2L-M
PEI <- readRDS("/data/Documents/MZM/AML/PEI/PEI_AML.Rds")
PEI <- subset(PEI, subset = clusters %in% c("HSC-like", "CMP-like", "CFUmono-like",
                                            "CD14 Mono", "CD16 Mono", "DC precursors", "DCs"))
PEI <- RunMonocle3(PEI, assay = "RNA",clusters = "clusters", close_loop = F)

PEI <- scgmt(PEI, signatures = IL10, method = "AUCell")
PEI <- scgmt(PEI, signatures = IL4, method = "AUCell")
PEI <- scgmt(PEI, signatures = IL13, method = "AUCell")
PEI <- scgmt(PEI, signatures = Ros, method = "AUCell")

ht1 <- DynamicHeatmap(PEI, lineages = "Monocle3_Pseudotime", pseudotime_palette = "Spectral", 
                      cell_annotation = "clusters", separate_annotation =  c("cytotrace",  "IL10_production", "IL4_production", "IL13_production", "ROS_biosynthetic"), 
                      features = "BLVRB", show_row_names = T, separate_annotation_palcolor = c4a("superfishel_stone", 5), 
                      heatmap_palette = "viridis")

print(ht1$plot)


###Figure2N
FeatureDimPlot(seob2, features = c("CD163","CD36","LILRB4", "S100A8"), reduction = "wnn.umap")


###Figure2O
a <- FetchData(seob2, vars = c("BLVRB", "CD163", "CD36", "LILRB4", "S100A8", "S100A9", "HLA-E"))
colnames <- colnames(a)
test1 <- data.frame(colnames)
for (i in 1:length(colnames)){  
  y=as.numeric(a[,"BLVRB"])
  test <- cor.test(as.numeric(a[,i]),y,type="spearman")    
  test1[i,2] <- test$estimate    
  test1[i,3] <- test$p.value  
}

colnames(test1)[2] <- "correlation"
colnames(test1)[3] <- "pvalue"
test1 <- test1[-1,]
test1 <- as.data.frame(t(test1))
test1$var <- rownames(test1)
test1 <- test1[-1,]
test1 <- test1[,c(7,1,2,3,4,5,6)]
colnames(test1) <- c("var", "CD163", "CD36", "LILRB4", "S100A8", "S100A9", "HLA-E")
str(test1)
library(ggradar)
ggradar(
  test1[1,],
  values.radar=c("-0.7","0","0.7"),
  grid.min=-0.5,
  grid.mid=0,
  grid.max=0.5,
  group.line.width=1,
  group.point.size=4,
  #gridline.min.linetype=1,
  #gridline.mid.linetype=1,
  #gridline.max.linetype=1,
  #group.colours=lcols,
  background.circle.colour="white",
  gridline.mid.colour="grey60",
  #legend.title="KSMPainter",
  legend.position="bottom")
