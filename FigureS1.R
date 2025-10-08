###FigureS1C
CellDimPlot(seob, group.by = "seurat_clusters", show_stat = F, reduction = "wnn.map", legend.position = "none", label = T, label_insitu = T)

###FigureS1D
CellDimPlot(seob, group.by = "T cells", show_stat = F, reduction = "wnn.map", legend.position = "none", label = T, label_insitu = T)

###FigureS1E-F
FeatureDimPlot(seob， feature = c("ADT-CD4", "ADT-CD8"), reduction = "wnn.umap", assay = "ADT")

###FigureS1G
marker1 <- c("CD8A","GNLY", "GZMB", "GZMH", "IFNG", "PRF1" ,"KLRG1","GZMK", "IKZF2", "TIGIT",
             "TRDC", "TRGC1","LEF1", "SELL", "TGFBR3", "SMAD3")
Idents(seob11) <- seob11$celltype
DotPlot(seob11, features = marker1) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))+
  labs(x=NULL, y=NULL)+guides(size=guide_legend(order = 3))+
  scale_color_gradientn(values = seq(0,1,0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

###FigureS1H
CellStatPlot(seob1, stat.by = "seurat_clusters", group.by = "group", label = T)

###FigureS1I
meta <- seob1@meta.data
meta$meta.cluster <- meta$celltype
meta$loc <- meta$group
out.prefix <- "./Fig_OR"


OR.immune.list <- do.tissueDist(cellInfo.tb=meta,
                                out.prefix=sprintf("%s.CD8T_cell",out.prefix),
                                pdf.width=4,pdf.height=8,verbose=1
)


a=OR.immune.list[["OR.dist.tb"]]
a <- as.data.frame(a)
rownames(a) <- a$rid
a <- a[,-1]
a <- na.omit(a)

b <- OR.immune.list$count.dist.melt.ext.tb[,c(1,2,6)]
b <- spread(b,key = "cid", value = "adj.p.value")
b <- data.frame(b[,-1],row.names = b$rid)
b <- b[rownames(a),]



col <- viridis(11,option = "D")
b = ifelse(b >= 0.05&(a>1.5|a<0.5), "",
           ifelse(b<0.0001&(a>1.5|a<0.5),"****",
                  ifelse(b<0.001&(a>1.5|a<0.5),"***",
                         ifelse(b<0.01&(a>1.5|a<0.5),"**",
                                ifelse(b < 0.05&(a>1.5|a<0.5),"*","")))))

bk=c(seq(0,0.99,by=0.01),seq(1,2,by=0.01))

pheatmap(a[,], border_color = "NA", fontsize = 9,cellheight = 12,cellwidth = 20,clustering_distance_rows="correlation",
         display_numbers = b,number_color="black",fontsize_number=10,
         cluster_col=F, cluster_rows=T, border= NULL, breaks=bk, 
         color = c(colorRampPalette(colors = col[1:6])(length(bk)/2),
                   colorRampPalette(colors = col[6:11])(length(bk)/2)))

###FigureS1J
library(data.table)
library(pbapply)
library(plyr)
library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)

rasMat <- read.csv("./aucell.csv",row.names = 1, check.names = F)# data.frame

colnames(rasMat) <- sub("(+)", "", colnames(rasMat), fixed = T)

cellinfo <- seob1@meta.data
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

PlotRegulonRank(rssMat, "NK like T cells")

###FigureS1K
regulon <- "BHLHE.gmt"

seob1 <- scgmt(seob1, signatures = regulon, method = "AUCell")

FeatureDimPlot(seob1, reduction = "wnn.umap", feature = "regulon")