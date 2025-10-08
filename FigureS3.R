###FigureS3A
CellDimPlot(PEI, group.by = "Clusters", show_stat = F, label = T)

###FigureS3B
FeatureDimPlot(PEI, features = c("ADT-CD34", "ADT-CD14", "ADT-CD3", "ADT-CD4", "ADT-CD56", "ADT-CD19"), show_stat = F, theme_use = "theme_blank", assay = "ADT")

###FIgureS3C
genes <- c("CD38", "CD34", "SPINK2", #"MEP", "HSC-like"
           "MPO", "CTSG", "ELANE", #"CMP-like"
           "LYZ", "S100A8", "S100A9", "CD14", #"CFUmono-like", "CD14 mono"
           "FCGR3A", #"CD16 mono"
           "VPREB3","VPREB1","CD79A", "IGKC", #"proB", "B"
           "JCHAIN", "MZB1", #"Plasma"
           "CD3D", "CD4", "IL7R", "TCF7", "CD8A",
           "GNLY", "NKG7",
           "LMNA","CST3", "FCER1A", "HLA-DQA1", "HLA-DPA1",
           "TCF4", "LILRA4")

Idents(pei) <- pei$clusters
DotPlot(pei, features = genes) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))+
  labs(x=NULL, y=NULL)+guides(size=guide_legend(order = 3))+
  scale_color_gradientn(values = seq(0,1,0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

###FigureS3D
FeatureDimPlot(PEI, features = "BLVRB", show_stat = F, theme_use = "theme_blank")

###FigureS3E
PEI <- RunDEtest(PEI, group_by = "new", fc.threshold = 1, only.pos = F)
PEI <- RunEnrichment(
  srt = PEI, group_by = "new", db = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome"), species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > log2(1) & p_val_adj < 0.05"
)
EnrichmentPlot(
  srt = PEI, group_by = "new", group_use = "BLVRB Mono",
  plot_type = "network"
)

pei <- RunGSEA(
  srt = PEI, group_by = "new", db = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05"
)
GSEAPlot(srt = PEI, group_by = "new", group_use = "BLVRB Mono", id_use = c("GO:0030099", "GO:0045639", "GO:1903131"), 
          line_color = c("#008F91", "#EE7C79", "#FFCD44"))

###FigureS3F
FeatureStatPlot(seob63, stat.by = "BLVRB", group.by = "celltype", sort = T, plot_type = "box")

###FigureS3E
FeatureStatPlot(seob81, stat.by = "BLVRB", group.by = "celltype", sort = T, plot_type = "box")

