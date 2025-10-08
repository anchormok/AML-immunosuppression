####FigureS2A
FeatureDimPlot(seob2， feature = c("ADT-CD117", "ADT-CD14", "ADT-CD11b"), reduction = "wnn.umap", assay = "ADT")

####FigureS2B
genes <- c("FLT3", "SPINK2", "CD38", #"HSC-like"
           "MPO", "CTSG", "ELANE", #"CMP-like"
           "LYZ", "S100A8", "S100A9", "CD14", #"CFUmono-like", "CD14 mono"
           "HLA-DQA1", "HLA-DPA1", #DCs
           "CD79A", "IGKC", #"B"
           "HBD", "ALAS2") #MEP
seob2$celltype <- factor(seob2$celltype, levels = c("HSC-like", "CMP-like", "CFUmono-like",
                                                    "CD14 Mono", "DCs", "B", "MEP"))

Idents(seob2) <- seob2$celltype
DotPlot(seob2, features = genes) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))+
  labs(x=NULL, y=NULL)+guides(size=guide_legend(order = 3))+
  scale_color_gradientn(values = seq(0,1,0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))

####FigureS2C-D
CellDimPlot(seobMy, group.by = c("Celltypes", "cancer"), show_stat = F, label = T)

####FigureS2E
CellStatPlot(seobMy, stat.by = "Celltypes", group.by = "cancer", label = T)

####FigureS2F
mye_marker <- c("LYZ", "CD14", "FCGR3A", 
                "MAFB","AIF1", "CSF1R", "CD163", "MRC1",  
                "S100A8","S100A9", "CD36",
                "CD1C", "CLEC10A")

seobMy$celltype <- factor(seobMy$celltype, levels = c("MDSC like-Mono-C1",
                                                    "MDSC like-Mono-C2",
                                                    "M2 like-Macro",
                                                    "Mono/Macro-C1",
                                                    "Mono/Macro-C2",
                                                    "Mono/Macro-C3",
                                                    "Mono/Macro-C4",
                                                    "CD16+ Mono-C1",
                                                    "CD16+ Mono-C2",
                                                    "CD16+ Mono-C3",
                                                    "CD16+ Mono-C4",
                                                    "DC-C1",
                                                    "DC-C2",
                                                    "Unknown"))

Idents(seobMy) <- seobMy$celltype
DotPlot(seobMy, features = mye_marker) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))+
  labs(x=NULL, y=NULL)+guides(size=guide_legend(order = 3))+
  scale_color_gradientn(values = seq(0,1,0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33')) + coord_flip()

####FigureS2G
CellCorHeatmap(srt_query = seob2, srt_ref = seobMy,
               query_group = "celltype", ref_group = "celltype",
               nlabel = 3, label_by = "row",
               show_row_names = T, show_column_names = T, heatmap_palette = "Spectral")
