##ZJU-cohort
###AML ADT
seob_list2 <- list()

samples2 <- list.files(path = ".")

for (sample in samples2) {
  
  scrna_data2 <- Read10X(
    data.dir = str_c("./", sample), gene.column = 1)
  
  seob2 <- CreateSeuratObject(
    counts = scrna_data2,
    project = sample,
    min.cells = 3,  
    min.features = 200)
  
  seob2[["percent.mt"]] <- PercentageFeatureSet(seob2, pattern = "^MT-")
  
  seob2[["percent.rb"]] <- PercentageFeatureSet(seob2, pattern = "^RP[SL]")
  
  HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(seob2))
  seob2[["percent.HB"]] <- PercentageFeatureSet(seob2, features = HB.genes)
  
  seob2[['sample']] <- sample
  
  seob_list2[[sample]] = seob2
}


seob <- merge(x = seob_list2[[1]], 
              y = seob_list2[-1], 
              add.cell.ids = names(seob_list2))

VlnPlot(seob, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB", "percent.rb"), 
        group.by  = "sample",
        log = T,
        pt.size = 0,
        raster = F)
seob <- subset(seob, subset = nFeature_RNA >300 & percent.mt < 10)

##remove doublets
library(scDblFinder)
library(BiocParallel)
table(seob$sample)
DefaultAssay(seob) <- "RNA"
seob3 <- as.SingleCellExperiment(seob)
seob3 <- scDblFinder(seob3, samples="sample", BPPARAM=MulticoreParam(30))

seob$df <- seob3$scDblFinder.class
table(seob$df)
seob <- subset(seob, subset = df == 'singlet')
seob$group <- ifelse(seob$sample %in% c("NML", "LZC", "FW"), "Mono", "Prim")


ADT <- c("ADT-CD3", "ADT-CD4", "ADT-CD8", "ADT-Fas",
         "ADT-CD27", "ADT-CD25", "ADT-CD69", "ADT-4-1BB",
         "ADT-CD45RA", "ADT-CD62L", "ADT-CCR7", "ADT-PD-1",
         "ADT-LAG-3", "ADT-Tim-3", "ADT-CD7", "ADT-CD14", "ADT-CD33",
         "ADT-HLA-ABC", "ADT-HLA-DR", "ADT-PD-L1", "ADT-CD70", "ADT-CD11b",
         "ADT-CD34", "ADT-CD117")


RNA <- rownames(seob)[!grepl("^ADT-", rownames(seob))]

seob_ADT <-seob[ADT,] 
matrix_ADT <- seob_ADT@assays$RNA@counts

seob_RNA <- seob[RNA,]
matrix_RNA <- seob_RNA@assays$RNA@counts

seob[["RNA"]] <- CreateAssayObject(counts = matrix_RNA)
seob[["ADT"]] <- CreateAssayObject(counts = matrix_ADT)

seob <- NormalizeData(seob, normalization.method = "LogNormalize")
seob <- ScaleData(seob)


seob <- FindVariableFeatures(seob, selection.method = 'vst', nfeatures = 3000)
seob <- RunPCA(seob)
seob <- harmony::RunHarmony(seob, group.by.vars = "sample")
seob <- seob %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.3) %>% 
  RunUMAP(reduction = "harmony", dims = 1:20)

DimPlot(seob, label = T)
DimPlot(seob, group.by = "sample")

DefaultAssay(seob) <- "ADT"
VariableFeatures(seob) <- rownames(seob[["ADT"]])
seob <- NormalizeData(seob, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = "apca")
seob <- FindMultiModalNeighbors(seob, reduction.list = list("harmony", "apca"),
                                dims.list = list(1:20, 1:20), modality.weight.name = "RNA.weight")
seob <- RunUMAP(seob, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                reduction.key = "wnnUMAP_")
seob <- FindClusters(seob, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = F)
DimPlot(seob, reduction = "wnn.umap")

###T cell recluster
seob1 <- subset(seob, subset = seurat_clusters %in% c("0", "1", "3", "6", "9", "10", "12"))
DefaultAssay(seob1) <- "RNA"
seob1 <- NormalizeData(seob1, normalization.method = "LogNormalize")
seob1 <- ScaleData(seob1)
seob1 <- FindVariableFeatures(seob1, selection.method = 'vst', nfeatures = 3000)
seob1 <- RunPCA(seob1)

ElbowPlot(seob1, ndims = 50)

seob1 <- seob1 %>% 
  FindNeighbors(reduction = "pca", dims = 1:10) %>% 
  FindClusters(resolution = 0.4) %>% 
  RunUMAP(reduction = "pca", dims = 1:10) %>% 
  identity()
DimPlot(seob1, label = T)

##ADT
DefaultAssay(seob1) <- "ADT"
VariableFeatures(seob1) <- rownames(seob1[["ADT"]])
seob1 <- NormalizeData(seob1, normalization.method = "CLR", margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = "apca")


seob1 <- FindMultiModalNeighbors(seob1, reduction.list = list("pca", "apca"),
                                 dims.list = list(1:10, 1:10), modality.weight.name = "RNA.weight")

seob1 <- RunUMAP(seob1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seob1 <- FindClusters(seob1, graph.name = "wsnn", algorithm = 3, resolution = 0.3, verbose = F)

DimPlot(seob1, reduction = "wnn.umap", label = T)



##tumor recluster
seob2 <- subset(seob, subset = seurat_clusters %nin% c("0", "1", "3", "6", "9", "10", "12"))
DefaultAssay(seob2) <- "RNA"
seob2 <- NormalizeData(seob2, normalization.method = "LogNormalize")
seob2 <- ScaleData(seob2)
seob2 <- FindVariableFeatures(seob2, selection.method = 'vst', nfeatures = 3000)
seob2 <- RunPCA(seob2)

ElbowPlot(seob2, ndims = 50)

seob2 <- harmony::RunHarmony(seob2, group.by.vars = "sample")
seob2 <- seob2 %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(resolution = 0.4) %>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  identity()
DimPlot(seob2, label = T)

##ADT
DefaultAssay(seob2) <- "ADT"
VariableFeatures(seob2) <- rownames(seob2[["ADT"]])
seob2 <- NormalizeData(seob2, normalization.method = "CLR", margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = "apca")


seob2 <- FindMultiModalNeighbors(seob2, reduction.list = list("harmony", "apca"),
                                 dims.list = list(1:10, 1:10), modality.weight.name = "RNA.weight")

seob2 <- RunUMAP(seob2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seob2 <- FindClusters(seob2, graph.name = "wsnn", algorithm = 3, resolution = 0.3, verbose = F)

DimPlot(seob2, reduction = "wnn.umap", label = T)

###pancer-myeloid cohort
seob_list2 <- list()

samples2 <- c('ESCA', 'LYM', 'KIDNEY', 'MYE', 'OV', 'PAAD', 'THCA', 'UCEC')

for (sample in samples2) {
  
  scrna_data2 <- read.csv(str_c(sample, ".csv"), row.names = 1)
  scrna_data2 <- as.data.frame(t(scrna_data2))
  meta <- read.csv(str_c(sample, "_metadata.csv"), row.names = 1)
  
  seob2 <- CreateSeuratObject(
    counts = scrna_data2,
    meta.data = meta,
    project = sample,
    min.cells = 3,  
    min.features = 200)
  
  seob2[["percent.mt"]] <- PercentageFeatureSet(seob2, pattern = "^MT-")
  
  seob2[["percent.rb"]] <- PercentageFeatureSet(seob2, pattern = "^RP[SL]")
  
  HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(seob2))
  seob2[["percent.HB"]] <- PercentageFeatureSet(seob2, features = HB.genes)
  
  seob2[['sample']] <- sample
  
  seob_list2[[sample]] = seob2
}

seob2[["percent.mt"]] <- PercentageFeatureSet(seob2, pattern = "^MT-")

VlnPlot(seob2, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB", "percent.rb"), 
        group.by  = "cancer",
        log = T,
        pt.size = 0,
        raster = F)


seob2 <- subset(seob2, subset = nFeature_RNA >300 & percent.mt < 10)

FeatureStatPlot(seob2, stat.by = c("CD14", "LYZ", "FCGR3A"), group.by = "seurat_clusters")
seobMy <- subset(seob2, subset = seurat_clusters %in% c("0", "1", "2", "3", "4", "5", "6", "8", "9", "10"))
seobMy <- NormalizeData(seobMy, normalization.method = "LogNormalize")
seobMy <- ScaleData(seobMy)
seobMy <- FindVariableFeatures(seobMy, selection.method = 'vst', nfeatures = 3000)

seobMy <- RunPCA(seobMy)
seobMy <- seobMy %>% RunHarmony(group.by.vars = "cancer")
ElbowPlot(seobMy, ndims = 50)

seobMy <- seobMy %>% 
  FindNeighbors(reduction = "harmony", dims = 1:16) %>% 
  FindClusters(resolution = 0.5) %>% 
  RunUMAP(reduction = "harmony", dims = 1:16) %>% 
  identity()




####GSE235063
seob_list2 <- list()

samples2 <- c('AML10_DX', 'AML10_REL', 'AML10_REM', 'AML11_DX', 'AML11_REL', 'AML11_REM',
              'AML12_DX', 'AML12_REL', 'AML12_REM', "AML13_DX", "AML13_REL", 'AML13_REM',
              'AML14_DX', 'AML14_REM', 'AML15_DX', 'AML15_REL', 'AML15_REM',
              'AML16_DX', 'AML16_REL', 'AML16_REM', 'AML17_DX', 'AML17_REL',
              'AML18_DX', 'AML18_REL', 'AML19_DX', 'AML19_REL',
              'AML1_DX', 'AML1_REM', 'AML20_DX', 'AML20_REM',
              'AML21_DX', 'AML21_REL', 'AML21_REM', 'AML22_DX', 'AML22_REL', 'AML22_REM',
              'AML23_DX', 'AML23_REL', 'AML23_REM', 'AML24_DX', 'AML24_REL', 'AML24_REM',
              'AML25_DX', 'AML25_REL', 'AML25_REM', 'AML26_DX', 'AML26_REL', 'AML26_REM',
              'AML27_DX', 'AML27_REL', 'AML27_REM', 'AML28_REL', 'AML28_REM',
              'AML2_DX', 'AML2_REL', 'AML2_REM', 'AML3_DX', 'AML3_REM',
              'AML4_DX', 'AML4_REL', 'AML5_DX', 'AML5_REL', 'AML5_REM',
              'AML6_DX', 'AML6_REL', 'AML6_REM', 'AML7_DX', 'AML7_REL', 'AML7_REM',
              'AML8_DX', 'AML8_REL', 'AML8_REM', 'AML9_DX', 'AML9_REL', 'AML9_REM')

for (sample in samples2) {
  
  scrna_data2 <- Read10X(
    data.dir = str_c("GSE235063/", sample), gene.column = 1)
  meta <- read.csv(str_c("GSE235063/",sample, "/meta.csv"))
  rownames(meta) <- meta$Cell_Barcode
  
  seob2 <- CreateSeuratObject(
    counts = scrna_data2,
    meta.data = meta,
    project = sample,
    min.cells = 3,  
    min.features = 200)
  
  seob2[["percent.mt"]] <- PercentageFeatureSet(seob2, pattern = "^MT-")
  
  seob2[["percent.rb"]] <- PercentageFeatureSet(seob2, pattern = "^RP[SL]")
  
  HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(seob2))
  seob2[["percent.HB"]] <- PercentageFeatureSet(seob2, features = HB.genes)
  
  seob2[['sample']] <- sample
  
  seob_list2[[sample]] = seob2
}

seob063 <- merge(x = seob_list2[[1]], 
                 y = seob_list2[-1], 
                 add.cell.ids = names(seob_list2))

seob2 <- subset(seob063, subset = nFeature_RNA >300 & percent.mt < 10)

library(scDblFinder)
library(BiocParallel)
table(seob2$sample)


DefaultAssay(seob2) <- "RNA"
seob3 <- as.SingleCellExperiment(seob2)
seob3 <- scDblFinder(seob3, samples="sample", BPPARAM=MulticoreParam(4))

seob2$df <- seob3$scDblFinder.class
table(seob2$df)
seob2 <- subset(seob2, subset = df == 'singlet')


###GSE185381
seob_list <- list()
samples <- dir(path = ".")

##total RNA
#ADT
samples <- samples[c(2:25,27:48)]
for (sample in samples) {
  
  scrna_data <- Read10X(
    data.dir = str_c("./", sample))
  seob <- CreateSeuratObject(
    counts = scrna_data$`Gene Expression`,
    project = sample)
  
  ADT <- CreateAssayObject(counts = scrna_data$`Antibody Capture`)
  
  seob[['ADT']] <- ADT
  
  metadata <- read.csv(str_c(sample, "/metadata.csv"), row.names = 1, check.names = F)
  rownames(metadata) <- metadata$cell
  rownames(metadata) <- colsplit(rownames(metadata), "[:]", names = c("c1", "c2"))$c2
  rownames(metadata) <- paste(rownames(metadata), sep = "-", "1")
  seob <- seob[,rownames(seob@meta.data) %in% rownames(metadata)] 
  seob$samples <- metadata$samples
  seob$malignant <- metadata$malignant
  seob$CNV <- metadata$CNV_pos
  seob$cluster <- metadata$clusters_res.2
  seob$celltype <- metadata$Cell_type_identity
  seob$broad <- metadata$Broad_cell_identity
  Umap <- as.matrix(metadata[,c('UMAP_1', 'UMAP_2')])
  seob[['umap']] <- CreateDimReducObject(embeddings = Umap)
  
  
  seob[["percent.mt"]] <- PercentageFeatureSet(seob, pattern = "^MT-")
  
  seob[["percent.rb"]] <- PercentageFeatureSet(seob, pattern = "^RP[SL]")
  
  HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(seob))
  seob[["percent.HB"]] <- PercentageFeatureSet(seob, features = HB.genes)
  
  seob[['batch']] <- sample
  
  seob_list[[sample]] = seob
}

seob <- merge(x = seob_list[[1]], 
              y = seob_list[-1], 
              add.cell.ids = names(seob_list))

umap_list <- list()
for (i in samples) {
  test <- seob_list[[i]]@reductions[['umap']]
  rownames(test@cell.embeddings) <- paste(i, sep = "_", rownames(test@cell.embeddings))
  umap_list[[i]] <- test
}

umap <- merge(x = umap_list[[1]],
              y = umap_list[-1])

seob[['umap']] <- umap

CellDimPlot(seob, group.by = "broad")

###GSE193891
seobc <- readRDS("/data/Documents/MZM/BLVRB/MAFB KO mouse/GSE193891_IM_mono.HT5-Control.seuratObject.rds")
seobk <- readRDS("/data/Documents/MZM/BLVRB/MAFB KO mouse/GSE193891_IM_mono.HT7-MAFb-KO.seuratObject.rds")

seobt <- merge(seobk, seobc)
seobt <- NormalizeData(seobt, normalization.method = "LogNormalize")
seobt <- ScaleData(seobt)


seobt <- FindVariableFeatures(seobt, selection.method = 'vst', nfeatures = 3000)
seobt <- RunPCA(seobt)

seobt <- seobt %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.4) %>% 
  RunUMAP(reduction = "pca", dims = 1:20)


