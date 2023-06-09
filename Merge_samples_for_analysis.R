library(harmony)
library(Seurat)
library(ggplot2)
library(data.table)
library(BiocManager)
library(magrittr)
library(dplyr)
library(stringr)
getwd()
setwd("E:/X/")
HDsxsBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/1-HC_SXS_BM/")
HDsxsBM<- CreateSeuratObject(counts = HDsxsBM, project = "HDbm", min.cells = 3, min.features = 200)
HDsxsBM
HDsxsPB<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/1-HC_SXS_PB/")
HDsxsPB<- CreateSeuratObject(counts = HDsxsPB, project = "HDpb", min.cells = 3, min.features = 200)
HDsxsPB
HDycBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/1-HC_YC_BM/")
HDycBM<- CreateSeuratObject(counts = HDycBM, project = "HDbm", min.cells = 3, min.features = 200)
HDycBM
HDycPB<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/1-HC_YC_PB/YC-PB-CDNA_matrix/")
HDycPB<- CreateSeuratObject(counts = HDycPB, project = "HDpb", min.cells = 3, min.features = 200)
HDycPB
HDlxBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/1-HD_LX_BM_1205/HD_LX_BM_1205_matrix_10X/")
HDlxBM<- CreateSeuratObject(counts = HDlxBM, project = "HDbm", min.cells = 3, min.features = 200)
HDlxBM
HDlxPB<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/1-HD_LX_PB_1205/")
HDlxPB<- CreateSeuratObject(counts = HDlxPB, project = "HDpb", min.cells = 3, min.features = 200)
HDlxPB
HDwllBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/1-HD_WLL_BM_1208/")
HDwllBM<- CreateSeuratObject(counts = HDwllBM, project = "HDbm", min.cells = 3, min.features = 200)
HDwllBM
HDwllPB<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/1-HD_WLL_PB_1208/")
HDwllPB<- CreateSeuratObject(counts = HDwllPB, project = "HDpb", min.cells = 3, min.features = 200)
HDwllPB
PNHphyBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/2-PNH_PHY_BM/PHY-BM_matrix/")
PNHphyBM<- CreateSeuratObject(counts = PNHphyBM, project = "PNHbm", min.cells = 3, min.features = 200)
PNHphyBM
PNHphyPB<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/2-PNH_PHY_PB/")
PNHphyPB<- CreateSeuratObject(counts = PNHphyPB, project = "PNHpb", min.cells = 3, min.features = 200)
PNHphyPB
PNHlzfBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/2-PNH_LZF_BM/")
PNHlzfBM<- CreateSeuratObject(counts = PNHlzfBM, project = "PNHbm", min.cells = 3, min.features = 200)
PNHlzfBM
PNHlhjPB<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/2-PNH_LHJ_PB/")
PNHlhjPB<- CreateSeuratObject(counts = PNHlhjPB, project = "PNHpb", min.cells = 3, min.features = 200)
PNHlhjPB
AAbyhBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/3-AA_BYH_BM/")
AAbyhBM<- CreateSeuratObject(counts = AAbyhBM, project = "AAbm", min.cells = 3, min.features = 200)
AAbyhBM
AAjyBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/3-AA_JY_BM/")
AAjyBM<- CreateSeuratObject(counts = AAjyBM, project = "AAbm", min.cells = 3, min.features = 200)
AAjyBM
AAlsBM<- read.table('E:/20220901_For_merged_CD34_subset/data/3-AA_LS_BM/AAbmLS_matrix.tsv/AAbmLS_matrix.tsv',
                     header = T , row.names=1, sep = '')
AAlsBM <- CreateSeuratObject(counts = AAlsBM,
                              min.cells = 3,
                              min.features = 200, 
                              project = "AAbm")
AAlsBM
AApfmPB<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/3-AA_PFM_PB/PFM-PB-WK_matrix/")
AApfmPB<- CreateSeuratObject(counts = AApfmPB, project = "AApb", min.cells = 3, min.features = 200)
AApfmPB
AAzbsBM<- read.table('E:/20220901_For_merged_CD34_subset/data/3-AA_ZBS_BM/AAbmZBS_matrix.tsv',
                    header = T , row.names=1, sep = '')
AAzbsBM <- CreateSeuratObject(counts = AAzbsBM,
                             min.cells = 3,
                             min.features = 200, 
                             project = "AAbm")
AAzbsBM
AAzfyPB<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/3-AA_ZFY_PB/")
AAzfyPB<- CreateSeuratObject(counts = AAzfyPB, project = "AApb", min.cells = 3, min.features = 200)
AAzfyPB
AAzfyBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/3-AA_ZFY_BM/")
AAzfyBM<- CreateSeuratObject(counts = AAzfyBM, project = "AAbm", min.cells = 3, min.features = 200)
AAzfyBM
AAzwjBM<- Read10X(data.dir = "E:/20220901_For_merged_CD34_subset/data/3-AA_ZWJ_BM/")
AAzwjBM<- CreateSeuratObject(counts = AAzwjBM, project = "AAbm", min.cells = 3, min.features = 200)
AAzwjBM
AAwydBM<- read.table('D:/Experiment/AA_scRNAseq/merge/merge/3-AA_WYD_BM/AAbmWYD1_matrix.tsv',
                     header = T , row.names=1, sep = '')
AAwydBM <- CreateSeuratObject(counts = AAwydBM,
                              min.cells = 3,
                              min.features = 200, 
                              project = "AAwydBM")
AAwydBM
AAwydPB<- read.table('D:/Experiment/AA_scRNAseq/merge/merge/3-AA_WYD_PB/AApbWYD2_matrix.tsv',
                     header = T , row.names=1, sep = '')
AAwydPB <- CreateSeuratObject(counts = AAwydPB,
                              min.cells = 3,
                              min.features = 200, 
                              project = "AAwydPB")
AAwylBM<- read.table('D:/Experiment/AA_scRNAseq/merge/merge/3-AA_WYL_BM/AAbmWYL1_matrix.tsv',
                     header = T , row.names=1, sep = '')
AAwylBM <- CreateSeuratObject(counts = AAwylBM,
                              min.cells = 3,
                              min.features = 200, 
                              project = "AAwylBM")
AAwylBM
AAwylPB<- read.table('D:/Experiment/AA_scRNAseq/merge/merge/3-AA_WYL_PB/AApbWYL2_matrix.tsv',
                     header = T , row.names=1, sep = '')
AAwylPB <- CreateSeuratObject(counts = AAwylPB,
                              min.cells = 3,
                              min.features = 200, 
                              project = "AAwylPB")
AAwyd0620BM <- Read10X(data.dir = "D:/Experiment/AA_scRNAseq/merge/merge/3-AA_WYD_BM_0620/AA_WYD-BM-0620_matrix_10X/")
AAwyd0620BM<- CreateSeuratObject(counts = AAwyd0620BM, project = "AAwyd0620BM", min.cells = 3, min.features = 200)
AAwyd0620BM
AAwyd0620PB <- Read10X(data.dir = "D:/Experiment/AA_scRNAseq/merge/merge/3-AA_WYD_PB_0620/AA_WYD-PB-0620_auto_matrix_10X/")
AAwyd0620PB<- CreateSeuratObject(counts = AAwyd0620PB, project = "AAwyd0620PB", min.cells = 3, min.features = 200)
AAwyd0620PB
AAwyd0827BM <- Read10X(data.dir = "D:/Experiment/AA_scRNAseq/merge/merge/3-AA_WYD_BM_0827/AA_WYD-BM-0827_auto_matrix_10X/")
AAwyd0827BM<- CreateSeuratObject(counts = AAwyd0827BM, project = "AAwyd0827BM", min.cells = 3, min.features = 200)
AAwyd0827BM
AAwyd0827PB <- Read10X(data.dir = "D:/Experiment/AA_scRNAseq/merge/merge/3-AA_WYD_PB_0828/AA_WYD-PB-0828_matrix_10X/")
AAwyd0827PB<- CreateSeuratObject(counts = AAwyd0827PB, project = "AAwyd0827PB", min.cells = 3, min.features = 200)
AAwyd0827PB
AAwylCD3 <- Read10X(data.dir = "D:/Experiment/AA_scRNAseq/merge/merge/3-AA_WYL_BM34_0526/AA_WYL_BM34_0526_auto_matrix_10X/")
AAwylCD3<- CreateSeuratObject(counts = AAwylCD3, project = "AAwylCD3", min.cells = 3, min.features = 200)
AAwylCD3

AAandHDandPBandBM <- merge(x = AAbyhBM, y = list(AAjyBM,AAlsBM,AApfmPB,AAwyd0620BM,
                                                 AAwyd0620PB,AAwyd0827BM,AAwyd0827PB,AAwydBM,AAwydPB,
                                                 AAwylBM,AAwylPB,AAwylCD3,AAzbsBM,AAzfyBM,AAzfyPB,AAzwjBM,HDlxBM,
                                                 HDlxPB,HDsxsBM,HDsxsPB,HDwllBM,HDwllPB,HDycBM,HDycPB,
                                                 PNHlhjPB,PNHlzfBM,PNHphyBM,PNHphyPB),
                           add.cell.ids = c("AAbyhBM","AAjyBM","AAlsBM",'AApfmPB',"AAwyd0620BM",
                                            "AAwyd0620PB","AAwyd0827BM","AAwyd0827PB","AAWYDbm","AAWYDpb",
                                            "AAWYLbm","AAWYLpb","AAwylCD3","AAzbsBM","AAzfyBM","AAzfyPB","AAzwjBM","HDlxBM",
                                            "HDlxPB","HDsxsBM","HDsxsPB","HDwllBM","HDwllPB","HDycBM","HDycPB",
                                            "PNHlhjPB","PNHlzfBM","PNHphyBM","PNHphyPB"
                           ))
AAandHDandPBandBM[["percent.mt"]] <- PercentageFeatureSet(AAandHDandPBandBM, pattern = "^MT-")##############
VlnPlot(AAandHDandPBandBM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)####
plot1 <- FeatureScatter(AAandHDandPBandBM, feature1 = "nCount_RNA", feature2 = "percent.mt")######
plot2 <- FeatureScatter(AAandHDandPBandBM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")####
plot3 <- FeatureScatter(AAandHDandPBandBM, feature1 = "nFeature_RNA", feature2 = "percent.mt")####
CombinePlots(plots = list(plot1, plot3, plot2))
library(stringr)
AAandHDandPBandBM <- subset(AAandHDandPBandBM, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 25)
AAandHDandPBandBM <- NormalizeData(AAandHDandPBandBM) %>% FindVariableFeatures() %>% ScaleData()%>%  RunPCA() %>% FindNeighbors(dims = 1:50) 
library(harmony)
library(Seurat)
library(ggplot2)
library(data.table)
library(BiocManager)
library(magrittr)
library(dplyr)
library(stringr)
system.time({AAandHDandPBandBM <- RunHarmony(AAandHDandPBandBM, group.by.vars = "orig.ident")})
AAandHDandPBandBM <- RunUMAP(AAandHDandPBandBM, reduction = "harmony", dims = 1:30)
AAandHDandPBandBM <- FindNeighbors(AAandHDandPBandBM, reduction = "harmony", dims = 1:50) %>% FindClusters(resolution=1.25)
DimPlot(AAandHDandPBandBM, reduction = "umap", label=T,pt.size=1)
DimPlot(AAandHDandPBandBM, reduction = "umap", label=T,split.by='orig.ident',pt.size=1)
DotPlot(AAandHDandPBandBM, features = c(
  "CD14","FCGR3A",
  "CD3D", "CD8A", "IL7R", 
  "MS4A1","CD19","IGHG","XBP1","MZB1",
  "GATA1","HBD","HBA1","HBA2",
  "CD34","MKI67","HOXA9",
  "CEACAM8","CD15","CD177","S100A8","S100A12",
  "CD1C","CLEC10A","HLA-DQA1","HLA-DPB1","HLA-DRA",
  "JCHAIN","IL3RA","IRF8","IRF7",
  "PPBP","PF4","GATA2",
  "APOE","VCAM1","DCN","COL1A2","COL14A1"
), cols = c("red"), dot.scale = 10) + RotatedAxis() +
  scale_colour_gradient(low = "lightblue", high = "darkred") +
  guides(color = guide_colorbar(title = 'Average Expression'))+ scale_size(range = c(0, 10))