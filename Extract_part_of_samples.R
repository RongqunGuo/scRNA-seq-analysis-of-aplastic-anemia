AAsample = c("AApbWYD2","AApbWYL2","AAbmZBS","AAbmLS","AAbmWYD1","AAbmWYL1")
OnlyAA <- subset(x = AAdata2, subset = orig.ident %in% AAsample)
DimPlot(OnlyAA, reduction = "tsne",label = TRUE, pt.size = 1)
DimPlot(OnlyAA, reduction = "tsne",label = TRUE, pt.size = 1,split.by="orig.ident")
OnlyAA
DotPlot(OnlyAA, features = c("TYMS","MKI67",#Proliferating T cells
                             "CD3D","CD3E","CD3G","TRAC","IL7R","CD4","CCR7","TCF7","LEF1",##naive
                             "AQP3","CD69",#central memory CD4+ T cells expressed high levels of CCR7
                             "CCR6", "CXCR6", "CCL5", "PRDM1",#effector memory CD4+ T cells
                             "IL2RA","FOXP3","CTLA4",#TH cells
                             "CD8A","CD8B","GZMK",#effector memory CD8+T cells
                             "GNLY","GZMB","PRF1",#cytotoxic CD8+ lymphocytes
                             "NCAM1","FCGR3A",#NCAM1=CD56/FCGR3A=CD16
                             "KLRF1", "KLRC1","KLRD1",#
                             "NKG7","CST7",
                             "CD19", "MS4A1", "IGHD", "IGHM", "IL4R", "TCL1A",#na??ve B cells,CD20 (MS4A1)
                             "CD27","CD38","IGHG",#memory B cells-CD27///immature B cells (B3) only expressing CD19 and CD20 (MS4A1)
                             "XBP1","MZB1",#plasma cells
                             "VCAN","S100A8","S100A9","CD14", "ITGAL","ITGAM","CSF3R","CSF1R","CX3CR1", "TYROBP","LYZ","S100A12",#CD14Mono
                             "FCN1","FCGR3B","CDKN1C","MS4A7",#CD16Mono
                             "CD1C","HLA-DPB1","HLA-DPA1","HLA-DQA1","ITGAX","CD1E","FCER1A","CLEC10A","FCGR2B",#DC
                             "CLEC9A",#CD1C+ cDC2 (M4)///CLEC9A+ cDC1///pDC (CLEC4C+CD123+)
                             "IL3RA","JCHAIN","IRF7", "TCF4","LILRA4","CLEC4C",#pDC"GZMB"
                             "PF4","PPBP","GP5","ITGA2B","NRGN","TUBB1","SPARC","RGS18","MYL9","GNG11",#Meg
                             "HBA1","HBA2", "HBB", "GYPA","KLF1","GATA1",#Ery
                             "MYB","GATA2","CD34","KIT","HOXA5","HOXA9","HOXA10","C7","FN1","FBN1"
), cols = c("red"), dot.scale = 10) + RotatedAxis() +
  scale_colour_gradient(low = "cyan3", high = "red") +
  guides(color = guide_colorbar(title = 'Average Expression'))+ scale_size(range = c(0, 10))
DotPlot(AAdata2, features = c("IRF1","IRF2","IRF3"
), cols = c("red"), split.by = "orig.ident",dot.scale = 10) + RotatedAxis() +
  scale_colour_gradient(low = "green", high = "red") +
  guides(color = guide_colorbar(title = 'Average Expression'))+ scale_size(range = c(0, 10))