AAandHDandPBandBM <- RenameIdents(object = AAandHDandPBandBM, 
                             '0' = "Mono",
                             "1"="T/NK","2"="T/NK",
                             "3"="B_lineage","4"="Eryth_lineage",
                             "5"="HSPC","6"="Eryth_lineage",
                             "7"="Mono","8"="Eryth_lineage",
                             "9"="B_lineage",
                             "10"="Granulocyte","11"="DC","12"="Mono","13"="B_lineage",
                             "14"="HSPC","15"="pDC","16"="B_lineage","17"="B_lineage",
                             "18"="Mono","19"="Meg","20"="Mono",
                             "21"="Eryth_lineage","22"="EC/FB","23"="DC",
                             "24"="Eryth_lineage","25"="B_lineage","26"="T/NK"
)
AAandHDandPBandBMpercentage <- table(AAandHDandPBandBM$orig.ident)
AAandHDandPBandBMpercentage
write.csv(AAandHDandPBandBMpercentage,file = "D:/Experiment/AA_scRNAseq/merge/AAandHDandPBandBMpercentage_number1.csv")
AAandHDandPBandBM$CellType <- Idents(AAandHDandPBandBM)
Idents(AAandHDandPBandBM) <- "replicate"
DimPlot(AAandHDandPBandBM, reduction = "umap")
Idents(AAandHDandPBandBM) <- "CellType"
table(Idents(AAandHDandPBandBM))
table(AAandHDandPBandBM$orig.ident)
prop.table(table(Idents(AAandHDandPBandBM)))
x <- table(Idents(AAandHDandPBandBM), AAandHDandPBandBM$orig.ident)
write.csv(x,file = "D:/number.csv")
y <- prop.table(table(Idents(AAandHDandPBandBM), AAandHDandPBandBM$orig.ident), margin = 2)
write.csv(y,file = "D:/prop.csv")