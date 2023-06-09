suppressor=read.table("suppressor.txt")[,1]pressor
suppressor=as.character(suppressor)
suppressor <- list(suppressor)

OR1<- AddModuleScore(
  object = AAandHDandPBandBM,
  features =suppressor,
  ctrl = 100, #Ĭame = 'suppressor'
)
colnames(OR1@meta.data)
colnames(OR1@meta.data)[8] <- 'suppressor' 
VlnPlot(OR1,features = 'suppressor')
median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}
Vlnlot(object = OR1, features = 'suppressor',split.by = "orig.ident") 
      + stat_summary(fun.y = median.stat, geom='point', size = 5, colour = "blue")