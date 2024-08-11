###### '

library(Seurat)
  #raw <- Read10X(data.dir = data.path.to.cellranger.outs)
  #raw <- CreateSeuratObject(counts = raw, project = "copycat.test", min.cells = 0, min.features = 0)
  table(sce$sub_celltype)
  sce.copykay<-sce[,sce$sub_celltype%in% c("Epithelial","B")]
  exp.rawdata <- as.matrix(sce.copykay@assays$RNA@counts)
  exp.rawdata[1:5,1:5]
  write.table(exp.rawdata, file="data/sce.exp.rawdata.txt", sep="\t", quote = FALSE, row.names = TRUE)
  
  normalcellnames<-colnames(sce[,sce$sub_celltype%in% c("B")])
  library(copykat)
  copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25,
                          KS.cut=0.1, sam.name="sc", 
                          distance="euclidean", 
                          norm.cell.names=normalcellnames,output.seg="FLASE", 
                          n.cores=50)
  ##
  pred.test <- data.frame(copykat.test$prediction)
  pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
  CNA.test <- data.frame(copykat.test$CNAmat)
  pred.test<-copykat.test$prediction %>% data.frame()
  table(pred.test$copykat.pred)
  write.csv(pred.test,file="data/copykay_epi.csv")
  sce.copykay$copykay<-NULL
  sce.copykay<-AddMetaData(sce.copykay,pred.test)
  
  metadata<-sce@meta.dat
  metadata$copykay[metadata$sub_celltype=="Epithelial"]<-sce.copykay[,sce.copykay$sub_celltype=="Epithelial"]$copykat.pred
  metadata$copykay[metadata$sub_celltype!="Epithelial"]<-"No"
  save(metadata,file = "data/totalscRNAdata.Rdata")
  
  table(metadata$copykay)
  
  
  
  ##############可视化
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)

  chr <- as.numeric(CNA.test$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black','grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR,CHR)
  #################################################
  rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1])
  com.preN <- pred.test$copykat.pred
  pred <- rbPal5(2)[as.numeric(factor(com.preN))]
  
  cells <- rbind(pred,pred)
  col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))
  
  heatmap.3(t(CNA.test[,4:ncol(CNA.test)]),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =50, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
  
  legend("topright", paste("pred.",names(table(com.preN)),sep=""), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[2:1], cex=0.6, bty="n")
  ####
  
  #dev.off()
  
  tumor.cells <- pred.test$cell.names[which(pred.test$copykat.pred=="aneuploid")]
  tumor.cells<-gsub("-",".",tumor.cells)
  tumor.mat <- CNA.test[, which(colnames(CNA.test) %in% tumor.cells)]
  hcc <- hclust(parallelDist::parDist(t(tumor.mat),threads =4, method = "euclidean"), method = "ward.D2")
  hc.umap <- cutree(hcc,2)
  
  rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
  subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
  cells <- rbind(subpop,subpop)
  
  heatmap.3(t(tumor.mat),dendrogram="r", distfun = function(x) parallelDist::parDist(x,threads =50, method = "euclidean"), hclustfun = function(x) hclust(x, method="ward.D2"),
            ColSideColors=chr1,RowSideColors=cells,Colv=NA, Rowv=TRUE,
            notecol="black",col=my_palette,breaks=col_breaks, key=TRUE,
            keysize=1, density.info="none", trace="none",
            cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
            symm=F,symkey=F,symbreaks=T,cex=1, cex.main=4, margins=c(10,10))
  
  legend("topright", c("c1","c2"), pch=15,col=RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex=0.9, bty='n')
  
  