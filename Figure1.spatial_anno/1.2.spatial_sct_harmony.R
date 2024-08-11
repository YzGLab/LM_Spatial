#
dir("spatial")
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)

# �ֱ��ȡÿ��10x�����Ľ���ļ���
if(F){ 
  samples=list.files('spatial/')[-1]
  samples
  library(Seurat)
  sceList = lapply(samples,function(pro){ 
    folder=file.path('spatial/',pro,'filtered_feature_bc_matrix') 
    CreateSeuratObject(counts = Read10X(folder), 
                       project = pro )
  })
  sce.big <- merge(sceList[[1]], 
                   y = c(sceList[[2]],sceList[[3]],sceList[[4]]), 
                   add.cell.ids = samples, 
                   project = "spatial_4")
  sce.big
  table(sce.big$orig.ident)
  save(sce.big,file = 'data/sce.big.merge.spatial_4.Rdata')
  
}
load( 'data/sce.big.merge.spatial_4.Rdata')# spatial 
####矫正一点样本的批次来源SCT+harmony 矫正（最推荐）
sce.big <- SCTransform(sce.big, vars.to.regress = "orig.ident", verbose = FALSE)
######################################
DefaultAssay(sce.big)  = 'SCT'
sce.big <- ScaleData(sce.big,verbose = FALSE) %>% RunPCA(pc.gene = sce.big@var.genes,npcs = 30,verbose = FALSE)  %>% RunHarmony(group.by.vars = "orig.ident")
sce<- FindNeighbors(sce.big, reduction = "harmony", dims = 1:30)
sce <- FindClusters(sce, resolution = 0.4)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30)
sce<- RunTSNE(sce, reduction = "harmony", dims = 1:30)
DimPlot(sce,reduction = "umap",label=T,split.by="orig.ident")
DimPlot(sce,reduction = "tsne",label=T,split.by="orig.ident")
DimPlot(sce,reduction = "umap",label=T)
saveRDS(sce, file = "data/Colon_HC_spatial.RDS")

#####
table(sce$seurat_clusters)

spatial_sce<-readRDS("data/Colon_HC_spatial.RDS")
############## 导出 计算compass 数据
library(readr)
exp_CRC<-NULL
exp_CRC<-c()
exp_CRC<-spatial_sce@assays$RNA@counts
exp_CRC<-as.data.frame(exp_CRC)
exp_CRC<-cbind(symbol=rownames(exp_CRC),exp_CRC)
exp_CRC[1:10,1:10]
exp_CRC<-as.data.frame(exp_CRC)
write_tsv(exp_CRC,file="/home/data/gaoyuzhen/Projects/LM_spatial/Compass/spatial_LM_CRC.tsv")
