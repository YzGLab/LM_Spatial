#
library(harmony)
library(Seurat)
#read 10x and merge them 
if(F){ 
  samples=list.files('scRNA/')
  samples
  library(Seurat)
  sceList = lapply(samples,function(pro){ 
    folder=file.path('scRNA/',pro,'sce.big_feature_bc_matrix') 
    CreateSeuratObject(counts = Read10X(folder), 
                       project = pro )
  })
  sce.big <- merge(sceList[[1]], 
                   y = c(sceList[[2]],sceList[[3]],sceList[[4]]), 
                   add.cell.ids = samples, 
                   project = "ls_4")
  sce.big
  table(sce.big$orig.ident)
  save(sce.big,file = 'sce.big.merge.ls_4.Rdata')
  
}
load('/home/data/gaoyuzhen/Projects/LM_spatial/data/sce.big.merge.ls_4.Rdata')

##RunHarmony
DefaultAssay(sce.big)<-"RNA"
sce.big <- NormalizeData(sce.big)
###Scale可以选择所有基因还是高变基因
sce.big <- FindVariableFeatures(sce.big, selection.method = "vst", nfeatures = 2000) %>% ScaleData(verbose = FALSE) %>% RunPCA(pc.gene = sce.big@var.genes,npcs = 30,verbose = FALSE)   
#sce.big <- ScaleData(sce.big,verbose = FALSE) %>% RunPCA(pc.gene = sce.big@var.genes,npcs = 20,verbose = FALSE)  %>% RunHarmony(group.by.vars = "orig.ident")
sce.big <- RunHarmony(sce.big, group.by.vars = "orig.ident")
############ 可视化
sce<- FindNeighbors(sce.big, reduction = "harmony", dims = 1:30)
sce <- FindClusters(sce, resolution = 0.3)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30)
sce<- RunTSNE(sce, reduction = "harmony", dims = 1:30)
DimPlot(sce,reduction = "umap",label=T,split.by="orig.ident")
DimPlot(sce,reduction = "tsne",label=T,split.by="orig.ident")
saveRDS(sce, file = "data/Colon_HC_combined.sce_sce.big.rds")
#### visulization
sce<-readRDS("data/Colon_HC_combined.sce_sce.big.rds")