#
dir("scRNA")
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)

# �ֱ��ȡÿ��10x�����Ľ���ļ���
if(F){ 
  samples=list.files('scRNA/')
  samples
  library(Seurat)
  sceList = lapply(samples,function(pro){ 
    folder=file.path('scRNA/',pro,'filtered_feature_bc_matrix') 
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
load( '/home/data/gaoyuzhen/Projects/LM_spatial/data/sce.big.merge.ls_4.Rdata')
####矫正一点样本的批次来源SCT+harmony 矫正（最推荐）
#sce.big <- SCTransform(sce.big, vars.to.regress = "orig.ident", verbose = FALSE)
#DefaultAssay(sce.big)  = 'SCT'
#sce.big <- ScaleData(sce.big,verbose = FALSE) %>% RunPCA(pc.gene = sce.big@var.genes,npcs = 20,verbose = FALSE)  %>% RunHarmony(group.by.vars = "orig.ident")

##筛选 线粒体基因和 核糖体基因
raw_sce=sce.big
raw_sce
rownames(raw_sce)[grepl('^MT-',rownames(raw_sce),ignore.case = T)]
rownames(raw_sce)[grepl('^RP[SL]',rownames(raw_sce),ignore.case = T)]
raw_sce[["percent.mt"]] <- PercentageFeatureSet(raw_sce, pattern = "^MT-")
fivenum(raw_sce[["percent.mt"]][,1])
rb.genes <- rownames(raw_sce)[grep("^RP[SL]",rownames(raw_sce))]
C<-GetAssayData(object = raw_sce, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
raw_sce <- AddMetaData(raw_sce, percent.ribo, col.name = "percent.ribo")
plot1 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw_sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
VlnPlot(raw_sce, features = c("percent.ribo", "percent.mt"), ncol = 2)
VlnPlot(raw_sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(raw_sce, features = c("percent.ribo", "nCount_RNA"), ncol = 2)
pro='merge'
VlnPlot(raw_sce, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(raw_sce, features = c("percent.ribo", "nCount_RNA"), ncol = 2)
raw_sce <- subset(raw_sce, subset = nFeature_RNA > 200 & nCount_RNA > 1000 & percent.mt < 20)
dim(raw_sce)
dim(sce)
####################################
#merge 合并初步查看聚类结果（效果不好，未进行矫正）
sce=raw_sce
sce <- NormalizeData(sce, normalization.method =  "LogNormalize", 
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", nfeatures = 2000) 
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
DimHeatmap(sce, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.2)
table(sce@meta.data$RNA_snn_res.0.2)
set.seed(123)
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
##

mycol <- c("#223D6C","#D20A13","#FFD121","#088247",
           "#58CDD9","#5D90BA","#431A3D","#11AA4D",
           "#91612D","#6E568C","#7A142C",
           "#E0367A","#D8D155","#64495D","#7CC767")
DimPlot(sce,reduction = "tsne",label=T,cols = mycol)
DimPlot(sce,reduction = "umap",label=T,split.by ='orig.ident')
dim(sce)

#DimPlot(sce,reduction = "tsne",label=T,group.by=,cols = mycol)
#sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25,  thresh.use = 0.25)
#write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
DoHeatmap(subset(sce, downsample = 1000),features=top10$gene,size=2,
          slot = 'scale.data', raster = F,group.colors = mycol) + 
  scale_fill_gradient2(
    low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
    mid = "white",
    high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
   midpoint = 0,
    guide = "colourbar",
    aesthetics = "fill"
  )
ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'))