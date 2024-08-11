library(Seurat)
library(ggplot2)
library(pbapply)###可视化 百分比进度
library(SCP)
##st RNA
sub_sce<-readRDS("data/spatial_combinend_removed_MT_RP.RDS")
####################

Idents(sub_sce)<-sub_sce$st_celltype
sce.markers1 <- FindAllMarkers(sub_sce[,sub_sce$set=="CRC"], 
                               only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
sce.markers2 <- FindAllMarkers(sub_sce[,sub_sce$set=="LM"], 
                               only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
#sce.markers.sig <- sce.markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10

sce.markers1$cluster<-paste0(sce.markers1$cluster,".CRC")
sce.markers2$cluster<-paste0(sce.markers2$cluster,".LM")
sce.markers<-rbind(sce.markers1,sce.markers2)################################
sce.cell_type_markers<-split(sce.markers$gene,sce.markers$cluster)###############合并两种差异
save(sce.cell_type_markers,file='data/spatial_differetgenes_clusters.Rdata')
#####################
###############scRNA-fib
##
library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)
library(scales)
library(tidyr)
library(tidyverse)
library(patchwork)
##############
sce<-readRDS("data/Colon_HC_combined.sce.rds")
load("data/totalscRNAdata.Rdata")
###
sce<-AddMetaData(sce,metadata)
mycol2<-c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F", "#E41A1C", '#FFD92F', "#C3BC3F", "#FF7F00", "#1B9E77", "#66A61E", "#F1788D")

mycol2<-c( '#A65628', '#984EA3','#FDB462','#BEAED4','#A6D854','#1B9E77','#FFD92F','#E5C494','#8DA0CB',
           '#E41A1C','#FF7F00','#377EB8','#E78AC3',
           '#FFFF33', '#A65628', '#F781BF', '#999999',
           '#8DD3C7', '#BEBADA', '#FB8072','#80B1D3','#B3DE69', '#FCCDE5')
###############################
table(sce$sub_celltype)
fib_sce<-sce[,sce$sub_celltype=="Fibroblasts"]


CellDimPlot(
  srt =sce, group.by =  "sub_sub_celltype",  stat_plot_size = 0.1,
  reduction = "UMAP", theme_use = "theme_blank", label = T,
  palcolor = mycol2,
)

fib_sce <- FindNeighbors(fib_sce, reduction = "harmony", dims = 1:30)
fib_sce <- FindClusters(fib_sce, resolution = 0.2)
fib_sce <- RunUMAP(fib_sce,reduction = "harmony",  dims = 1:20)
fib_sce <- RunUMAP(fib_sce,reduction = "harmony", dims = 1:30)
##################
colnames(fib_sce@meta.data)
DefaultAssay(fib_sce)<-"integrated"
Idents(fib_sce)<-fib_sce$integrated_snn_res.0.2
#########################################
sce.markers<- FindAllMarkers(fib_sce, 
                                only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.x.sig<- sce.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)# top 10 

CellDimPlot(
  srt =fib_sce, group.by =  "seurat_clusters",  stat_plot_size = 0.1,
  reduction = "UMAP", theme_use = "theme_blank", label = T,split.by = "set",
  palcolor = mycol2,
)

ggsave(file = file.path(Figure4,"celltype_immune_CRC_LM_split.pdf"),width = 9.5,height = 7.03)
##lable them in marker genes
sce.markers.x.sig.10<- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)# top 10
sce.markers.x.sig.10

###


FeatureDimPlot(
  srt = sce, features = c("COL1A1","ABCC9", "MYH11","THBS2", "PKHD1","TOP2A"),ncol=3,
  reduction = "UMAP", theme_use = "theme_blank"
)

table(sce$sub_sub_celltype,sce$immunecells)
############

metadata<-fib_sce@meta.data
colnames(metadata)
metadata$sub_celltype_genes<-c()
############################################################################################################
metadata$sub_celltype_genes[metadata$seurat_clusters=="0"]<-"Fib(C1)-ABCC9"
metadata$sub_celltype_genes[metadata$seurat_clusters=="1"]<-"Fib(C2)-MYH11"
metadata$sub_celltype_genes[metadata$seurat_clusters=="2"]<-"Fib(C3)-THBS2"
metadata$sub_celltype_genes[metadata$seurat_clusters=="3"]<-"Fib(C4)-SPP1"
metadata$sub_celltype_genes[metadata$seurat_clusters=="4"]<-"Fib(C5)-TOP2A"
fib_sce<-AddMetaData(fib_sce,metadata)

fib_sce$sub_celltype_genes<-factor(fib_sce$sub_celltype_genes,
                                        levels=c("Fib(C1)-ABCC9","Fib(C2)-MYH11", "Fib(C3)-THBS2", 
                                                 "Fib(C4)-SPP1", "Fib(C5)-TOP2A")
)

CellDimPlot(
  srt =fib_sce, group.by =  "sub_celltype_genes",  stat_plot_size = 0.1,
  reduction = "UMAP", theme_use = "theme_blank", label = T,split.by = "set",
  palcolor = mycol2,
)
ggsave(file = file.path(Figure4,"celltype_immune_CRC_LM_split.pdf"),width = 9.5,height = 7.03)
FeatureDimPlot(
  srt = fib_sce, features = c("COL1A1","ABCC9", "MYH11","THBS2", "CXCL8","TOP2A"),ncol=3,
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave(file = file.path(Figure4,"celltype_immune_CRC_LM_split_markergenes.pdf"),width = 9.5,height = 7.03)


saveRDS(fib_sce,file="data/fib_sce_scRNA.Rds")

FeatureDimPlot(
  srt = fib_sce, features = c("COL1A1","ABCC9", "MYH11","THBS2", "CXCL8","TOP2A"),ncol=3,
compare_features = TRUE, label = TRUE, label_insitu = TRUE,
add_density=TRUE,
palcolor=mycol2,
#label_segment_color =mycol2,
lineages_palcolor = mycol2,
reduction = "UMAP", theme_use = "theme_blank"
)

#########################################
sce.markers<- FindAllMarkers(fib_sce, 
                             only.pos = TRUE, 
                             min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.x.sig<- sce.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)# top 10 
write.csv(sce.markers,file.path(Figure4,"Fibroblasts_markergenes_split.csv"))
write.csv(sce.markers.x.sig,file=file.path(Figure4,"Fibroblasts_markergenes_split_sce.markers_sig.csv"))
#############################
fib.sce.markers.x.sig<-split(sce.markers.x.sig$gene,sce.markers.x.sig$cluster)
###存储文件
save(sce.markers.x.sig,fib.sce.markers.x.sig,file=file.path(Figure4,"Fibroblasts_markergenes_split_markers_sig.Rdata"))

#####################
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
########
DefaultAssay(st_sce)<-"SCT"
#############评分 看表达差异  
st_sce.copy<-AddModuleScore(st_sce,features=fib.sce.markers.x.sig,name=names(fib.sce.markers.x.sig))

##################


stcol<-c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#f5e801', '#08f7f0')
colnames(st_sce.copy@meta.data)[c(32:36)]<-c("Fib(C1)-ABCC9","Fib(C2)-MYH11", "Fib(C3)-THBS2", 
                                             "Fib(C4)-SPP1", "Fib(C5)-TOP2A")
VlnPlot(st_sce.copy,features = c("Fib(C1)-ABCC9","Fib(C2)-MYH11", "Fib(C3)-THBS2", 
                                 "Fib(C4)-SPP1", "Fib(C5)-TOP2A"),
        #stack = T, 
        #fill.by='ident', 
        cols =stcol, pt.size=0)
table(st_sce.copy$st_celltype)
VlnPlot(st_sce.copy[,st_sce.copy$st_celltype=="Fibroblast"],features = c("Fib(C1)-ABCC9","Fib(C2)-MYH11", "Fib(C3)-THBS2", 
                                                                         "Fib(C4)-SPP1", "Fib(C5)-TOP2A"),
        group.by = "st_celltype",split.by = "set",ncol=5,
        #stack = T, 
        #fill.by='ident', 
        cols =c('#1f77b4','#d62728'), pt.size=0)
ggsave(file = file.path(Figure4,"Fibroblast_score_CRC_LM_split_markergenes.pdf"),width =11.5,height = 3.53)
######################
library(raster)
fig4a <-DotPlot(fib_sce,features = unique(sce.markers.x.sig.10$gene),group.by = "sub_celltype_genes") + 
  #rotate() +
  scale_color_gradientn(colours=brewer.pal(n=8, name="PuBuGn"), guide = "colourbar") +
  theme_minimal_grid() + theme(axis.text.x=element_text(angle=30, hjust=1))+coord_flip()
fig4a

ggsave(file.path(Figure4,"Fibroblasts_markergenes_split.pdf"),width =5.2, height =8.5)
###
##组织偏好性分析
library("sscVis")
library(plyr)
library(ldply)
library(RColorBrewer)
library(data.table)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
#加载函数及数据
source("/home/data/gaoyuzhen/Projects/LiverCancer/Metabolism_code/test_function.R")
source("/home/data/gaoyuzhen/Projects/LiverCancer/Metabolism_code/draw_analysis.R")
#############################################################################################
colnames(fib_sce@meta.data)
Tumor <- fib_sce
Tumor$sub_sub_celltype=factor(Tumor$sub_sub_celltype,
                              levels=c("C0-Fibroblasts", "C1-Fibroblasts","C2-Fibroblasts","C3-Fibroblasts",
                                       "C4-Fibroblasts", "C5-Fibroblasts",  "C6-Fibroblasts", "C8-Fibroblasts", 
                                       "C7-Fibroblasts"))

#数据分析
A <- do.tissueDist(cellInfo.tb = Tumor@meta.data,#这里的输入文件需要的是整理好的有分组和细胞类型的metadata文件
                   out.prefix = "/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/",#设定导出的文件名，自己设置
                   pdf.width = 5,#热图的宽设置
                   pdf.height = 8,#热图的高度设置
                   verbose=1,#设置为1文件以list形式存储
                   meta.cluster = 'sub_celltype_genes',#这里是细胞类型，也可以是seurat_clusters，名称没有要求，就是你细胞类型的列名
                   loc = 'set',#这里就是分组，metadata中分组的列名，至于命名没有要求
                   z.hi=2,
                   "Fibroblasts")#热图legend最大值，根据实际情况自己设置
###