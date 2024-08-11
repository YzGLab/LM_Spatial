#####merge them in sce。 
###########
library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)
##
############
stcol<-c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#f5e801', '#08f7f0')
st_sce<-readRDS("data/Colon_HC_spatial.RDS")


sample<-"LM_2" 
sample<-"LM_1" 
sample<-"C_1" 
sample<-"C_2" 
##############
#list.files("/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/stRNA_scMetabolismuniondata")
load(paste0("~/Projects/LM_spatial/ResultsLM/Figure7_sp/",sample,"sample_classification.Rdata"))
sp<-readRDS(paste0("/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/stRNA_scMetabolismuniondata/",sample,"_both.rds"))
##########################
#匹配数据 
########################################################################################################
sp_matraix<-sp$metadata
rownames(sp_matraix)<-sp_matraix$mz
#############################
table(st_sce$orig.ident)
sp_matraix_filter<-sp_matraix[,c(15:ncol(sp_matraix))]
colnames(sp_matraix_filter)<-paste0(colnames(sp_matraix_filter),"-1")
colnames(sp_matraix_filter)<-gsub("\\.","_",colnames(sp_matraix_filter))

###############################################导入空间代谢组学分群数据
sp_cluster<-data.frame(sample_classification)
#sp_cluster<-cbind(sample=sp$samplename$transname,clusters)
sp_cluster$Sample<-paste0(sp_cluster$Sample,"-1")
sp_cluster$Sample<-gsub("\\.","_",sp_cluster$Sample)
rownames(sp_cluster)<-sp_cluster$Sample
colnames(sp_cluster)[c(2:ncol(sp_cluster))]<-c(3,6,8,12)
colnames(sp_cluster)[c(2:ncol(sp_cluster))]<-paste0("cluster",colnames(sp_cluster)[c(2:ncol(sp_cluster))])

#################匹配样本 
sub_sce<-st_sce[,st_sce$orig.ident==sample]
sub_sce<-sub_sce[,colnames(sub_sce) %in% sp_cluster$Sample]

#############建立差异分析的assay 
sp_matraix_filter<-sp_matraix_filter[,colnames(sp_matraix_filter) %in% sp_cluster$Sample]
sp_matraix_assay <- CreateAssayObject(counts =sp_matraix_filter)
sub_sce[["metabolits"]] <- sp_matraix_assay
##########添加信息
sub_sce<-AddMetaData(sub_sce,metadata = sp_cluster)


sub_sce.list<-list()
sample_names <- c("LM_1", "LM_2", "C_1", "C_2")

for (sample in sample_names) {
  # 读取数据
  print(sample)
  load(paste0("~/Projects/LM_spatial/ResultsLM/Figure7_sp/",sample,"sample_classification.Rdata"))
  sp<-readRDS(paste0("/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/stRNA_scMetabolismuniondata/",sample,"_both.rds"))
  sp_matraix<-sp$metadata
  rownames(sp_matraix)<-sp_matraix$mz
  #############################
  table(st_sce$orig.ident)
  sp_matraix_filter<-sp_matraix[,c(15:ncol(sp_matraix))]
  colnames(sp_matraix_filter)<-paste0(colnames(sp_matraix_filter),"-1")
  colnames(sp_matraix_filter)<-gsub("\\.","_",colnames(sp_matraix_filter))
  
  ###############################################导入空间代谢组学分群数据
  sp_cluster<-data.frame(sample_classification)
  #sp_cluster<-cbind(sample=sp$samplename$transname,clusters)
  sp_cluster$Sample<-paste0(sp_cluster$Sample,"-1")
  sp_cluster$Sample<-gsub("\\.","_",sp_cluster$Sample)
  rownames(sp_cluster)<-sp_cluster$Sample
  colnames(sp_cluster)[c(2:ncol(sp_cluster))]<-c(3,6,8,12)
  colnames(sp_cluster)[c(2:ncol(sp_cluster))]<-paste0("cluster",colnames(sp_cluster)[c(2:ncol(sp_cluster))])
  
  #################匹配样本 
  sub_sce<-st_sce[,st_sce$orig.ident==sample]
  sub_sce<-sub_sce[,colnames(sub_sce) %in% sp_cluster$Sample]
  
  #############建立差异分析的assay 
  sp_matraix_filter<-sp_matraix_filter[,colnames(sp_matraix_filter) %in% sp_cluster$Sample]
  sp_matraix_assay <- CreateAssayObject(counts =sp_matraix_filter)
  sub_sce[["metabolits"]] <- sp_matraix_assay
  ##########添加信息
  sub_sce<-AddMetaData(sub_sce,metadata = sp_cluster)
  sub_sce.list[[sample]]<-sub_sce
}


sce.big <- merge(sub_sce.list[[1]], 
                 y = c(sub_sce.list[[2]],sub_sce.list[[3]],sub_sce.list[[4]]), 
                 #add.cell.ids = samples, 
                 project = "spatial_4")
table(sce.big$orig.ident)
DefaultAssay(sce.big)  = 'metabolits'
sce.big <- ScaleData(sce.big,verbose = FALSE)  
sce.big<- FindVariableFeatures(sce.big) 
sce.big<-RunPCA( sce.big,npcs = 30,verbose = FALSE)  
sce.big <- RunHarmony(sce.big,group.by.vars = "orig.ident")
library(harmony)
sce<- FindNeighbors(sce.big,reduction = "harmony", dims = 1:20)
sce <- FindClusters(sce,reduction = "harmony",resolution = 0.4)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:30)
sce<- RunTSNE(sce, reduction = "harmony",dims = 1:30)
saveRDS(sce,file="data/metabolits_scp_st_sp_sce.RDS")


sce<-readRDS("data/metabolits_scp_st_sp_sce.RDS")
  ###############
library(SCP)
colnames(sce@meta.data)
DimPlot(sce,reduction = "pca",label=T,group.by="orig.ident",)
DimPlot(sce,reduction = "umap",group.by="orig.ident",label=T)
DimPlot(sce,reduction = "umap",label=T,split.by="set",group.by="st_celltype",cols = stcol1)
DimPlot(sce,reduction = "umap",label=T)
stcol1<-c()
#names(stcol1)<- c('Normal Epi', 'Tumor', 'Fibroblast', 'Hepatocytes', 'Lamina propria', 'B/Plasma', 'Monocyte')
stcol1 <- c(
            'Tumor' = '#ff7f0e', 
            'Normal Epi' = '#1f77b4', 
            'Hepatocytes' = '#d62728', 
            'Fibroblast' = '#2ca02c', 
            'B/Plasma' = '#f5e801',
            'Monocyte' = '#08f7f0',
            'Lamina propria' = '#9467bd'
            )

stcol1<- setNames(c('#1f77b4','#ff7f0e','#2ca02c',  '#d62728', '#9467bd','#f5e801','#08f7f0'),
                       c('Normal Epi', 'Tumor', 'Fibroblast', 'Hepatocytes', 'Lamina propria', 'B/Plasma', 'Monocyte'))


CellDimPlot(
  srt = sce,
  split.by = "set",
  #add_mark = TRUE,
  cells.highlight = TRUE,
  stat_plot_size = 0.1,
  group.by="st_celltype",
  palcolor = stcol1,
  label = T,
  reduction = "UMAP",
  theme_use = "theme_blank"
)

ggsave(file = file.path(Figure7,"ST_Celltype_CRC_LM_split_sp.pdf"),width = 11.5,height = 6)
