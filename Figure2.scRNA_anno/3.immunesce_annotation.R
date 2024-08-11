####
###导入数据
library(Seurat)
library(ggplot2)
immunesce<-readRDS("data/immunesce.rds")
#Nonimmune cells with >6000 or <200 genes or >40% mitochondrial
#genes were discarded. Immune cells with >4000 or <200
#genes or >25% mitochondrial genes were filtered out.
###check the map site
#######查看情况
DefaultAssay(immunesce) <- "integrated"
immunesce<- FindNeighbors(immunesce, reduction = "harmony", dims = 1:30)
immunesce <- FindClusters(immunesce, resolution =1)
##############################第一种方式注释 已有数据注释SCINA
library(readxl)
library(SCINA)
library(tidyr)
DefaultAssay(immunesce)<-"RNA"
exp <- as.matrix(GetAssayData(immunesce))
signatures<-read_excel("/home/data/gaoyuzhen/Projects/LM_spatial/人和小鼠单细胞marker基因数据库.xlsx",3)
signatures<-signatures[signatures$organ %in% c("Intestine"),]
table(signatures$coarseCellType)
data_long <- signatures %>%
  separate_rows(canonicalMarkers, sep = ",")
table(data_long$coarseCellType)
data_long<-data_long[data_long$coarseCellType %in% c("Immune cells","T cells","Macrophages","Granulocytes","Monocytes","Neurons","B cells","Innate lymphoid cells"),]
#data_long<-data_long[data_long$coarseCellType %in% c("Immune cells"),]
table(data_long$cellTypes)
signatures<-data_long [,c(6,7)]
colnames(signatures)<-c("cluster","gene")
signatures<-signatures[signatures$gene %in% rownames(exp) ,]
signatures<-split(signatures$gene,signatures$cluster)
lengths <- sapply(signatures, length)
# 筛选长度大于10的元素
#signatures <- signatures[lengths >1]
signatures<-signatures[-1]
names(signatures)
signatures<-signatures[-c(1,2,54,55,56,57,58,31,32,33,37,34,35,36)]
signatures<-signatures[!names(signatures) %in% c("Effector T cells", "Eosinophils" ,"Erythrocytes" ,"Exhausted T cells" ,
                            "Gamma delta T cells","Granulocytes", "Group 1 innate lymphoid cells",
                            "Group 3 innate lymphoid cells","T cells","Neurons", "Neutrophils")]

results_immune = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, 
                sensitivity_cutoff = 1, 
                rm_overlap=FALSE,
                #rm_overlap=TRUE, 
                allow_unknown=TRUE, log_file='SCINA.log')
##################################################################################################################################
SCINA_immune<-data.frame(immunesce$seurat_clusters,results_immune$cell_labels)
table(SCINA_immune$immunesce.seurat_clusters,SCINA_immune$results_immune.cell_labels)
####归类
SCINA_immune$results_immune.cell_labels<-ifelse(SCINA_immune$results_immune.cell_labels %like% "CD4","CD4",
                                                ifelse(SCINA_immune$results_immune.cell_labels %like% "CD8","CD8",
                                                       ifelse(SCINA_immune$results_immune.cell_labels %like% "helper","CD4",
                                                              ifelse(SCINA_immune$results_immune.cell_labels %like% "B cells","B",SCINA_immune$results_immune.cell_labels))))
SCINA_immune$results_immune.cell_labels<-ifelse(SCINA_immune$results_immune.cell_labels %like% "CD4","CD4",
                                                ifelse(SCINA_immune$results_immune.cell_labels %like% "Mononuclear phagocytes","Myeloid",
                                                       ifelse(SCINA_immune$results_immune.cell_labels %like% "Plasmacytoid dendritic cells","Myeloid",
                                                              ifelse(SCINA_immune$results_immune.cell_labels %like% "Innate lymphoid cells","B",SCINA_immune$results_immune.cell_labels))))

table(SCINA_immune$immunesce.seurat_clusters,SCINA_immune$results_immune.cell_labels)
###############
SCINA_immune$celltype[SCINA_immune$immunesce.seurat_clusters %in% c(3,5,12,17,18)]="CD8"
SCINA_immune$celltype[SCINA_immune$immunesce.seurat_clusters %in% c(15,18)]="MAIT"
SCINA_immune$celltype[SCINA_immune$immunesce.seurat_clusters %in% c(5,12)]="NK"
#SCINA_immune$celltype[SCINA_immune$immunesce.seurat_clusters %in% c(6,8,11)]="NK"
SCINA_immune$celltype[SCINA_immune$immunesce.seurat_clusters %in% c(7,10)]="B"
SCINA_immune$celltype[SCINA_immune$immunesce.seurat_clusters %in% c(20)]="Plasma"
SCINA_immune$celltype[SCINA_immune$immunesce.seurat_clusters %in% c(13,16)]="Myeloid"
SCINA_immune$celltype[SCINA_immune$immunesce.seurat_clusters %in% c(0,1,2,4,9,14,19)]="CD4+"
SCINA_immune$celltype[SCINA_immune$immunesce.seurat_clusters %in% c(4,19)]="Treg"
######################
SCINA_immune<-SCINA_immune[,c(2,3)]
colnames(SCINA_immune)<-c("SCINA","celltype")
immunesce<-AddMetaData(immunesce,SCINA_immune)
immunesce$SCINA<-NULL
DimPlot(immunesce,reduction = "umap",label=T,group.by = "celltype")
DimPlot(immunesce,reduction = "umap",label=T)


VlnPlot(immunesce,features = "CD3G",group.by = "celltype2")
VlnPlot(immunesce,features = "CD8A",group.by = "celltype2")
VlnPlot(immunesce,features = "CD4",group.by = "celltype2")
VlnPlot(immunesce,features = "KLRF1",group.by = "celltype2")
VlnPlot(immunesce,features = "SLC4A10",group.by = "celltype2")

VlnPlot(immunesce,features = "CD4")
VlnPlot(immunesce,features = "CD8A")
VlnPlot(immunesce,features = "KLRF1")

VlnPlot(immunesce,features = "FOXP3")
VlnPlot(immunesce,features = "CD79A")
VlnPlot(immunesce,features = "MZB1")
VlnPlot(immunesce,features = "LYZ")
VlnPlot(immunesce,features = "CD14")
VlnPlot(immunesce,features = "SLC4A10")
###############第二种方式注释 single R 注释
library(SingleR)
#HumanRNA<-HumanPrimaryCellAtlasData()
HumanRNA<-DatabaseImmuneCellExpressionData()
sce_for_SingleR <- as.SingleCellExperiment(immunesce)
cluster_sce<-immunesce$seurat_clusters
#humanImmu <- ImmGenData()
#pred.humanImmu <- SingleR(test = sce_for_SingleR, ref = humanImmu, labels =humanImmu$label.main)
#mouseRNA <- MouseRNAseqData()
pred.sce <- SingleR(test = sce_for_SingleR, ref = HumanRNA, labels =HumanRNA$label.main,clusters=factor(cluster_sce))
#pred.humanRNA <- SingleR(test = sce_for_SingleR, ref = list(BP=Blue.se, HPCA=hpca.se), labels = list(Blue.se$label.main, hpca.se$label.main)) 
#table(pbmc.hesc$labels,meta$seurat_clusters)
table(pred.sce$labels)
DimPlot(immunesce,reduction = "umap",label=T,group.by = "Major_celltype")
DimPlot(immunesce,reduction = "umap",label=T)
#saveRDS(immunesce, file = "data/immunesce.rds")
#################
library(zellkonverter)
library(Seurat)
library(hdf5r)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(paletteer) 
library(harmony)
library(cowplot)
library(ggpubr)
library(AUCell)
library(SeuratDisk)
exp <- data.frame(immunesce@assays$integrated@scale.data)
meta <-immunesce@meta.data
identical(rownames(meta), colnames(exp))
data <- CreateSeuratObject(counts = exp)
data@meta.data <- meta
data$Subset<-data$seurat_clusters
table(data$Subset)
data$Subset <- as.character(data$Subset)
SaveH5Seurat(data,filename = "immunesce.h5Seurat", overwrite = TRUE)
Convert("immunesce.h5Seurat",dest = "h5ad", overwrite = TRUE)
##############################第三种方式注释 自建模型python celltypist注释
###########注释分析 python 
## celltypist python 运行
exp <- data.frame(sce@assays$RNA@counts)
meta <-sce@meta.data
exp<-exp[rownames(exp) %in% rownames(GSE225857scesub),]
identical(rownames(meta), colnames(exp))
data <- CreateSeuratObject(counts = exp)
data@meta.data <- meta
SaveH5Seurat(data,filename = "LM.h5Seurat", overwrite = TRUE)
Convert("LM.h5Seurat",dest = "h5ad", overwrite = TRUE)
##############
metadata<-read.csv("/home/data/gaoyuzhen/Projects/LM_spatial/adata_obs_with_predictions_major_immune_match.csv",row.names = 1)
rownames(metadata)<-gsub("\\.","-",rownames(metadata))
metadata$sci_adv_celltype[metadata$Major_celltype!="T/ILC"]<-metadata$Mini_celltype[metadata$Major_celltype!="T/ILC"]
metadata$sci_adv_celltype[metadata$Mini_celltype=="MAIT"]<-metadata$Mini_celltype[metadata$Mini_celltype=="MAIT"]
metadata<-metadata[,-c(1:9)]
metadata<-metadata[,-c(7:21)]
immunesce<-AddMetaData(immunesce,metadata)
#######
table(immunesce$seurat_clusters,immunesce$majority_voting)
##
VlnPlot(immunesce,features = "CD8A")
VlnPlot(immunesce,features = "CD8A",group.by ="sci_adv_celltype" )
VlnPlot(immunesce,features = "SLC4A10")

immunesce$sci_adv_celltype[immunesce$seurat_clusters %in% c(0,1,3,5,10,13)]="CD4"
immunesce$sci_adv_celltype[immunesce$seurat_clusters %in% c(2,4,7,9,12)]="CD8"
immunesce$sci_adv_celltype[immunesce$seurat_clusters %in% c(9)]="MAIT"
immunesce$sci_adv_celltype[immunesce$seurat_clusters %in% c(6)]="Myeloid"
immunesce$sci_adv_celltype[immunesce$seurat_clusters %in% c(8,11)]="B"
immunesce$sci_adv_celltype[immunesce$seurat_clusters %in% c(14)]="Plasma"
immunesce$sci_adv_celltype[immunesce$majority_voting %like% "Treg"]="Treg"
immunesce$sci_adv_celltype[immunesce$Mini_celltype=="NK"]="NK"
table(immunesce$seurat_clusters,immunesce$sci_adv_celltype)
table(immunesce$sci_adv_celltype)
#################################################################################################
p1<-DimPlot(immunesce,reduction = "umap",label=T,group.by = "sci_adv_celltype")
p2<-DimPlot(immunesce,reduction = "umap",label=T)
p1+p2


#########################################
sce.markers <- FindAllMarkers(immunesce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.sig<- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)# top 10 
write.csv(sce.markers.sig,file = "results/sce.markers.sig.csv")
table(sce$seurat_clusters)
########################################
VlnPlot(immunesce,features ="FOXP3",group.by = "sub_celltype")

######################################################################
####合并到原始文件中等所有数据注释完毕后合并
############################################
metadata<-immunesce@meta.data
colnames(metadata)
immunecellsmetadata<-metadata[,-c(2:9,21:34)]
immunecellsmetadata<-immunecellsmetadata[,-12]
colnames(immunecellsmetadata)[12]<-"sub_celltype"
##############################################
load("data/totalscRNAdata.Rdata")
metadata<-sce@meta.data
nonimmunecellsmetadata<-metadata[metadata$immunecells=="Nonimmunecells",]
colnames(nonimmunecellsmetadata)
nonimmunecellsmetadata<-nonimmunecellsmetadata[,-c(13,14)]
metadata<-rbind(immunecellsmetadata,nonimmunecellsmetadata)
table(metadata$immunecells,metadata$sub_celltype)
save(metadata,file ="data/totalscRNAdata.Rdata")

table(metadata$seurat_clusters,metadata$sub_celltype)

############## 导出 计算compass 数据
library(readr)
exp_CRC<-c()
exp_CRC<-immunesce@assays$RNA@counts
exp_CRC<-as.matrix(exp_CRC)
exp_CRC<-cbind(symbol=rownames(exp_CRC),exp_CRC)
exp_CRC[1:5,1:5]
exp_CRC<-as.data.frame(exp_CRC)
write_tsv(exp_CRC,file="/home/data/gaoyuzhen/Projects/LM_spatial/Compass/exp_LM_CRC_immunecells.tsv")
