###导入所有数据
library(Seurat)
sce<-readRDS("data/Colon_HC_combined.sce_filtered.rds")
immunesce<-readRDS("data/immunesce.RDS")


###############查看初步情况并生成 nonimmunesce
nonimmunesce<-sce[,!colnames(sce) %in% colnames(immunesce)]
DefaultAssay(nonimmunesce) <- "integrated"
nonimmunesce<- FindNeighbors(nonimmunesce, reduction = "harmony", dims = 1:30)
nonimmunesce <- FindClusters(nonimmunesce, resolution =0.5)
nonimmunesce <- RunUMAP(nonimmunesce, reduction = "harmony", dims = 1:30)
nonimmunesce<- RunTSNE(nonimmunesce, reduction = "harmony", dims = 1:30)
DimPlot(nonimmunesce,reduction = "umap",label=T,split.by="orig.ident")
DimPlot(nonimmunesce,reduction = "tsne",label=T,split.by="orig.ident")
saveRDS(nonimmunesce, file = "data/nonimmunesce.rds")
#########################################################################################################################
#######################导入数据
nonimmunesce<-readRDS("data/nonimmunesce.rds")
###############第二种方式注释 single R 注释##single R 不是很准 重新弄 
library(SingleR)
HumanRNA<-HumanPrimaryCellAtlasData()
#HumanRNA<-DatabaseImmuneCellExpressionData()
sce_for_SingleR <- GetAssayData(nonimmunesce, slot="data")
#sce_for_SingleR <- as.SingleCellExperiment(sce)
cluster_sce<-nonimmunesce$seurat_clusters
#humanImmu <- ImmGenData()
#pred.humanImmu <- SingleR(test = sce_for_SingleR, ref = humanImmu, labels =humanImmu$label.main)
#mouseRNA <- MouseRNAseqData()
pred.sce <- SingleR(test = sce_for_SingleR, ref = HumanRNA, labels =HumanRNA$label.main,clusters=factor(cluster_sce))
#pred.humanRNA <- SingleR(test = sce_for_SingleR, ref = list(BP=Blue.se, HPCA=hpca.se), labels = list(Blue.se$label.main, hpca.se$label.main)) 
#table(pbmc.hesc$labels,meta$seurat_clusters)
table(pred.sce$labels)
celltype <-data.frame(ClusterID=rownames(pred.sce),   pred.sce$labels)    
nonimmunesce[['celltype']]<-celltype[,2][match(Idents(nonimmunesce), celltype$ClusterID)]
DimPlot(nonimmunesce,reduction = "umap",label=T,split.by="orig.ident")
DimPlot(nonimmunesce,reduction = "umap",label=T,split.by="orig.ident",group.by = "celltype")
VlnPlot(nonimmunesce,features = "COL1A1",pt.size = 0)
VlnPlot(nonimmunesce,features = "ABCC9",pt.size = 0)
colnames(celltype)<-c("seurat_clusters","celltype")
celltype$celltype[celltype$seurat_clusters %in% c(2,3,7,9,10)]<-"Stomal cells"
nonimmunesce[['celltype']]<-celltype[,2][match(Idents(nonimmunesce), celltype$seurat_clusters)]

##############可视化
DimPlot(nonimmunesce,reduction = "umap",label=T,split.by="orig.ident",group.by = "celltype")
###################第二种方式注释 已有数据注释SCINA
nonimmunesce<-readRDS("data/nonimmunesce.rds")
library('SCINA')
library('preprocessCore')
#exp = nonimmunesce@assays$RNA@data
#exp<-as.matrix(exp)
#exp_raw=log(exp+1)
#exp[]=normalize.quantiles(exp_raw)
#exp[1:5,1:5]
exp <- as.matrix(GetAssayData(nonimmunesce))
# Or .csv examples
library(readxl)
signatures<-read_excel("/home/data/gaoyuzhen/Projects/LM_spatial/人和小鼠单细胞marker基因数据库.xlsx",3)
signatures<-signatures[signatures$organ %in% c("Intestine"),]
table(signatures$coarseCellType)
library(tidyr)
data_long <- signatures %>%
  separate_rows(canonicalMarkers, sep = ",")
table(data_long$coarseCellType)
data_long<-data_long[data_long$coarseCellType %in% c("Endothelial cells","Epithelial cells","Stem cells","Stromal cells","Glial cells","Neuroendocrine cells"),]
table(data_long$cellTypes)
signatures<-data_long [,c(6,7)]
colnames(signatures)<-c("cluster","gene")
signatures<-signatures[signatures$gene %in% rownames(exp) ,]
signatures<-split(signatures$gene,signatures$cluster)
lengths <- sapply(signatures, length)
# 筛选长度大于10的元素
signatures <- signatures[lengths > 5]
signatures<-signatures[-15]
#########################################
results = SCINA(exp, signatures, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.999, sensitivity_cutoff = 1, 
                rm_overlap=FALSE,
                #rm_overlap=TRUE, 
                allow_unknown=TRUE, log_file='SCINA.log')

SCINA<-data.frame(nonimmunesce$seurat_clusters,results$cell_labels)
table(SCINA$nonimmunesce.seurat_clusters,SCINA$results.cell_labels)
SCINA$celltype[SCINA$nonimmunesce.seurat_clusters %in% c(0,1,5,6,8,10,12)]="Endothelial"
SCINA$celltype[SCINA$nonimmunesce.seurat_clusters %in% c(2,3,7,9)]="Fibroblasts"
SCINA$celltype[SCINA$nonimmunesce.seurat_clusters %in% c(4,11)]="Epithelial"
#SCINA$celltype[SCINA$nonimmunesce.seurat_clusters %in% c(9)]="Unknown"
SCINA$celltype[SCINA$nonimmunesce.seurat_clusters %in% c(13)]="Neuroendocrine cells"
SCINA$celltype[SCINA$nonimmunesce.seurat_clusters %in% c(10) & SCINA$results.cell_labels=="Fibroblasts"]="Fibroblasts"
SCINA<-SCINA[,c(2,3)]
colnames(SCINA)<-c("SCINA","celltype")
nonimmunesce<-AddMetaData(nonimmunesce,SCINA)
nonimmunesce$SCINA<-NULL
##可视化
DimPlot(nonimmunesce,reduction = "umap",label=T)
DimPlot(nonimmunesce,reduction = "umap",label=T,group.by = "celltype")
nonimmunesce$set<-ifelse(nonimmunesce$orig.ident %in% c("chang1","chang2"),"CRC","LM")
DimPlot(nonimmunesce,reduction = "umap",label=T,group.by = "celltype",split.by = "set")
saveRDS(nonimmunesce, file = "data/nonimmunesce.rds")
mycol<-c('#1B9E77','#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0','#F0027F', '#BF5B17', '#E41A1C', 
         '#D95F02','#7570B3','#E7298A', '#66A61E','#E6AB02','#A6761D','#666666',
         '#377EB8','#4DAF4A', '#984EA3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999',
         '#8DD3C7', '#FFFFB3','#BEBADA', '#FB8072','#80B1D3','#FDB462','#B3DE69', '#FCCDE5',
         '#D9D9D9','#BC80BD','#CCEBC5', '#FFED6F','#66C2A5', '#FC8D62','#8DA0CB','#E78AC3',
         '#A6D854','#FFD92F', '#E5C494', '#B3B3B3','#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
         '#FB9A99', '#E31A1C','#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928')
p.SCINA2 = DimPlot(nonimmunesce, reduction = "umap",split.by="orig.ident",cols = mycol,
                   group.by = "SCINA",
                   label = TRUE, label.size = 3, repel = TRUE)
p.SCINA2
p0<-DimPlot(nonimmunesce,reduction = "umap",label=T,split.by="orig.ident",group.by = "celltype")
p0+p.SCINA2
###

##############################第三种方式注释 自建模型python celltypist注释
exp <- data.frame(nonimmunesce@assays$RNA@counts)
meta <-nonimmunesce@meta.data
identical(rownames(meta), colnames(exp))
data <- CreateSeuratObject(counts = exp)
data@meta.data <- meta
data$Subset<-data$seurat_clusters
table(data$Subset)
data$Subset <- as.character(data$Subset)
SaveH5Seurat(data,filename = "nonimmunesce.h5Seurat", overwrite = TRUE)
Convert("nonimmunesce.h5Seurat",dest = "h5ad", overwrite = TRUE)
#################################差异基因可视化
library(scales)
library(tidyr)
library(tidyverse)
Idents(nonimmunesce)<-nonimmunesce$seurat_clusters
sce.markers <- FindAllMarkers(nonimmunesce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.sig <- sce.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)# top 10 
sce.markers.sig.9<-sce.markers.sig[sce.markers.sig$cluster==9,7]
sce.markers.sig.9$gene
genes_list <- paste(sce.markers.sig.9$gene, collapse = ",")
print(genes_list)
#######################
DotPlot(sce,features = unique(top10$gene),group.by = "metabolism_cluster") + RotatedAxis()# +coord_flip()
#ggsave("fibresults/DotPlot_withtogenes.pdf",width = 8,height = 3.68)
#
sce.markers.sig <- sce.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)# top 10 
sce.markers.sig<-sce.markers.sig[sce.markers.sig$p_val_adj<0.05,]
nonimmunesce$SCINA <- results$cell_labels


#sce_fj<-subset(x = sce_fj, downsample = 5000)
subsce<-nonimmunesce[,nonimmunesce$celltype=="Endothelial_cells"]
############compass
nonimmunesce<-readRDS("data/nonimmunesce.rds")
library(readr) ## for compass
exp_CRC<-c()
exp_CRC<-nonimmunesce@assays$RNA@counts
exp_CRC<-as.matrix(exp_CRC)
exp_CRC<-cbind(symbol=rownames(exp_CRC),exp_CRC)
exp_CRC[1:5,1:5]
exp_CRC<-as.data.frame(exp_CRC)
write_tsv(exp_CRC,file="/home/data/gaoyuzhen/Projects/LM_spatial/Compass/exp_LM_CRC_nonimmunecells.tsv")
