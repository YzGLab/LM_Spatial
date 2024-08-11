####
#read count 
GSE225857sce<-fread("/home/data/gaoyuzhen/Projects/LM_spatial/data/GSE225857_sciadv/GSM7058754_immune_counts.txt.gz")
GSE225857sce<-data.frame(GSE225857sce)
GSE225857sce[1:5,1:5]
rownames(GSE225857sce)<-GSE225857sce[,1]
GSE225857sce<-GSE225857sce[,-1]

##################################################################
GSE225857sce<-CreateSeuratObject(counts =GSE225857sce, project = "GSE225857sce" )
GSE225857metadata<-read.delim2("/home/data/gaoyuzhen/Projects/LM_spatial/data/GSE225857_sciadv/GSM7058754_immune_meta.txt.gz",row.names = 1)
rownames(GSE225857metadata)<-gsub("\\-",".",rownames(GSE225857metadata))
table(GSE225857metadata$cluster)
table(GSE225857metadata$orig.ident)
colnames(GSE225857sce)
GSE225857sce<-AddMetaData(GSE225857sce,GSE225857metadata)

GSE225857sce<-GSE225857sce[,GSE225857sce$orig.ident!="NA"]
exp<-fread("/home/data/gaoyuzhen/Projects/LM_spatial/data/GSE225857 /GSM7058754_immune_counts.txt.gz") %>% data.frame()
exp <- data.frame(GSE225857sce@assays$RNA@counts)
meta <-GSE225857sce@meta.data
identical(rownames(meta), colnames(exp))
#GSE225857sce <- CreateSeuratObject(counts = GSE225857sce)
table(GSE225857sce$organs)
GSE225857scesub<-GSE225857sce[,GSE225857sce$organs %in% c("CCL","LCL")]

table(GSE225857scesub$cluster)
GSE225857scesub$celltype<-str_split(GSE225857scesub$cluster,"_") %>% sapply(., "[[", 1) 
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
library(stringr)
# 自定义函数
get_element <- function(x) {
  if (length(x) >= 2) {
    return(x[2])  # 如果列表有两个或更多元素，返回第二个元素
  } else {
    return(x[1])  # 如果列表只有一个元素，返回第一个元素
  }
}

# 应用自定义函数
GSE225857scesub$celltype <- sapply(str_split(GSE225857scesub$cluster, "_"), get_element)
GSE225857scesub
table(GSE225857scesub$celltype)
??SaveH5Seurat
SaveH5Seurat(GSE225857scesub,filename = "LM_GSE225857.h5Seurat", overwrite = TRUE)
Convert("LM_GSE225857.h5Seurat",dest = "h5ad", overwrite = TRUE)
