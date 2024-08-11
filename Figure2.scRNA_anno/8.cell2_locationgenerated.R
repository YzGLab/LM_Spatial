###
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
library(GEOquery)
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(biomaRt)
library(GSVA)


## 单细胞参考矩阵
## all for Leading edge sample

sce$Subset<-sce$Major_celltype
table(sce$Subset)
c2lLM<-sce[,sce$set=="LM"]
exp <- data.frame(c2lLM@assays$RNA@counts)
meta <- c2lLM@meta.data
identical(rownames(meta), colnames(exp))
data <- CreateSeuratObject(counts = exp)
data@meta.data <- meta
table(data$Subset)
data$Subset <- as.character(data$Subset)
SaveH5Seurat(data,filename = "LM_scRNA_for_C2L.h5Seurat", overwrite = TRUE)
Convert("LM_scRNA_for_C2L.h5Seurat",dest = "h5ad", overwrite = TRUE)

#########################
myexp<-exp
myexp<-as.matrix(myexp)
hsmart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart = "ensembl",mirror = 'asia')
mapping <- getBM(
  attributes = c("entrezgene_id",'ensembl_gene_id', 'ensembl_transcript_id','hgnc_symbol'), 
  filters = "hgnc_symbol",###########根据指标id,基因名称或者ensembl 等筛选 
  values = rownames(myexp),
  mart = hsmart
)
mapping<-unique(mapping)
save(mapping,file = "mapping_gene.Rdata")
mapping<-mapping[,c(2,4)]
mapping<-unique(mapping)
probe2symbol<-mapping
colnames(probe2symbol)<-c("symbol",'GeneID-2')
rownames(probe2symbol)<-probe2symbol[,1]

myexp<-data.frame(myexp)
#myexp <- myexp %>% 
 # rownames_to_column(var="probeset") %>% 
  #inner_join(probe2symbol,by="probeset") %>% 
 # dplyr::select(-probeset) %>% 
 # dplyr::select(symbol,everything()) %>% 
 # mutate(rowMean =rowMeans(.[,-1])) %>% 
 # filter(symbol != "NA") %>% 
  #arrange(desc(rowMean)) %>% 
  #distinct(symbol,.keep_all = T) %>% 
  #dplyr::select(-rowMean) %>% 
 # column_to_rownames(var = "symbol")


####
library(Seurat)
library(ggplot2)
library(ggsci)
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
##spatial  tiss distubution
LM_data <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir ='/home/data/gaoyuzhen/Projects/LM_spatial/spatial/LM_1',
  filename = "filtered_feature_bc_matrix.h5",
  #assay = "Spatial", # specify name of the initial assay
  #slice = "spatial", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)

SaveH5Seurat(LM_data,filename = "LM_spatial_for_C2L.h5Seurat", overwrite = TRUE)
Convert("LM_spatial_for_C2L.h5Seurat",dest = "h5ad", overwrite = TRUE)
