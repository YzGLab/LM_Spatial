### compareGroups+moonBook+table1sci
#spatial for lm
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
Colon_data <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir ='/home/data/gaoyuzhen/Projects/LM_spatial/spatial/aggr',
  filename = "filtered_feature_bc_matrix.h5",
  #assay = "Spatial", # specify name of the initial assay
  #slice = "spatial", # specify name of the stored image
  filter.matrix = TRUE, 
  image = NULL,
  to.upper = FALSE
)

dim(x = Colon_data) # the number of features (genes) by samples (spots)
nrow(x = Colon_data) # the number of features
ncol(x = Colon_data) # the number of samples
#Extract the feature and sample names using rownames or colnames.
head(x = rownames(Colon_data), n = 5)
tail(x = colnames(Colon_data), n = 5)
##


colnames(Colon_data[[]]) # automatically calculated while creating the Seurat object
## [1] "orig.ident"       "nCount_Spatial"   "nFeature_Spatial"
head(Colon_data@meta.data)
Colon_data$nCount_Spatial[1:3]
## AAACAAGTATCTCCCA-1 AAACACCAATAACTGC-1 AAACAGCTTTCAGAAG-1 
##              14786               9009               7853
# nFeature_Spatial: the number of unique genes in each sample
sum(Colon_data$nFeature_Spatial ==  colSums(Colon_data@assays$Spatial@counts > 0))
## [1] 3289
# nCount_Spatial: the total number of detected molecules in each sample
sum(Colon_data$nCount_Spatial ==  colSums(Colon_data@assays$Spatial@counts))
## [1] 3289
#2.2.4 objects(e.g. Assay) together with feature-level metadata
# A vector of names of associated objects can be had with the names function
# These can be passed to the double [[ extract operator to pull them from the Seurat object
names(x = Colon_data)
Colon_data[['Spatial']] # equivalent to: Colon_data@assays$Spatial
Colon_data[['slice1']] # equivalent to: Colon_data@images$slice1
Colon_data@assays$Spatial@counts[5:10, 1:3]
Colon_data[['Spatial']]@meta.features
head(Colon_data[['Spatial']][[]])

Colon_data[["percent.mt"]] <- PercentageFeatureSet(Colon_data,pattern = "^MT-")

VlnPlot(
  Colon_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot1 <- FeatureScatter(
  Colon_data, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(
  Colon_data, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
  NoLegend()
plot1 + plot2

Colon_subset <- subset(
  Colon_data, 
  subset = nFeature_Spatial < 8000 & nFeature_Spatial > 10 & 
    nCount_Spatial < 50000 & percent.mt < 40)

print(paste("Filter out", ncol(Colon_data) - ncol(Colon_subset), 
            "samples because of the outlier QC metrics, with", ncol(Colon_subset),
            "samples left."))

SpatialFeaturePlot(
  Colon_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "bottom")  

Colon_norm <- SCTransform(Colon_subset, assay = "Spatial", verbose = FALSE)
names(Colon_norm)
Colon_obj <- RunPCA(Colon_norm, assay = "SCT", verbose = FALSE)
# compute K nearest neighbors (KNN)
Colon_obj <- FindNeighbors(Colon_obj, reduction = "pca", dims = 1:30)
# Leiden algorithm for community detection
Colon_obj <- FindClusters(Colon_obj, verbose = FALSE)
# PCA result is the default UMAP input, use dimensions 1:30 as input features
Colon_obj <- RunUMAP(Colon_obj, reduction = "pca", dims = 1:30)
Colon_obj <- RunTSNE(Colon_obj, reduction = "pca", dims = 1:30)

plot3 <- DimPlot(Colon_obj, reduction = "umap", label = TRUE) + NoLegend()
plot4 <- SpatialDimPlot(Colon_obj, label = TRUE, label.size = 3) + NoLegend()
plot3 + plot4

Colon_obj@reductions
# identity class of each sample
table(Colon_obj@active.ident)
# find all markers of cluster 10
cluster10_markers <- FindMarkers(Colon_obj, ident.1 = 1, min.pct = 0.25)
head(cluster10_markers, n = 5)

VlnPlot(Colon_obj, features = c("EFEMP1"))
SpatialFeaturePlot(object = Colon_obj, 
                   features = rownames(cluster10_markers)[1:3], 
                   alpha = c(0.1, 1), ncol = 3)
SpatialFeaturePlot(object = Colon_obj, 
                   features = c("TIMP1","S100A10","STC1","VCX3A","CD8A"), 
                   alpha = c(0.1, 1), ncol = 3)

################################################
DimPlot(Colon_obj, reduction = "tsne", cols = pal, pt.size = 2, label = TRUE)
# 切片图像上映射分群

#paletteer_d("ggsci::default_nejm")
#paletteer_d("ggsci::lanonc_lancet")
##paletteer_d("ggsci::default_jama")
#paletteer_d("ggsci::default_jco")
#paletteer_d("ggsci::nrc_npg")
#paletteer_d("ggsci::category20_d3")
#paletteer_d("ggsci::category20b_d3")
pal <- c(paletteer_d("ggsci::category20_d3")[c(1:7,9:14,16,17,19,20)],paletteer_d("ggsci::category20b_d3"))
col <- pal[1:13]
names(col) <- levels(Colon_obj$seurat_clusters)
SpatialDimPlot(Colon_obj, label = TRUE, label.size =3,image.alpha = 0.8,crop = F) + scale_fill_manual(values = col)
p1 <- SpatialDimPlot(Colon_obj, label = TRUE, label.size = 3) + scale_fill_manual(values = col)




load("~/Projects/ImmunePro/human.immune.CIBERSORT.RData")
human.immune.CIBERSORT
human.immune.CIBERSORT<-split(human.immune.CIBERSORT$V2,human.immune.CIBERSORT$V1)
names(human.immune.CIBERSORT)
human.immune.CIBERSORT<-human.immune.CIBERSORT[c(10,16)]
names(human.immune.CIBERSORT)
Colon_obj<-AddModuleScore(Colon_obj,features =human.immune.CIBERSORT,col.name =names(human.immune.CIBERSORT))
names(human.immune.CIBERSORT)

colnames( Colon_obj@meta.data)
SpatialFeaturePlot(object = Colon_obj, 
                   features = c("S100A10" ), 
                   #alpha = c(0.1, 1),
                   ncol = 2)

SpatialFeaturePlot(object = Colon_obj, 
                   features = c("percent.mt"),image.alpha = 0.8,crop = F,#alpha = c(0.1, 6),
                   ncol = 1) +  ggplot2::scale_fill_gradient2(limits = c(0.0, 10), breaks = c(0.0, 6, 10),  low ="black", mid = "white", high = "red", midpoint =3)

#ggsave(file.path(Figure3,"spatital_fibroblast_Colonsamples_orignial.pdf"),width=4,height =5)

p1<-SpatialFeaturePlot(object = Colon_obj, 
                       features = c("Fibroblast_C0"), 
                       #cols=colorRampPalette(brewer.pal(n = 4, name ="YlOrRd"))(10),
                       #alpha = c(0.1, 6),
                       ncol = 1) +  ggplot2::scale_fill_gradient2(limits = c(0.0, 15), breaks = c(0.0,6, 15),   low ="black", mid = "white", high = "red", midpoint =6)
p2<-SpatialFeaturePlot(object = Colon_obj, 
                       features = c("Fibroblast_C1"), 
                       #cols=colorRampPalette(brewer.pal(n = 4, name ="YlOrRd"))(10),
                       #alpha = c(0.1, 6),
                       ncol = 1) +  ggplot2::scale_fill_gradient2(limits = c(0.0, 15), breaks = c(0.0, 6, 15), low ="black", mid = "white", high = "red", midpoint = 6)
p3<-SpatialFeaturePlot(object = Colon_obj, 
                       features = c("Fibroblast_C2"   
                       ), 
                       #cols=colorRampPalette(brewer.pal(n = 4, name ="YlOrRd"))(10),
                       #alpha = c(0.1, 6),
                       ncol = 1) +  ggplot2::scale_fill_gradient2(limits = c(0.0, 15), breaks = c(0.0,8, 15),low = "black", mid = "white", high = "red", midpoint =6)
p4<-SpatialFeaturePlot(object = Colon_obj, 
                       features = c("Fibroblast_C3"), 
                       #cols=colorRampPalette(brewer.pal(n = 4, name ="YlOrRd"))(10),
                       #alpha = c(0.1, 6),
                       ncol = 1) +  ggplot2::scale_fill_gradient2(limits = c(0.0, 15), breaks = c(0.0,6, 15), low ="black", mid = "white", high = "red", midpoint =6)
#ggplot2::scale_fill_continuous(limits = c(0.0,10), breaks = c(0.0, 5, 10))
p1
p1/p2/p3/p4
#ggsave(file.path(Figure3,"spatital_fibroblast_Colonsamples.pdf"),width=4,height =15)
p1+p2+p3+p4
#scale_color_gradient2(low="#1B9E77FF",mid="white",high="#D51113" )


###
Colon_obj_markers <- FindAllMarkers(Colon_obj, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
Colon_obj_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

top10 <- Colon_obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Colon_obj, features = top10$gene) + NoLegend()




