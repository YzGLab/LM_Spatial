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
  data.dir ='/home/data/gaoyuzhen/Projects/LM_spatial/spatial/C_1',
  filename = "filtered_feature_bc_matrix.h5",
  #assay = "Spatial", # specify name of the initial assay
  #slice = "spatial", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)
dim(x = Colon_data) # the number of features (genes) by samples (spots)
nrow(x = Colon_data) # the number of features
ncol(x = Colon_data) # the number of samples
#Extract the feature and sample names using rownames or colnames.
head(x = rownames(Colon_data), n = 5)
## [1] "Xkr4"    "Gm1992"  "Gm19938" "Gm37381" "Rp1"
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
?SpatialDimPlot



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
ggsave(file.path(Figure3,"spatital_fibroblast_Colonsamples.pdf"),width=4,height =15)

p1+p2+p3+p4





p1+p2+p4
#scale_color_gradient2(low="#1B9E77FF",mid="white",high="#D51113" )

?SpatialFeaturePlot


####
Colon_obj_markers <- FindAllMarkers(Colon_obj, only.pos = TRUE, min.pct = 0.25, 
                                   logfc.threshold = 0.25)
Colon_obj_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

top10 <- Colon_obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Colon_obj, features = top10$gene) + NoLegend()




#remotes::install_github("RubD/Giotto") 
###/home/data/gaoyuzhen/Projects/LiverCancer/data/HCC-scRNA/ST data/Colon/HCC-2P
#install.packages("magick")
library(Giotto)
#results_folder = '/path/to/directory/'
results_folder = '/home/data/gaoyuzhen/Projects/LiverCancer/Results_LC//'
# set python path to NULL if you want to automatically install (only the 1st time) and use the giotto miniconda environment
python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}

############
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE)

## provide path to visium folder
data_path = '/home/data/gaoyuzhen/Projects/ImmunePro/Data/CancerDiscovery_livermeta/ST/ST_liver1/'

## directly from visium folder


visium_LM<-createGiottoVisiumObject( visium_dir ='/home/data/gaoyuzhen/Projects/ImmunePro/Data/CancerDiscovery_livermeta/ST/ST_liver1',
                                     expr_data =  "filter",
                                     #gene_column_index = 1,
                                     h5_visium_path ="/home/data/gaoyuzhen/Projects/ImmunePro/Data/CancerDiscovery_livermeta/ST/ST_liver1/filtered_feature_bc_matrix.h5",
                                     h5_gene_ids = "symbols",
                                     h5_tissue_positions_path = "/home/data/gaoyuzhen/Projects/ImmunePro/Data/CancerDiscovery_livermeta/ST/ST_liver1/tissue_positions_list.csv",
                                     h5_image_png_path ='/home/data/gaoyuzhen/Projects/ImmunePro/Data/CancerDiscovery_livermeta/ST/ST_liver1/tissue_lowres_image.png',
                                     png_name = "tissue_lowres_image.png",
                                     #do_manual_adj = TRUE,
                                     xmax_adj = 0,
                                     xmin_adj = 0,
                                     ymax_adj = 0,
                                     ymin_adj = 0,
                                     #scale_factor = NULL,
                                     instructions = instrs
                                     #cores = NA,
                                     # verbose = TRUE
)






## update and align background image
# problem: image is not perfectly aligned
spatPlot(gobject = visium_LM, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
         save_param = list(save_name = '2_a_spatplot_image'))
##
# check name
showGiottoImageNames(visium_LM) # "image" is the default name
# adjust parameters to align image (iterative approach)
visium_LM = updateGiottoImage(visium_LM, image_name = 'image',
                              xmax_adj = 400, xmin_adj = 600,
                              ymax_adj = 2100, ymin_adj = 2400)

# now it's aligned
spatPlot(gobject = visium_LM, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
         save_param = list(save_name = '2_b_spatplot_image_adjusted'))
##
## check metadata
pDataDT(visium_LM)

## compare in tissue with provided jpg
spatPlot(gobject = visium_LM, cell_color = 'in_tissue', point_size = 2,
         cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
         save_param = list(save_name = '2_c_in_tissue'))
##
## subset on spots that were covered by tissue
metadata = pDataDT(visium_LM)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_LM = subsetGiotto(visium_LM, cell_ids = in_tissue_barcodes)

## filter
visium_LM <- filterGiotto(gobject = visium_LM,
                          expression_threshold = 1,
                          gene_det_in_min_cells = 10,
                          min_det_genes_per_cell = 500,
                          expression_values = c('raw'),
                          verbose = T)

## normalize
visium_LM <- normalizeGiotto(gobject = visium_LM, scalefactor = 10000, verbose = T)

## add gene & cell statistics
visium_LM <- addStatistics(gobject = visium_LM)
###
"GPR171" %in% visium_LM@gene_ID


## visualize
patPlot2D(gobject = visium_LM, show_image = T, point_alpha = 0.7,
          ave_param = list(save_name = '2_d_spatial_locations'))
spatPlot2D(gobject = visium_LM, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_genes', color_as_factor = F,
           save_param = list(save_name = '2_e_nr_genes'))
## highly variable genes (HVG)
visium_LM <- calculateHVG(gobject = visium_LM,
                          save_param = list(save_name = '3_a_HVGplot'))
#
## run PCA on expression values (default)
gene_metadata = fDataDT(visium_LM)
featgenes = gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID

visium_LM <- runPCA(gobject = visium_LM,
                    genes_to_use = featgenes,
                    scale_unit = F, center = T,
                    method="factominer")

screePlot(visium_LM, ncp = 30, save_param = list(save_name = '3_b_screeplot'))
plotPCA(gobject = visium_LM,
        save_param = list(save_name = '3_c_PCA_reduction'))
## run UMAP and tSNE on PCA space (default)
visium_LM <- runUMAP(visium_LM, dimensions_to_use = 1:10)
plotUMAP(gobject = visium_LM,
         save_param = list(save_name = '3_d_UMAP_reduction'))
visium_LM <- runtSNE(visium_LM, dimensions_to_use = 1:10)
plotTSNE(gobject = visium_LM,
         save_param = list(save_name = '3_e_tSNE_reduction'))
#
## sNN network (default)
visium_LM <- createNearestNetwork(gobject = visium_LM, dimensions_to_use = 1:10, k = 15)
## Leiden clustering
visium_LM <- doLeidenCluster(gobject = visium_LM, resolution = 0.2, n_iterations = 1000)
plotUMAP(gobject = visium_LM,
         cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5,
         save_param = list(save_name = '4_a_UMAP_leiden'))
#
# expression and spatial
spatDimPlot(gobject = visium_LM, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5,
            save_param = list(save_name = '5_a_covis_leiden'))
spatDimPlot(gobject = visium_LM, cell_color = 'nr_genes', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5,
            save_param = list(save_name = '5_b_nr_genes'))




###
DG_subset = subsetGiottoLocs(visium_LM,
                             x_max = 6500, x_min = 3000,
                             y_max = -2500, y_min = -5500,
                             return_gobject = TRUE)

spatDimPlot(gobject = DG_subset,
            cell_color = 'leiden_clus', spat_point_size = 5,
            save_param = list(save_name = '5_c_DEG_subset'))
##
#6. Cell-Type Marker Gene Detection
#6.1 Gini
gini_markers_subclusters = findMarkers_one_vs_all(gobject = visium_LM,
                                                  method = 'gini',
                                                  expression_values = 'normalized',
                                                  cluster_column = 'leiden_clus',
                                                  min_genes = 20,
                                                  min_expr_gini_score = 0.5,
                                                  min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes

# violinplot
violinPlot(visium_LM, genes = unique(topgenes_gini), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = list(save_name = '6_a_violinplot_gini', base_width = 5, base_height = 10))
##
violinPlot(visium_LM, genes = c("SLC22A2" , "CST1"), cluster_column = 'leiden_clus',
           strip_text = 8, strip_position = 'right',
           save_param = list(save_name = '6_a_violinplot_gini', base_width = 5, base_height = 10))

# cluster heatmap
plotMetaDataHeatmap(visium_LM, selected_genes = topgenes_gini,
                    metadata_cols = c('leiden_clus'),
                    x_text_size = 10, y_text_size = 10,
                    save_param = list(save_name = '6_b_metaheatmap_gini'))
##
# umap plots
dimGenePlot2D(visium_LM, expression_values = 'scaled',
              genes = gini_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes,
              cow_n_col = 3, point_size = 1,
              save_param = list(save_name = '6_c_gini_umap', base_width = 8, base_height = 5))


dimGenePlot2D(visium_LM, expression_values = 'scaled',
              genes = c("CD8A","GPR171"),
              cow_n_col = 2, point_size = 1,
              save_param = list(save_name = '6_c_gini_umap_GPR171', base_width = 8, base_height = 5))


##
BiocManager::install('scran')
#6.2 Scran
scran_markers_subclusters = findMarkers_one_vs_all(gobject = visium_LM,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes

# violinplot
violinPlot(visium_LM, genes = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(save_name = '6_d_violinplot_scran', base_width = 5))

# cluster heatmap
plotMetaDataHeatmap(visium_LM, selected_genes = topgenes_scran,
                    metadata_cols = c('leiden_clus'))

# umap plots
dimGenePlot(visium_LM, expression_values = 'scaled',
            genes = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes,
            cow_n_col = 3, point_size = 1,
            save_param = list(save_name = '6_f_scran_umap', base_width = 8, base_height = 5))

dimGenePlot(visium_LM, expression_values = 'scaled',
            genes =c("CD8A","GPR171"),
            cow_n_col = 3, point_size = 1,
            save_param = list(save_name = '6_f_scran_umap_GPR171', base_width = 8, base_height = 5))


#####
sign_matrix_path = system.file("extdata", "sig_matrix.txt", package = 'Giotto')
# cell-type annotation


##
load("~/Projects/ImmunePro/human.immune.CIBERSORT.RData")
human.immune.CIBERSORT<-split(human.immune.CIBERSORT$V2,human.immune.CIBERSORT$V1)

#"CD8Tcells" "CD8Tex" 


human.immune.CIBERSORT<-split(sce.markers_GPR$gene,sce.markers_GPR$cluster)

names(human.immune.CIBERSORT)<-c("low","high")

human.immune.CIBERSORT<-human.immune.CIBERSORT[c(6,7)]

sig_matrix = makeSignMatrixPAGE(sign_names = names(human.immune.CIBERSORT),
                                sign_list =human.immune.CIBERSORT)
names(human.immune.CIBERSORT)
# 1.3 enrichment test with PAGE

# runSpatialEnrich() can also be used as a wrapper for all currently provided enrichment options
visium_LM = runPAGEEnrich(gobject = visium_LM, sign_matrix = sig_matrix)
scorelist<-visium_LM@spatial_enrichment$PAGE
# 1.4 heatmap of enrichment versus annotation (e.g. clustering result)
cell_types = colnames(sig_matrix)
plotMetaDataCellsHeatmap(gobject = visium_LM,
                         metadata_cols = 'leiden_clus',
                         value_cols = cell_types,
                         spat_enr_names = 'PAGE',
                         x_text_size = 8,
                         y_text_size = 8,
                         save_param = list(save_name="7_a_metaheatmap"))
# 1.5 visualizations
cell_types_subset = colnames(sig_matrix)[1:10]
spatCellPlot(gobject = visium_LM,
             spat_enr_names = 'PAGE',
             cell_annotation_values = cell_types_subset,
             cow_n_col = 4,coord_fix_ratio = NULL, point_size = 0.75,
             save_param = list(save_name="7_b_spatcellplot_1"))
?spatCellPlot

spatCellPlot(gobject = visium_LM,
             spat_enr_names = 'PAGE',
             cell_annotation_values = cell_types_subset,
             cell_color_gradient(c("blue", "white", "red")),
             cow_n_col = 2,coord_fix_ratio = NULL, point_size = 1.75,
             save_param = list(save_name="7_b_spatcellplot_1")) 
#scale_fill_gradientn(colors = c("blue", "white", "red"), na.value = "white")

cell_types_subset = colnames(sig_matrix)[11:20]

spatCellPlot(gobject = visium_LM, spat_enr_names = 'PAGE',
             cell_annotation_values = cell_types_subset, cow_n_col = 4,
             coord_fix_ratio = 1, point_size = 0.75,
             save_param = list(save_name="7_c_spatcellplot_2"))

spatDimCellPlot(gobject = visium_LM,
                spat_enr_names = 'PAGE',
                cell_annotation_values = c( "CD8Tcells","CD8Tex","Treg"),
                cow_n_col = 1, spat_point_size = 1,
                plot_alignment = 'horizontal',
                save_param = list(save_name="7_d_spatDimCellPlot", base_width=7, base_height=10))
##8. Spatial Grid
visium_LM <- createSpatialGrid(gobject = visium_LM,
                               sdimx_stepsize = 400,
                               sdimy_stepsize = 400,
                               minimum_padding = 0)
spatPlot(visium_LM, cell_color = 'leiden_clus', show_grid = T,
         grid_color = 'red', spatial_grid_name = 'spatial_grid',
         save_param = list(save_name = '8_grid'))
#
visium_LM <- createSpatialNetwork(gobject = visium_LM,
                                  method = 'kNN', k = 5,
                                  maximum_distance_knn = 400,
                                  name = 'spatial_network')

showNetworks(visium_LM)

spatPlot(gobject = visium_LM, show_network = T,
         network_color = 'blue', spatial_network_name = 'spatial_network',
         save_param = list(save_name = '9_a_knn_network'))
#10. Spatial Genes and Patterns
#10.1 Spatial Genes
## kmeans binarization
kmtest = binSpect(visium_LM, calc_hub = T, hub_min_int = 5,
                  spatial_network_name = 'spatial_network')
spatGenePlot(visium_LM, expression_values = 'scaled',
             genes = kmtest$genes[1:6], cow_n_col = 2, point_size = 1.5,
             save_param = list(save_name = '10_a_spatial_genes_km'))
spatGenePlot(visium_LM, expression_values = 'scaled', point_border_stroke =0.1,gradient_midpoint=0.5,
             genes = "GPR171", cow_n_col = 1, point_size = 1,
             save_param = list(save_name = '10_a_spatial_genes_km'))

