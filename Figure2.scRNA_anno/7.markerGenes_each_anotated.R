####
library(Seurat)
sce<-readRDS("data/Colon_HC_combined.sce.rds")
load("data/totalscRNAdata.Rdata")
#####
table(sce@meta.data$sub_sub_celltype)
metadata<-sce@meta.data
## 根据差异基因进行标记~ 可以后面重新标记

############################################################################################################
metadata$sub_celltype_genes<-NULL
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-B"]<-"B(C0)-FBXO32"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-B"]<-"B(C1)-S100A4"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C2-B"]<-"B(C2)-IGHM"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C3-B"]<-"B(C3)-FCRL1"
############################################################################################################
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-Plasma"]<-"Plasma(C0)-CRADD"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-Plasma"]<-"Plasma(C1)-ASTL"
############################################################################################################
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-CD4"]<-"CD4(C0)-SEMA4A"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-CD4"]<-"CD4(C1)-IKZF2"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C2-CD4"]<-"CD4(C2)-S100A4"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C3-CD4"]<-"CD4(C3)-S100A6"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C4-CD4"]<-"CD4(C4)-CCL5"
############################################################################################################
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-Treg"]<-"Treg(C0)-SLC16A10"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-Treg"]<-"Treg(C1)-LTB"
############################################################################################################
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-CD8"]<-"CD8(C0)-GZMA"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-CD8"]<-"CD8(C1)-CCL3L1"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C2-CD8"]<-"CD8(C2)-APOLD1"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C3-CD8"]<-"CD8(C3)-CCL4L2"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C4-CD8"]<-"CD8(C4)-IMMP2L"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C5-CD8"]<-"CD8(C5)-GNLY"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C6-CD8"]<-"CD8(C6)-IER3"
############################################################################################################

metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-NK"]<-"NK(C0)-HDAC9"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-NK"]<-"NK(C1)-GZMB"
############################################################################################################
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-MAIT"]<-"MAIT(C0)-MEGF9"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-MAIT"]<-"MAIT(C1)-KLRB1"
############################################################################################################
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-Myeloid"]<-"Myeloid(C0)-S100A9"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-Myeloid"]<-"Myeloid(C1)-HES1"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C2-Myeloid"]<-"Myeloid(C2)-CD74"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C3-Myeloid"]<-"Myeloid(C3)-AFF3"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C4-Myeloid"]<-"Myeloid(C4)-TNF"
############################################################################################################
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-Endothelial"]<-"End(C0)-SLC9A3R2"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-Endothelial"]<-"End(C1)-ROBO2"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C2-Endothelial"]<-"End(C2)-CCL14"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C3-Endothelial"]<-"End(C3)-FLRT2"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C4-Endothelial"]<-"End(C4)-IFI27"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C5-Endothelial"]<-"End(C5)-UGDH"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C6-Endothelial"]<-"End(C6)-HMGB2"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C7-Endothelial"]<-"End(C7)-COL1A1"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C8-Endothelial"]<-"End(C8)-PLAC9"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C9-Endothelial"]<-"End(C9)-TFPI"
############################################################################################################

metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-Epithelial"]<-"Epi(C0)-SHISA6"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-Epithelial"]<-"Epi(C1)-FABP1"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C2-Epithelial"]<-"Epi(C2)-ZNF331"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C3-Epithelial"]<-"Epi(C3)-FTH1"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C4-Epithelial"]<-"Epi(C4)-CACNB2"

############################################################################################################
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C0-Fibroblasts"]<-"Fib(C0)-CCDC102B"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C1-Fibroblasts"]<-"Fib(C1)-ADIRF"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C2-Fibroblasts"]<-"Fib(C2)-PLA2G5"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C3-Fibroblasts"]<-"Fib(C3)-DCN"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C4-Fibroblasts"]<-"Fib(C4)-HMCN1"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C5-Fibroblasts"]<-"Fib(C5)-ONECUT1"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C6-Fibroblasts"]<-"Fib(C6)-TNF"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C7-Fibroblasts"]<-"Fib(C7)-MCTP1"
metadata$sub_celltype_genes[metadata$sub_sub_celltype=="C8-Fibroblasts"]<-"Fib(C7)-TOP2A"

# Extract the main cell type (e.g., CD4, CD8) from the metadata$sub_celltype_genes
main_cell_types <- sapply(strsplit(as.character(metadata$sub_celltype_genes), "-"), `[`, 1)
# Define the order of the main cell types
ordered_cell_types <- c("B", "Plasma", "CD4", "Treg", "CD8", "NK", "MAIT", "Myeloid", "End", "Epi", "Fib")
# Create a factor with the specified order
main_cell_types_factor <- factor(main_cell_types, levels = ordered_cell_types)
# Sort metadata$sub_celltype_genes by this factor
metadata <- metadata[order(main_cell_types_factor), ]
################################################################

sce<-AddMetaData(sce,metadata)
colnames(sce@meta.data)
CellDimPlot(
  srt =sce, group.by =  "seurat_clusters", stat_plot_size = 0.1,split.by ="immunecells",
  reduction = "UMAP", theme_use = "theme_blank",show_stat = TRUE,
  palcolor = mycol2,
)

