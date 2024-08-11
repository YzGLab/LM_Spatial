
sce<-readRDS("data/metabolits_scp_st_sp_sce.RDS")
metadata<-sce@meta.data
metabolites_obj<-sce@assays$metabolits@data %>% data.frame(.,check.names = FALSE)
# Libraries
library(Seurat)
library(dplyr)

##t淘汰的方法 单细胞方法。不适合代谢组学
setwd("/home/data/gaoyuzhen/Projects/LM_spatial")
# Function to filter and find markers for each combination
find_markers <- function(sce, sample, celltype1, celltype2) {
  sub_sce <- sce[, sce@meta.data$orig.ident == sample & sce@meta.data$st_celltype %in% c(celltype1, celltype2)]
  Idents(sub_sce) <- sub_sce$st_celltype
  markers <- FindAllMarkers(sub_sce, 
                            #only.pos = TRUE, 
                            test.use = "wilcox",
                            min.pct = 0.25, logfc.threshold = 0.25)
  return(markers)
}

# List of cell type combinations
celltype_combinations <- list(
  c("Fibroblast", "Hepatocytes"),
  #c("Fibroblast", "Tumor"),
  c("Fibroblast", "Normal Epi")
)

# Initialize list to store markers for each sample
all_samples_markers <- list()

# List of samples
samples <- c("LM_1", "LM_2")

# Loop through samples and find markers for each combination
for (sample in samples) {
  all_markers <- list()
  for (comb in celltype_combinations) {
    markers <- find_markers(sce, sample, comb[1], comb[2])
    markers_sig <- markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
    markers_sig <- markers_sig[markers_sig$cluster == "Fibroblast", ]
    all_markers[[paste(comb, collapse = "_vs_")]] <- markers_sig
  }
  
  # Extract and cross-filter genes from the markers
  common_genes <- Reduce(intersect, lapply(all_markers, function(x) unique(x$gene)))
  
  # Filter markers to only include common genes
  all_markers_filtered <- lapply(all_markers, function(x) x[x$gene %in% common_genes, ])
  
  # Combine the filtered marker lists into one data frame
  all_markers_filtered_combined <- do.call(rbind, lapply(names(all_markers_filtered), function(name) {
    df <- all_markers_filtered[[name]]
    df$comparison <- name
    return(df)
  }))
  
  # Filter for adjusted p-value < 0.001
  all_markers_filtered_combined_LM <- all_markers_filtered_combined[all_markers_filtered_combined$p_val_adj < 0.01, ]
  
  # Add sample name
  all_markers_filtered_combined_LM$sample_name <- sample
  
  # Store results for the current sample
  all_samples_markers[[sample]] <- all_markers_filtered_combined_LM
}

# Combine results from all samples
combined_LM <- do.call(rbind, all_samples_markers)

# Display the combined results
combined_LM
