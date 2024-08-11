
####cell pathway
library(data.table)
#library(GEOquery)
library(dplyr)
#library(Seurat)
library(ggplot2)
library(umap)
library(hypeR)
library(msigdbr)
library(tidyverse)
library(pheatmap)
library(clustree)
library(survival)
library(survminer)
library(GSVA)
library(clusterProfiler)
library(enrichplot)
##load data
sub_sce<-readRDS("data/spatial_combinend_removed_MT_RP.RDS")
####################

Idents(sub_sce)<-sub_sce$st_celltype
sce.markers1 <- FindAllMarkers(sub_sce[,sub_sce$set=="CRC"], 
                              only.pos = TRUE, 
                              min.pct = 0.55, logfc.threshold = 0.55)
sce.markers2 <- FindAllMarkers(sub_sce[,sub_sce$set=="LM"], 
                               only.pos = TRUE, 
                               min.pct = 0.55, logfc.threshold = 0.55)
#sce.markers.sig <- sce.markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10

sce.markers1$cluster<-paste0("CRC_",sce.markers1$cluster)
sce.markers2$cluster<-paste0("LM_",sce.markers2$cluster)
sce.markers<-rbind(sce.markers1,sce.markers2)################################
sce.cell_type_markers<-split(sce.markers$gene,sce.markers$cluster)###############合并两种差异

save(sce.cell_type_markers,file='data/spatial_differetgenes_clusters.Rdata')
##pathway analysis
Myenrich <- function(genes, category = c("kegg", "gobp"), 
                     geneid = c("SYMBOL", "ENTREZID", "ENSEMBL", "UNIPROT")){
  library(clusterProfiler)
  category <- match.arg(category)
  geneid <- match.arg(geneid)
  if (geneid != "ENTREZID"){
    genes <- bitr(genes, fromType = geneid,toType ="ENTREZID",OrgDb="org.Hs.eg.db") %>% .[, 2] %>% as.character()
  }else{genes <- genes}
  if (category == "kegg"){
    enrich <- enrichKEGG(genes, organism = "hsa",keyType = "kegg",
                         pAdjustMethod = "BH",pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.01, maxGSSize = 5000)
  }
  if (category == "gobp"){
    enrich <- enrichGO(genes, OrgDb="org.Hs.eg.db",ont= "BP",pAdjustMethod = "BH",
                       pvalueCutoff  = 0.0001,qvalueCutoff  = 0.0001, readable = TRUE)
  }
  return(enrich)
}

###########################
ProSig <- sce.cell_type_markers####设置基因集
Sigbp <- lapply(ProSig, Myenrich, category = "gobp", geneid = "SYMBOL")
Sigkegg <- lapply(ProSig, Myenrich, category = "kegg", geneid = "SYMBOL")

data.frame(Sigbp[1])
??clusterProfiler

edox1 <- pairwise_termsim(Sigbp[[6]])
edox2 <- pairwise_termsim(Sigbp[[12]])
p1 <- treeplot(edox1, hclust_method = "average")
p2 <- treeplot(edox2, hclust_method = "average")
p1+p2
ggsave(file.path(Figure2,"CRC_Tumor_vs_LM_Tumor_pathway_GO.pdf"),
       width = 11.5, height =4.73,
       bg = "transparent", dpi = 300, useDingbats=FALSE)

p1<-emapplot(edox1, hclust_method = "average")
p2<-emapplot(edox2, hclust_method = "average")
p1+p2
ggsave(file.path(Figure2,"CRC_Tumor_vs_LM_Tumor_pathway_GO_network.pdf"),
       width = 12.5, height =8.73,
       bg = "transparent", dpi = 300, useDingbats=FALSE)
emapplot(edox1, layout="kk")
emapplot(edox2, cex_category=1)



############################plot enrichmenet map for Go 未运行
names(Sigbp)<-gsub('/',"_",names(Sigbp))
singele_bp_enrichmap<-function(x){
  print(x)
  edox2 <- pairwise_termsim(Sigbp[[x]])
  p1<-emapplot(edox2, cex_category=1.5)
  ggsave(file.path(Figure3,paste0(x,"GO_enrichmap.pdf")),width =7,height=7.33)
}
sapply(names(Sigbp)[-c(5,10)],singele_bp_enrichmap)
###########################


#################################################plot enrichmenet map for KEGG 
names(Sigkegg)<-gsub('/',"_",names(Sigkegg))
singele_enrichmap<-function(x){
  print(x)
  edox2 <- pairwise_termsim(Sigkegg[[x]])
  p1<-emapplot(edox2, cex_category=1.5)
  ggsave(file.path(Figure3,paste0(x,"KEGG_enrichmap.pdf")),width =7,height=7.33)
}
sapply(names(Sigkegg)[-c(5,10)],singele_enrichmap)




# Load necessary libraries
library(clusterProfiler)
library(ggplot2)
library(DOSE)
library(enrichplot)
# Save Sigbp and Sigkegg results
save(Sigbp, file = 'data/Sigbp_results.RData')
save(Sigkegg, file = 'data/Sigkegg_results.RData')

# Function to create dot plots for top 20 GO terms
plot_top20_GO_dotplot <- function(enrich_result, title) {
  if (is.null(enrich_result) || length(enrich_result) == 0) {
    print(paste("No enrichment results for", title))
    return(NULL)
  }
  
  dotplot(enrich_result, x = "GeneRatio", color = "p.adjust", showCategory = 20, title = title) + 
    theme_minimal() +
    labs(title = title)
}

# Function to save enrichment results to CSV
save_enrichment_to_csv <- function(enrich_result, filename) {
  if (!is.null(enrich_result) && nrow(as.data.frame(enrich_result)) > 0) {
    write.csv(as.data.frame(enrich_result), file = filename, row.names = FALSE)
  } else {
    print(paste("No enrichment results to save for", filename))
  }
}
# Plot and save the top 20 GO terms for each cell type in Sigbp



for (name in names(Sigbp)[-c(5,10)]) {
  pdf(file.path(Figure3,paste0(name,"Sigbp_top20_dotplots.pdf")))
  plot_result <- plot_top20_GO_dotplot(Sigbp[[name]], paste("Top 20 GO BP for", name))
  if (!is.null(plot_result)) {
    print(plot_result)
  }
  dev.off()
  # Save enrichment results to CSV
  save_enrichment_to_csv(Sigbp[[name]], file.path(Figure3,paste0(name, "_results.csv")))
}


for (name in names(Sigkegg)[-c(5,10,11)]) {
  pdf(file.path(Figure3,paste0(name,"Sigkegg_top20_dotplots.pdf")))
  plot_result <- plot_top20_GO_dotplot(Sigkegg[[name]], paste("Top 20 Sigkegg for", name))
  if (!is.null(plot_result)) {
    print(plot_result)
  }
  dev.off()
  # Save enrichment results to CSV
  save_enrichment_to_csv(Sigkegg[[name]], file.path(Figure3,paste0(name, "Sigkegg_results.csv")))
}

