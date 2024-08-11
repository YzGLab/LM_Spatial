##fib LM vs CRC
###diff of lM vs CRC for st _celltypes
library(SCP)
library(Seurat)
###############
load("data/stmetadata.Rdata")
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
#######################################################
sub_sce<-readRDS("data/spatial_combinend_removed_MT_RP.RDS")
####################

Idents(sub_sce)<-sub_sce$st_celltype
sce.markers1 <- FindAllMarkers(sub_sce[,sub_sce$set=="CRC"], 
                               only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
sce.markers2 <- FindAllMarkers(sub_sce[,sub_sce$set=="LM"], 
                               only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
#sce.markers.sig <- sce.markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10

sce.markers1$cluster<-paste0(sce.markers1$cluster,".CRC")
sce.markers2$cluster<-paste0(sce.markers2$cluster,".LM")
sce.markers<-rbind(sce.markers1,sce.markers2)################################
sce.cell_type_markers<-split(sce.markers$gene,sce.markers$cluster)###############合并两种差异
save(sce.cell_type_markers,file='data/spatial_differetgenes_clusters.Rdata')


###
table(sub_sce$st_celltype)
fib_sub_sce<-sub_sce[,sub_sce$st_celltype=="Fibroblast"]

Idents(fib_sub_sce)<-fib_sub_sce$set
sce.markers.fib<- FindAllMarkers(fib_sub_sce, 
                               only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
################
sce.markers.fib.sig <- sce.markers.fib %>% group_by(cluster) %>% top_n(n =100, wt = avg_log2FC)# top 10

######################
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
                       pvalueCutoff  = 0.01,qvalueCutoff  = 0.01, readable = TRUE)
  }
  return(enrich)
}
###########################
sce.markers_celltype_LM_vs_CRC <- split(sce.markers.fib.sig$gene, sce.markers.fib.sig$cluster)
ProSig <- sce.markers_celltype_LM_vs_CRC####设置基因集

Sigbp <- lapply(ProSig, Myenrich, category = "gobp", geneid = "SYMBOL")
Sigkegg <- lapply(ProSig, Myenrich, category = "kegg", geneid = "SYMBOL")
#
pdf(file=file.path(Figure5,"Cluster_Go-up.pdf"),width = 6,height = 8)
dotplot(Sigbp[[1]],showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free') + scale_color_continuous(low='purple', high='green')+ aes(shape=GeneRatio > 0.1)
dev.off()

pdf(file=file.path(Figure5,"stRNA_CRC_fib_Cluster_Go-net.pdf"),width = 9,height = 8)
Sigbp1<- pairwise_termsim(Sigbp[[1]])
emapplot(Sigbp1, showCategory = 10,color = "p.adjust")+ scale_color_continuous(low='purple', high='green')
dev.off()

pdf(file=file.path(Figure5,"stRNA_LM_fib_Cluster_Go-net.pdf"),width = 9,height = 8)
Sigbp2<- pairwise_termsim(Sigbp[[2]])
emapplot(Sigbp2, showCategory = 10,color = "p.adjust")+ scale_color_continuous(low='purple', high='green')
dev.off()



pdf(file=file.path(Figure5,"Cluster_Go-circ.pdf"),width = 15,height = 8)
cnetplot(Sigbp2, showCategory = 15,categorySize="count",circular = FALSE, colorEdge = TRUE,cex_gene = 0.4,cex_label_gene =0.5)
dev.off()



#### gsoap first show the pathway 
###
names(Sigbp)
##展示所有 富集信息
x<-Sigkegg[[2]]
# Convert to data frame and set BP description as rownames
x = as.data.frame(x, row.names = x$Description)

l = gsoap_layout(x,
                 genes = 'geneID',
                 pvalues = 'p.adjust',
                 projection = 'tsne',
                 scale.factor = 0.8,
                 no.clusters = 4)

# idx = which(l$cluster == 'Cluster 1')
idx1 = l %>% 
  tibble::rownames_to_column("feature") %>% 
  group_by(cluster) %>% 
  slice_max(order_by = significance, n = 3)
idx = which(rownames(l) %in% idx1$feature)

library(ggrepel)
library(ggforce)
p = gsoap_plot(l,
               as.alpha = 'significance',
               as.color = 'cluster',
               which.labels = idx,
               viridis.option = 'plasma', 
               # viridis.option = 'mako',  ### mac
               #viridis.option = 'cividis',
               viridis.direction = 1,
               viridis.range = c(.2, .8),
               size.guide.loc = c(1., 1.),
               label.fontsize = 10)
p
ggsave(file=file.path(Figure5,"gsoap_plot_LM_fib_stRNA.pdf"),width = 9,height = 9)
###CRC
##展示所有 富集信息
x<-Sigbp[[1]]
# Convert to data frame and set BP description as rownames
x = as.data.frame(x, row.names = x$Description)

l = gsoap_layout(x,
                 genes = 'geneID',
                 pvalues = 'p.adjust',
                 projection = 'tsne',
                 scale.factor = 0.8,
                 no.clusters = 1)

# idx = which(l$cluster == 'Cluster 1')
idx1 = l %>% 
  tibble::rownames_to_column("feature") %>% 
  group_by(cluster) %>% 
  slice_max(order_by = significance, n = 3)
idx = which(rownames(l) %in% idx1$feature)

library(ggrepel)
library(ggforce)
p = gsoap_plot(l,
               as.alpha = 'significance',
               as.color = 'cluster',
               which.labels = idx,
               viridis.option = 'plasma', 
               # viridis.option = 'mako',  ### mac
               #viridis.option = 'cividis',
               viridis.direction = 1,
               viridis.range = c(.2, .8),
               size.guide.loc = c(1., 1.),
               label.fontsize = 10)
p
ggsave(file=file.path(Figure5,"gsoap_plot_LM_fib_stRNA.pdf"),width = 9,height = 9)
