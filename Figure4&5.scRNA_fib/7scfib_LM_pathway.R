### pathway of fib4 and st fib clusters 
######################
load(file.path(Figure4,"Fibroblasts_markergenes_split_markers_sig.Rdata"))
fib_sce<-readRDS("data/fib_sce_scRNA.Rds")
############################
table(fib_sce$sub_celltype_genes)
Idents(fib_sce)<-fib_sce$sub_celltype_genes
#########################################
sce.markers<- FindAllMarkers(fib_sce, 
                             only.pos = TRUE, 
                             min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.x.sig<-sce.markers[sce.markers<0.001,]
sce.markers.x.sig<- sce.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)# top 10 
write.csv(sce.markers,file.path(Figure4,"Fibroblasts_markergenes_split.csv"))
write.csv(sce.markers.x.sig,file=file.path(Figure4,"Fibroblasts_markergenes_split_sce.markers_sig.csv"))
###########################################################################################################
sce.markers.x.sig.10<- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)# top 10 
sce.markers.x.sig.10[sce.markers.x.sig.10$cluster=="Fib(C4)-CXCL8",]$gene
######可视化 差异分析结果

sub_sce.x <- RunDEtest(srt =fib_sce[,fib_sce$sub_celltype_genes=="Fib(C4)-CXCL8"], group_by = "set", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = sub_sce.x , group_by = "set",nlabel = 10,
            features_label= c("PKHD1",  "DCDC2",  "ACSM3" , "ELF3"  , "CFTR"   ,"SORBS2", "CXCL8" , "BICC1",  "SPP1"  , "SLC4A4"),
            #x_metric = "avg_log2FC", 
            DE_threshold = "abs(avg_log2FC) > 0.25 & p_val_adj < 0.05") & scale_color_gradientn(colours=brewer.pal(n=8, name="PuBuGn"), guide = "colourbar") 
ggsave(file=file.path(Figure4,"VolcanoPlot_CRC_LM.pdf"),width = 12,height = 7)


###


sub_sce.y<- RunDEtest(srt =fib_sce, group_by = "sub_celltype_genes", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = sub_sce.y , group_by = "sub_celltype_genes",nlabel = 10,
            #features_label= c("PKHD1",  "DCDC2",  "ACSM3" , "ELF3"  , "CFTR"   ,"SORBS2", "CXCL8" , "BICC1",  "SPP1"  , "SLC4A4"),
            #x_metric = "avg_log2FC", 
            DE_threshold = "abs(avg_log2FC) > 0.25 & p_val_adj < 0.05") & scale_color_gradientn(colours=brewer.pal(n=8, name="PuBuGn"), guide = "colourbar") 
ggsave(file=file.path(Figure4,"VolcanoPlot_fib.pdf"),width = 14,height = 9)




degs<-sub_sce.x@tools$DEtest_sub_celltype_genes$AllMarkers_wilcox
DoHeatmap(subset(fib_sce, downsample = 1000),
          size = 2.5,
          features=sce.markers.x.sig.10$gene,
          #slot = 'counts', 
          raster = T,
          #group.colors =stcol
          ) + 
  scale_fill_gradient(
    low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
    #mid = "white",
    high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
    #midpoint = 1,
    guide = "colourbar",
    aesthetics = "fill"
  )


##
Idents(fib_sce)<-fib_sce$set
sce.markers.group<- FindAllMarkers(fib_sce, 
                             only.pos = TRUE, 
                             min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.group.sig<-sce.markers.group[sce.markers.group$p_val_adj<0.001,]
sce.markers.group.sig<- sce.markers.group %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)# top 10 
write.csv(sce.markers.group.sig,file=file.path(Figure4,"Fibroblasts_markergenes_split_sce.markers.group_LM.csv"))




a1<-sce.markers.x.sig[sce.markers.x.sig$cluster=="Fib(C4)-CXCL8",]

table(sce.markers.x.sig$cluster)

markergeness<-intersect(a,b)
markergenes<-intersect(sce.markers.x.sig[sce.markers.x.sig$cluster=="Fib(C4)-CXCL8",]$gene,sce.markers.group[sce.markers.group$cluster=="LM",]$gene)
markergenes



a<-unique(sce.markers.x.sig$gene[sce.markers.x.sig$cluster=="Fib(C4)-CXCL8"])
b<-sce.markers.group.sig[sce.markers.group.sig$cluster=="LM",]$gene
# Remove NAs
a <- na.omit(a)
b <- na.omit(b)
fib_genes_list<-list(a,b)
names(fib_genes_list)<-c("scRNA_Fib_C4","scRNA_Fib_LM")


ggvenn(fib_genes_list)
# 创建维恩图
library(VennDiagram)
# 创建维恩图
fib_genes_list <- list(scRNA_Fib_C4 = a, scRNA_Fib_LM = b)

# Plot Venn diagram
venn.plot <- venn.diagram(
  x = fib_genes_list,
  category.names = names(fib_genes_list),
  filename = NULL,
  output = TRUE,
  imagetype = "pdf",
  height = 480, 
  width = 480, 
  resolution = 300,
  lwd = 2,
  fill = c( '#2ca02c', "red"),
  alpha = 0.5,
  label.col = "black",
  cex = 1.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "outer",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.dist = c(0.05, 0.05),  # Two values for two sets
  cat.pos = c(-20, 20)  # Positions for two categories
)

# 保存图像
pdf("LM_fib_marker_VennDiagram.pdf")
grid.draw(venn.plot)
dev.off()

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
                       pvalueCutoff  = 0.05,qvalueCutoff  = 0.05, readable = TRUE)
  }
  return(enrich)
}

ProSig <- new_cell_markers_genes####设置基因集
#Sigbp <- lapply(ProSig, Myenrich, category = "gobp", geneid = "SYMBOL")
#Sigkegg <- lapply(ProSig, Myenrich, category = "kegg", geneid = "SYMBOL")

library(ggnewscale)
library("clusterProfiler")
options(connectionObserver = NULL)
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")



Sigbp<-Myenrich(markergenes,category = "gobp", geneid = "SYMBOL")

pdf(file=file.path(Figure4,"Cluster_Go-up.pdf"),width = 6,height = 8)
dotplot(Sigbp,showCategory = 10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free') + scale_color_continuous(low='purple', high='green')+ aes(shape=GeneRatio > 0.1)
dev.off()

pdf(file=file.path(Figure4,"Cluster_Go-net.pdf"),width = 9,height = 8)
Sigbp2<- pairwise_termsim(Sigbp)
emapplot(Sigbp2, showCategory = 10,color = "p.adjust")+ scale_color_continuous(low='purple', high='green')
dev.off()

pdf(file=file.path(Figure4,"Cluster_Go-circ.pdf"),width = 15,height = 8)
cnetplot(Sigbp2, showCategory = 15,categorySize="count",circular = FALSE, colorEdge = TRUE,cex_gene = 0.4,cex_label_gene =0.5)
dev.off()

pdf(file=file.path(Figure4,"Cluster_Go-circ.pdf"),width = 8,height = 5)
setReadable(Sigbp2, 'org.Hs.eg.db', 'ENTREZID') 
test<-c(grep("gly",Sigbp$Description,value=TRUE),grep("glu",Sigbp$Description,value=TRUE))
#test<-c("ATP generation from ADP","ATP metabolic process","glycolytic process","generation of precursor metabolites and energy")
#cnetplot(Go2, showCategory = test, categorySize="count",circular = TRUE, colorEdge = TRUE)
cnetplot(Sigbp2, showCategory = test, categorySize="count",circular = FALSE, colorEdge = TRUE,cex_gene = 0.4,cex_label_gene =0.5)
dev.off()


write.csv(x,file=file.path(Figure4,"fib_stRNA.csv"))


##展示所有 富集信息
x<-Sigbp
# Convert to data frame and set BP description as rownames
x = as.data.frame(x, row.names = x$Description)

l = gsoap_layout(x,
                 genes = 'geneID',
                 pvalues = 'p.adjust',
                 projection = 'tsne',
                 scale.factor = 0.8,
                 no.clusters = 5)

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
#??gsoap_plot


