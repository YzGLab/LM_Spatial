#分群热图
library(ggSCvis)
library(scales)
library(tidyr)
library(tidyverse)
library(patchwork)
DefaultAssay(sub_sce)<- "metabolits"


##
#######
sub_sce$cluster3<-factor(sub_sce$cluster3)
Idents(sub_sce)<-sub_sce$cluster3
sce.markers.x <- FindAllMarkers(sub_sce, 
                                only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.sig <- sce.markers.x %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10

write.csv(sce.markers.x,file=file.path(Figure7,paste0(sample,"cluster3.csv")))##
###
sub_sce$cluster6<-factor(sub_sce$cluster6)
Idents(sub_sce)<-sub_sce$cluster6
sce.markers.x <- FindAllMarkers(sub_sce, 
                                only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.sig <- sce.markers.x %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10

write.csv(sce.markers.x,file=file.path(Figure7,paste0(sample,"cluster6.csv")))
############
#######
sub_sce$cluster8<-factor(sub_sce$cluster8)
Idents(sub_sce)<-sub_sce$cluster8
sce.markers.x <- FindAllMarkers(sub_sce, 
                                only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.sig <- sce.markers.x %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10

write.csv(sce.markers.x,file=file.path(Figure7,paste0(sample,"cluster8.csv")))
#######


##
library(data.table)

DefaultAssay(sce)



sub_sce<-sce[,sce@meta.data$orig.ident=='LM_2']
#############

Idents(sub_sce)<-sub_sce$st_celltype
sce.markers.celltype <- FindAllMarkers(sub_sce, 
                                       only.pos = TRUE, test.use = "wilcox",
                                       #idents.1="Fibroblast",idents.2="Hepatocytes",
                                       min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.sig <- sce.markers.celltype %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC)# top 10
sce.markers.sig <-sce.markers.sig [sce.markers.sig$p_val_adj<0.001,]

write.csv(sce.markers.celltype,file=file.path(Figure8,paste0(sample,"_full_st_celltype.csv")))



DoHeatmap(#subset(sub_sce, downsample = 1000),
  sub_sce,
  size = 2.5,
  features=sce.markers.sig$gene,
  slot = 'counts', raster = T,group.colors =stcol1) + 
  scale_fill_gradient(
    low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
    #mid = "white",
    high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
    #midpoint = 1,
    guide = "colourbar",
    aesthetics = "fill"
  )

###########GroupHeatmap
library(SCP)
names(stcol1)<- c('Normal Epi', 'Tumor', 'Fibroblast', 'Hepatocytes', 'Lamina propria', 'B/Plasma', 'Monocyte')
stcol1 <- c('Normal Epi' = '#1f77b4', 
            'Tumor' = '#ff7f0e', 
            'Fibroblast' = '#2ca02c', 
            'Hepatocytes' = '#d62728', 
            'Lamina propria' = '#9467bd', 
            'B/Plasma' = '#f5e801', 
            'Monocyte' = '#08f7f0')

ht1<-GroupHeatmap(sub_sce, feature_split = sce.markers.sig$cluster,
                  feature_split_palcolor = stcol1,
                  #cell_annotation_palcolor=stcol1,
                  features = sce.markers.sig$gene,
                  heatmap_palette = "YlOrRd",
                  fill.by="st_celltype",
                  #group.by ="cluster12",
                  group.by ="st_celltype",
                  add_dot = TRUE, 
                  add_bg = TRUE,
                  features_fontsize = 6,label_size = 6,keys_fontsize = c(4, 6),
                  nlabel = 0, 
                  show_row_names = TRUE)
ht1$plot
ggsave(file = file.path(Figure7,paste0(sample,"sp_GroupHeatmap_top_clusters_metabolits_10.pdf")),width = 6,height =12.03)
ht2<-GroupHeatmap(sub_sce, feature_split = sce.markers.sig$cluster,
                  feature_split_palcolor = stcol1,
                  cell_annotation_palcolor=stcol1,
                  features = sce.markers.sig$gene,
                  heatmap_palette = "YlOrRd",
                  #fill.by="st_celltype",
                  #group.by ="cluster6",
                  group.by ="st_celltype",
                  #exp_method =   "log1p",
                  #add_dot = TRUE, 
                  add_bg = TRUE,
                  terms_fontsize = 6,
                  exp_cutoff = 50,
                  flip = TRUE,
                  nlabel = 0, 
                  show_row_names = TRUE)
ht2$plot
ggsave(file = file.path(Figure7,paste0(sample,"sp_stcelltype_GroupHeatmap_topmetabolits_10.pdf")),width = 6,height =12.03)




##代谢通路