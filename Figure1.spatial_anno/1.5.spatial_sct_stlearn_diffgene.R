###diff of lM vs CRC for st _celltypes
library(SCP)
library(Seurat)
library(tidydr)
library(data.table)
library(tidyverse)
###############
load("data/stmetadata.Rdata")
st_sce<-readRDS("data/Colon_HC_spatial.RDS")

#update metadata information from leiden SEM clusters
table(stmetadata$st_celltype)
#stmetadata<-st_sce@meta.data
stmetadata$st_celltype<-stmetadata$leiden
stmetadata$st_celltype[stmetadata$leiden=="0"]="Normal_Epi/Endo"#
stmetadata$st_celltype[stmetadata$leiden=="1"]="Tumor"# ok
stmetadata$st_celltype[stmetadata$leiden=="2"]="Fibroblast"# OK
stmetadata$st_celltype[stmetadata$leiden=="3"]="Hepatocytes"#OK
stmetadata$st_celltype[stmetadata$leiden=="4"]="Lamina propria"##OK
stmetadata$st_celltype[stmetadata$leiden=="5"]="B/Plasma"
stmetadata$st_celltype[stmetadata$leiden=="6"]="Monocyte"# "ADM"
st_sce<-AddMetaData(st_sce,stmetadata)
save(stmetadata,file="data/stmetadata.Rdata")

##########
DefaultAssay(st_sce)<-"SCT"
# 假设'seurat_obj'是您的Seurat对象
# 首先，找出包含'MT-'或'RP'的基因
mt_genes <- grep("^MT-", rownames(st_sce), value = TRUE)
rp_genes <- grep("^RP[L|S]", rownames(st_sce), value = TRUE)
# 合并MT和RP基因列表
genes_to_remove <- union(mt_genes, rp_genes)
# 从Seurat对象中移除这些基因
sub_sce<- subset(st_sce, features = setdiff(rownames(st_sce), genes_to_remove))##########################
saveRDS(sub_sce,file="data/spatial_combinend_removed_MT_RP.RDS")
### Define a function to generate differential expression plots between LM and CRC
diffplot<-function(x){
  sub_sce.x <- RunDEtest(srt =sub_sce[,sub_sce$st_celltype==x], group_by = "set", fc.threshold = 1, only.pos = FALSE)
  VolcanoPlot(srt = sub_sce.x , group_by = "set")
  x<-gsub('/',"_",x)
  ggsave(file = file.path(Figure2,paste0("CRC_LM_patients_diffgenes_",x,".pdf")),width =10.97,height = 3)####
  sce.markers.x <- FindAllMarkers(sub_sce.x, 
                                          only.pos = TRUE, 
                                          min.pct = 0.25, logfc.threshold = 0.25)
  
  write.csv(sce.markers.x,file =file.path(Figure2,paste0("CRC_LM_patients_diffgenes_",x,".csv")) )
}

sapply(c( 'Normal_Epi/Endo', 'Tumor', 'Fibroblast', 'B/Plasma', 'Monocyte'),diffplot)
###################################################################################################################
DefaultAssay(sub_sce)<- "SCT"
Idents(sub_sce)<-sub_sce$st_celltype
sce.markers <- FindAllMarkers(sub_sce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.x<-sce.markers[sce.markers$p_val_adj<0.001,]

sce.markers.sig <- sce.markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10
write.csv(sce.markers.x,file=file.path(Figure1,"sce.markers_all_celltypes.csv"))
write.csv(sce.markers.sig,file=file.path(Figure1,"sce.markers_all_celltypes_top10.csv"))
write.csv(stmetadata,file=file.path(Figure1,"stmetadata_stRNA.csv"))



# Generate heatmap for selected markers
DoHeatmap(subset(sub_sce, downsample = 1000),
          size = 2.5,
          features=sce.markers.sig$gene,
          slot = 'counts', raster = T,group.colors =stcol) + 
  scale_fill_gradient(
    low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
    #mid = "white",
    high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
    #midpoint = 1,
    guide = "colourbar",
    aesthetics = "fill"
  )

# Generate feature statistics plot (box plot)
FeatureStatPlot(sub_sce,
                stat.by = sce.markers.sig$gene,
                plot_type = "box",
                group.by = "st_celltype",
                bg.by ="st_celltype", 
                stack = TRUE, 
                flip = TRUE,
                 legend.direction = "horizontal"
                
)
###########GroupHeatmap
names(stcol1)<- c('Normal_Epi/Endo', 'Tumor', 'Fibroblast', 'Hepatocytes', 'Lamina propria', 'B/Plasma', 'Monocyte')
stcol1 <- c('Normal_Epi/Endo' = '#1f77b4', 
           'Tumor' = '#ff7f0e', 
           'Fibroblast' = '#2ca02c', 
           'Hepatocytes' = '#d62728', 
           'Lamina propria' = '#9467bd', 
           'B/Plasma' = '#f5e801', 
           'Monocyte' = '#08f7f0')
ht1<-GroupHeatmap(sub_sce, feature_split = sce.markers.sig$cluster,feature_split_palcolor = stcol1,cell_annotation_palcolor=stcol1,
                    features = sce.markers.sig$gene,heatmap_palette = "YlOrRd",fill.by="st_celltype",
                    group.by ="st_celltype",add_dot = TRUE, add_bg = TRUE, nlabel = 0,row_names_rot = 0,
                  row_title_rot = 90, show_row_names = TRUE)
ht1

ggsave(file = file.path(Figure1,"CRC_LM_markergenes_GroupHeatmap_topgenes10.pdf"),width = 5.8,height =12.03)


FeatureStatPlot(sub_sce,
                 stat.by = sce.markers.sig$gene,
                 plot_type = "dot",
                 group.by = "st_celltype",
                 bg.by ="st_celltype", 
                 stack = TRUE, 
                 flip = TRUE
) %>% panel_fix_overall(width = 8, height = 5) # As the plot is created by combining, we can adjust the overall height and width directly.



VlnPlot(sub_sce, features=sce.markers.sig$gene, 
        stack = T, fill.by='ident', cols =stcol, pt.size=0) 


###########
p_lm<-DoHeatmap(subset(st_sce[,st_sce$set=="LM"], downsample = 1000),
                size = 2.5,
                features=c("CD4", "CD8A", "FOXP3","KLRF1","SLC4A10","CD14","CD79A","MZB1","PECAM1","EPCAM","COL1A1"),
                slot = 'counts', raster = T,group.colors = mycol2) + 
  scale_fill_gradient(
    low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
    #mid = "white",
    high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
    #midpoint = 1,
    guide = "colourbar",
    aesthetics = "fill"
  )+NoLegend()
p_CRC<-DoHeatmap(subset(st_sce[,st_sce$set=="CRC"], downsample = 1000),
                 size = 2.5,
                 features=c("CD4", "CD8A", "FOXP3","KLRF1","SLC4A10","CD14","CD79A","MZB1","PECAM1","EPCAM","COL1A1"),
                 slot = 'counts', raster = T,group.colors = mycol2) + 
  scale_fill_gradient(
    low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
    #mid = "white",
    high = rev(c("#E41A1C", '#ef8a62', '#fddbc7')),
    #midpoint = 1,
    guide = "colourbar",
    aesthetics = "fill"
  )+NoLegend()
p_CRC+p_lm
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_markergenes_combined.pdf"),width = 6,height = 5.03)
