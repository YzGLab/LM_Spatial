###scp for map. 
##
library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)

##############
sce<-readRDS("data/Colon_HC_combined.sce.rds")
load("data/totalscRNAdata.Rdata")
######################################################
sce<-AddMetaData(sce,metadata)
mycol2<-c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F", "#E41A1C", '#FFD92F', "#C3BC3F", "#FF7F00", "#1B9E77", "#66A61E", "#F1788D")

mycol2<-c( '#A65628', '#984EA3','#FDB462','#BEAED4','#A6D854','#1B9E77','#FFD92F','#E5C494','#8DA0CB',
           '#E41A1C','#FF7F00','#377EB8','#E78AC3',
           '#FFFF33', '#A65628', '#F781BF', '#999999',
           '#8DD3C7', '#BEBADA', '#FB8072','#80B1D3','#B3DE69', '#FCCDE5')
##########################################################################################################################
p1<-DimPlot(sce,reduction = "umap",label=T,group.by = "sub_celltype")
p2<-DimPlot(sce,reduction = "umap",label=T)
p1+p2
#################try
colnames(sce@meta.data)
ggscplot(object = sce) +
  geom_scPoint2(aes(color = sub_celltype_genes,
                    cluster = sub_celltype_genes,
                    cluster_anno = sub_celltype),
                show.legend = F,
                add_label = F,
                add_legend = T,
                lgd_x = 1.2) +
  theme_sc(r = 0.3)

#######################################
DefaultAssay(sce) <- "integrated"
sce<- FindNeighbors(sce, reduction = "harmony", dims = 1:20)
sce <- FindClusters(sce, resolution =1)
sce <- RunUMAP(sce,reduction = "harmony", dims = 1:20)
sce<- RunTSNE(sce, reduction = "harmony", dims = 1:20)
DimPlot(sce,reduction = "umap",label=T)
DimPlot(sce,reduction = "umap",label=T,group.by="orig.ident")
DimPlot(sce,reduction = "umap",label=T,split.by="orig.ident")
DimPlot(sce,reduction = "tsne",label=T,split.by="orig.ident")
################################### intital checking of cell types

library(SCP)
sce@meta.data[,c(2:25)]<-NULL
sce<-AddMetaData(sce,metadata)
sce$Patients<-ifelse(sce$orig.ident =="gan1","LM1",ifelse(sce$orig.ident=="gan2","LM2",ifelse(sce$orig.ident=="chang1","CRC1","CRC2")))
sce$set<-ifelse(sce$orig.ident %in% c("chang1","chang2"),"CRC","LM")
sce$set<-factor(sce$set,levels = c("CRC","LM"))#
sce$sub_celltype<-factor(sce$sub_celltype,levels = c("CD4","CD8", "Treg", "NK","MAIT","Myeloid","B","Plasma", "Endothelial", "Epithelial" , "Fibroblasts"))
metadata<-sce@meta.data
####################################
CellDimPlot(
  srt =sce, group.by =  "seurat_clusters", stat_plot_size = 0.1,
  reduction = "UMAP", theme_use = "theme_blank",
  palcolor = mycol2,
)
DimPlot(sce,reduction = "umap",label=T,group.by="orig.ident")


CellDimPlot(
  srt =sce, group.by =  "sub_celltype",  stat_plot_size = 0.1,
  reduction = "UMAP", theme_use = "theme_blank",show_stat = T, label = TRUE,
  palcolor = mycol2,
)
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_total.pdf"),width = 9.5,height = 5.03)

CellDimPlot(
  srt =sce, group.by =  "set",  stat_plot_size = 0.1,
  reduction = "UMAP", theme_use = "theme_blank",
  palcolor =c('#E41A1C','#377EB8'),
)
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM.pdf"),width = 9.5,height = 5.03)


CellDimPlot(
  srt =sce, group.by =  "Patients",  stat_plot_size = 0.1,split.by = "immunecells",
  reduction = "UMAP", theme_use = "theme_blank",
  palcolor =c('#E41A1C','#377EB8'),
)

ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_Patients.pdf"),width = 11.5,height = 8.03)

CellDimPlot(
  srt =sce, group.by =  "sub_celltype",  stat_plot_size = 0.1,split.by = "Patients",
  reduction = "UMAP", theme_use = "theme_blank",
  palcolor = mycol2,
)
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_Patients_split.pdf"),width = 11.5,height =9)
######比例：
CellStatPlot(sce,stat.by ="sub_celltype", group.by ="set",
             plot_type = "dot",palcolor = mycol2)
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_Patients_probability.pdf"),width = 3.8,height = 5)
####
######比例：
CellStatPlot(sce,stat.by ="sub_celltype", group.by ="set",
             plot_type = "dot",palcolor = mycol2)
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_source_probability.pdf"),width = 3.8,height = 8)

###table ratio
prop.table(table(sce$sub_celltype,sce$Patients),1)
library(moonBook) 
library(autoReg)
library(tidyverse)
library(rrtable)
metadata<-sce@meta.data
baseline <- gaze(Patients~sub_celltype,data=metadata) %>% 
  myft()
baseline
table2docx(baseline)
################################## gene plot_deg 
DefaultAssay(sce)<-"integrated"
Idents(sce)<-sce$sub_celltype
sce.markers <- FindAllMarkers(sce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sce.markers ,file = file.path(Figure1,"full_scemarkers.csv"))
sce.markers.sig <- sce.markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10 
print(sce)
##
library(scRNAtoolVis)
sce_vis <- RunDEtest(srt = sce, group_by = "sub_celltype", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = sce_vis, group_by = "sub_celltype")
###############
DimPlot(sce,reduction = "umap",label=T)
###########################################check the marker genes in scRNA 
VlnPlot(sce,features = "CD8A",pt.size = 0)
VlnPlot(sce,features = "CD4",pt.size = 0)

VlnPlot(sce,features = "FOXP3",pt.size = 0)
VlnPlot(sce,features = "CD79A",pt.size = 0)
VlnPlot(sce,features = "SLC4A10",pt.size = 0 )
VlnPlot(sce,features = "KLRF1",pt.size = 0)

VlnPlot(sce,features ="EPCAM",pt.size = 0 )##"Epi cells
VlnPlot(sce,features ="PECAM1",pt.size = 0)##"Endo cells

FeaturePlot(sce,features ="CD3E" )##"T
FeaturePlot(sce,features ="CD4" )
FeaturePlot(sce,features ="CD8A" )

FeaturePlot(sce,features ="FOXP3" )##"Treg cells
FeaturePlot(sce,features ="CD79A" )##"B cells

FeaturePlot(sce,features ="SLC4A10" )##"MAIT cells
FeaturePlot(sce,features ="KLRF1" )##"NK cells

FeaturePlot(sce,features ="LYZ" )##"Myeloid cells
FeaturePlot(sce,features ="FCGR3A" )##"Myeloid cells

FeaturePlot(sce,features ="MZB1" )##"Plasma cells
FeaturePlot(sce,features ="FCGR3B" )##"Plasma cells

FeaturePlot(sce,features ="EPCAM" )##"Epi cells
FeaturePlot(sce,features ="CDH1" )##"Epi cells

FeaturePlot(sce,features ="PECAM1" )##"Endo cells
FeaturePlot(sce,features ="CDH5" )##"Endo cells

FeaturePlot(sce,features ="COL1A1" )##"fib
#########
##

DefaultAssay(sce)<-"RNA"
FeatureDimPlot(
  srt = sce, features = c("CD4", "CD8A", "FOXP3","KLRF1","SLC4A10","LYZ","CD79A","MZB1","PECAM1","EPCAM","COL1A1"),ncol=3,
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_cellmarkers_combined_features.pdf"),width = 11.58,height =9.5)

FeatureDimPlot(
  srt = sce, features = c("CD4", "CD8A", "FOXP3","KLRF1","SLC4A10","LYZ","CD79A","MZB1","PECAM1","EPCAM","COL1A1"),
  compare_features = TRUE, label = TRUE, label_insitu = TRUE,
  add_density=TRUE,
  palcolor=mycol2,
  #label_segment_color =mycol2,
  lineages_palcolor = mycol2,
  reduction = "UMAP", theme_use = "theme_blank"
)
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_cellmarkers.pdf"),width = 11.58,height =8.03)

###
CellStatPlot(sce,stat.by ="sub_celltype", group.by =c("Patients","set"),
             plot_type = "dot",palcolor = mycol2)
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_patients_precent_combined.pdf"),width = 9.5,height = 5.03)

CellStatPlot(sce, stat.by ="sub_celltype", group.by = "set",  stat_type = "percent", plot_type = "rose",palcolor = mycol2)
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_patients_precent_combined_rose.pdf"),width = 9.5,height = 5.03)
###
library(scales)
library(tidyr)
library(tidyverse)
library(patchwork)
colnames(sce@meta.data)
data_plotB <- table(sce[,sce$immunecells!="immunecells"]@meta.data$Patients,sce[,sce$immunecells!="immunecells"]@meta.data$sub_celltype) %>% melt()
colnames(data_plotB) <- c("Cell_subtype", "NMFcluster","Number")

data_plotB$NMFcluster<-as.factor(data_plotB$NMFcluster)

pC1 <- ggplot(data = data_plotB, aes(x = Cell_subtype, y = Number, fill =NMFcluster )) +
  geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="fill")+
  #scale_fill_manual(values=c("#1B9E77", "#377EB8","#E41A1C")) +
  scale_fill_manual(values=mycol2) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Average number")+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  #theme(axis.title.y=element_blank(),
  #  axis.ticks.y=element_blank(),
  # axis.text.y=element_blank()  )+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) + 
  coord_flip()+
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
  scale_y_reverse()+theme(legend.position="none")
#

pC1

data_plotC <- table(sce[,sce$immunecells=="immunecells"]@meta.data$Patients,sce[,sce$immunecells=="immunecells"]@meta.data$sub_celltype) %>% melt()
colnames(data_plotC) <- c("Cell_subtype", "NMFcluster","Number")
data_plotC$NMFcluster<-as.factor(data_plotC$NMFcluster)

pC2 <- ggplot(data = data_plotC, aes(x = Cell_subtype, y = Number, fill =NMFcluster )) +
  geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="fill")+
  #scale_fill_manual(values=c("#1B9E77", "#377EB8","#E41A1C")) +
  scale_fill_manual(values=mycol2) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Average number")+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  #theme(axis.title.y=element_blank(),
  #  axis.ticks.y=element_blank(),
  # axis.text.y=element_blank()  )+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6)) + 
  coord_flip()+
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "left") +
  scale_y_reverse()#theme(legend.position="none")
#

pC2
pC1+pC2
pC <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,3),guides = 'collect')
pC
ggsave(file.path(Figure1,"celltype_CRC_LM_number.pdf"),width =7.98,height = 3.82 )
#############################


pC2<- ggplot(data = data_plotC, aes(x =NMFcluster , y = Number, fill =Cell_subtype )) +
  geom_bar(stat = "identity", width=0.8,aes(group=NMFcluster),position="fill")+
  scale_fill_manual(values=c("#1B9E77", "#377EB8","#E41A1C")) +
  theme_bw()+coord_flip()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  scale_y_continuous(labels = percent)+  ####用来将y轴移动位置
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.8)) +theme(legend.position="right")
# + scale_x_discrete(position = "top") 
#让横轴上的标签倾斜45度
pC2 
#############################
library(scRNAtoolVis)
clusterCornerAxes(object = sce,reduction = 'umap',
                  clusterCol = "sub_celltype",
                  #noSplit = T,
                  stripCol=mycol2,
                  arrowType = "open", noSplit =T,
                  cornerTextSize = 3.5,
                  themebg = 'bwCorner',
                  addCircle = TRUE,
                  cicAlpha = 0.1,
                  nbin = 20) 
h2<-AverageHeatmap(object =sce[,sce$set=="CRC"],
                   markerGene = c("CD4", "CD8A", "FOXP3","KLRF1","SLC4A10","LYZ","CD79A","MZB1","PECAM1","EPCAM","COL1A1"),
                   border = FALSE,
                   group.by = "sub_celltype",
                   annoCol = TRUE,
                   myanCol=c("CD4"="#A65628" ,"CD8"="#984EA3" ,"Treg"="#FDB462","NK"="#BEAED4","MAIT"= "#A6D854","Myeloid"= "#1B9E77",
                             "B"= "#FFD92F","Plasma"= "#E5C494","Endothelial"= "#8DA0CB","Epithelial"= "#E41A1C","Fibroblasts"= "#FF7F00"),
                   htCol = c('#377EB8', "white","#E41A1C"))
h2

##
table(sce$sub_celltype)
mycol2
h3<-AverageHeatmap(object =sce[,sce$set=="LM"],
                   markerGene = c("CD4", "CD8A", "FOXP3","KLRF1","SLC4A10","LYZ","CD79A","MZB1","PECAM1","EPCAM","COL1A1"),
                   border = FALSE,
                   group.by = "sub_celltype",
                   annoCol = TRUE,
                   myanCol=c("CD4"="#A65628" ,"CD8"="#984EA3" ,"Treg"="#FDB462","NK"="#BEAED4","MAIT"= "#A6D854","Myeloid"= "#1B9E77",
                             "B"= "#FFD92F","Plasma"= "#E5C494","Endothelial"= "#8DA0CB","Epithelial"= "#E41A1C","Fibroblasts"= "#FF7F00"),
                   htCol = c('#377EB8', "white","#E41A1C"))
h3
h2+h3
pdf(file.path(Figure1,"celltype_immune_CRC_LM_markergenes_combined_heatmap.pdf"),width = 6,height = 5.03)
h2+h3
dev.off()


ecols <- setNames(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"), c("0", 1:11))
VlnPlot(sce[,sce$set=="LM"], features=c("CD4", "CD8A", "FOXP3","KLRF1","SLC4A10","CD14","CD79A","MZB1","PECAM1","EPCAM","COL1A1"), 
                 stack = T, fill.by='ident', cols =mycol2, pt.size=0) +
  NoLegend() 
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_markergenes_combined_vlnplot_LM.pdf"),width = 6,height = 5.03)

VlnPlot(sce[,sce$set=="CRC"], features=c("CD4", "CD8A", "FOXP3","KLRF1","SLC4A10","CD14","CD79A","MZB1","PECAM1","EPCAM","COL1A1"), 
           stack = T, fill.by='ident', cols =mycol2, pt.size=0) +
  NoLegend() 
ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_markergenes_combined_vlnplot_CRC.pdf"),width = 6,height = 5.03)
countexp.Seurat<-NULL
metabolism_assay<-NULL


ggsave(file = file.path(Figure1,"celltype_immune_CRC_LM_markergenes_combined.pdf"),width = 6,height = 5.03)
###########
p_lm<-DoHeatmap(subset(sce[,sce$set=="LM"], downsample = 1000),
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
p_CRC<-DoHeatmap(subset(sce[,sce$set=="CRC"], downsample = 1000),
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
############################### figure1 
###
###########

