###高变代谢物的寻找
#spatial 
library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)

#save(sce.cell_type_markers,sce.markers.epi,sce.markers.epi.sig,file='data/spatial_scRNA_differetgenes_clusters.Rdata')
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
load('data/spatial_scRNA_differetgenes_clusters.Rdata')# from before differen analysis
DefaultAssay(st_sce)<-"SCT"
print(st_sce)
metadata<-st_sce@meta.data
#######################
DefaultAssay(st_sce)<-"metabolism"
#st_sce<-st_sce[,st_sce$orig.ident %in% c("LM_2","C_2")]

st_sce@assays$metabolism@data[1:5,1:5]
#tumor vs other cells  in LM
#tumor vs other cells  in LM  
sub_st_sce<-st_sce[,st_sce$set=="LM" & st_sce$st_celltype %in% c("Normal Epi","Tumor") ]

Idents(sub_st_sce)<-sub_st_sce$st_celltype
#all.markers <- FindAllMarkers(object = sub_st_sce)
st_tumor_LM.markers<- FindAllMarkers(sub_st_sce, slot="counts",
                                  only.pos = TRUE, 
                                  logfc.threshold = 0.1,min.pct=0.1)
table(st_tumor_LM.markers$cluster)
st_tumor_LM.markers.sig <- st_tumor_LM.markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10
#### 显示无明显代谢物 
##
#######################################################################
######scp pipline
#sub_sc_sce<- RunDEtest(srt =sub_sc_sce, group_by = "sub_celltype_genes", fc.threshold = 1, only.pos = FALSE)
#VolcanoPlot(srt = sub_sc_sce, group_by = "sub_celltype_genes")


#tumor vs other cells  in CRC
st_sce@assays$metabolism@data[1:5,1:5]
sub_st_sce<-st_sce[,st_sce$set=="CRC" ]
Idents(sub_st_sce)<-sub_st_sce$st_celltype
st_tumor_CRC.markers<- FindAllMarkers(sub_st_sce, slot="counts",
                                  only.pos = TRUE, 
                                  min.pct = 0.25, logfc.threshold = 0.25)

table(st_tumor_CRC.markers$cluster)
st_tumor_CRC.markers.sig <- st_tumor_CRC.markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC)# top 10
#######################################################################

### LM tumor vs CRC tumor 
#tumor vs other cells  in LM
sub_st_sce<-st_sce[,st_sce$st_celltype=="Tumor"]
Idents(sub_st_sce)<-sub_st_sce$set
st_tumor_LM_vs_CRC.markers<- FindAllMarkers(sub_st_sce, slot="counts",
                                      only.pos = TRUE, 
                                      min.pct = 0.25, logfc.threshold = 0.25)
table(st_tumor_LM_vs_CRC.markers$cluster)
st_tumor_LM_vs_CRC.markers.sig <- st_tumor_LM_vs_CRC.markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC)# top 10
##########################
### LM_nontumor vs CRC_Notumor 
#tumor vs other cells  in LM

table(st_sce$st_celltype,st_sce$set)
sub_st_sce<-st_sce[,st_sce$st_celltype %in% c("Normal Epi")]
table(sub_st_sce$set)
Idents(sub_st_sce)<-sub_st_sce$set
??FindAllMarkers
DefaultAssay(sub_st_sce)<-"metabolism"
st_notumor_LM_vs_CRC.markers<- FindAllMarkers(sub_st_sce,  slot="counts",
                                            only.pos = TRUE, 
                                            min.pct = 0.15, logfc.threshold = 0.1)
table(st_notumor_LM_vs_CRC.markers$cluster)
st_notumor_LM_vs_CRC.markers.sig <- st_notumor_LM_vs_CRC.markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC)# top 10
###交叉
st_tumor_CRC.markers
intersect(st_tumor_CRC.markers$gene,st_tumor_LM_vs_CRC.markers$gene)
