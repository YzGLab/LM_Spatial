##
library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)

#save(sce.cell_type_markers,sce.markers.epi,sce.markers.epi.sig,file='data/spatial_scRNA_differetgenes_clusters.Rdata')
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
##############################
#load('data/spatial_differetgenes_clusters.Rdata')# from before differen analysis

##############
Idents(sub_sce)<-sub_sce$st_celltype
sce.markers1 <- FindAllMarkers(sub_sce[,sub_sce$set=="CRC"], 
                               only.pos = TRUE, 
                               min.pct = 0.5, logfc.threshold = 0.5)
sce.markers2 <- FindAllMarkers(sub_sce[,sub_sce$set=="LM"], 
                               only.pos = TRUE, 
                               min.pct = 0.5, logfc.threshold = 0.5)
#sce.markers.sig <- sce.markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10

sce.markers1$cluster<-paste0("CRC_",sce.markers1$cluster)
sce.markers2$cluster<-paste0("LM_",sce.markers2$cluster)
sce.markers<-rbind(sce.markers1,sce.markers2)################################
sce.cell_type_markers<-split(sce.markers$gene,sce.markers$cluster)###############合并两种差异

######计算单细胞差异基因 

nonimmunecells_markers<-read.csv("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure1_overall/nonimmunecells_markergenes_list.csv")


library(dplyr)
library(stringr)
nonimmunecells_markers$cluster<-factor(nonimmunecells_markers$cluster,
                                  levels=c("C0-Endothelial", "C1-Endothelial",  "C2-Endothelial", 
                                           "C3-Endothelial","C4-Endothelial", "C5-Endothelial",
                                           "C6-Endothelial","C7-Endothelial","C8-Endothelial","C9-Endothelial",
                                           ###
                                           "C0-Epithelial", "C1-Epithelial","C2-Epithelial", "C3-Epithelial","C4-Epithelial",
                                           ############
                                           "C0-Fibroblasts", "C1-Fibroblasts","C2-Fibroblasts","C3-Fibroblasts", "C4-Fibroblasts", "C5-Fibroblasts",  "C6-Fibroblasts", "C7-Fibroblasts", 
                                           "C8-Fibroblasts"))
nonimmunecells_markers.list<-split(nonimmunecells_markers$gene,nonimmunecells_markers$cluster)
## non-tumor cells
#######转化top gene top 20 
sce.markers.sig_st_sc<-c(sce.cell_type_markers,nonimmunecells_markers.list)
####使用addmodulescore计算空间转录组学分数
DefaultAssay(st_sce)<-"SCT"
st_sce.copy<-AddModuleScore(st_sce,features=sce.markers.sig_st_sc,name=names(sce.markers.sig_st_sc))


###########cor of cell correlation
st_sc_metadata<-st_sce.copy@meta.data
colnames(st_sc_metadata)[c(32:67)]<-names(sce.markers.sig_st_sc)
####
library(corrplot)
library(RColorBrewer)
library(ggcorrplot)
library(ggthemes)
library(Hmisc)

display.brewer.all(type="seq")
brewer.pal.info[brewer.pal.info$category == "seq",]

correlation<-rcorr(as.matrix(st_sc_metadata[,c(32:67)]))
correlation_r<-correlation$r
correlation_p_values<-correlation$P
write.csv(correlation_r,file.path(Figure3,"corrplot_non_immune.csv"))
write.csv(correlation_p_values,file.path(Figure3,"corrplot_non_immune_pvalues.csv"))


cor.martix<-cor(st_sc_metadata[,c(32:67)])
cor.martix<-cor.martix[-c(1:22),c(1:12)]
corrplot(cor.martix)
corrplot(as.matrix(cor.martix), 
         method="circle", 
         #type="upper", 
         #col=colorRampPalette(c("yellow2","goldenrod","green"))(8), 
         col=brewer.pal(n=8, name="YlOrRd"), 
         tl.col="black", 
         tl.srt=45, 
         #diag=FALSE,
         #tl.pos='l',
         #p.mat = as.matrix(correlation_p_values[-c(1:12),c(1:12)]), 
         #sig.level = 0.00000000001, 
         insig = "label_sig")
pdf(file.path(Figure3,"corrplot_non_immune.pdf"),width =9.82, height = 10.74)
print(corrplot(as.matrix(cor.martix), 
               method="circle", 
               #type="upper", 
               #col=colorRampPalette(c("yellow2","goldenrod","green"))(8), 
               col=brewer.pal(n=8, name="PuBuGn"), 
               tl.col="black", 
               tl.srt=45, 
               #diag=FALSE,
               #tl.pos='l',
               #p.mat = as.matrix(correlation_p_values[-c(1:12),c(1:12)]), 
               #sig.level = 0.00000000001, 
               insig = "label_sig")
)
dev.off()
#######################
cor.martix.non_immune<-cor.martix


##immune cells
######计算单免疫细胞差异基因 

immunecells_markers<-read.csv("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure1_overall/immunecells_markergenes_list.csv")
st_sce<-readRDS("data/Colon_HC_spatial.RDS")

library(dplyr)
library(stringr)
immunecells_markers$cluster<-factor(immunecells_markers$cluster,
                                   levels=c("C0-CD4", "C1-CD4", "C2-CD4", "C3-CD4", "C4-CD4",
                                            "C0-Treg", "C1-Treg", 
                                            "C0-CD8", "C1-CD8", "C2-CD8", "C3-CD8", "C4-CD8", "C5-CD8", "C6-CD8",
                                            "C0-NK", "C1-NK",
                                            "C0-MAIT", "C1-MAIT",
                                            "C0-Myeloid", "C1-Myeloid", "C2-Myeloid", "C3-Myeloid", "C4-Myeloid", 
                                            "C0-B", "C1-B", "C2-B", "C3-B", "C0-Plasma", "C1-Plasma"))

immunecells_markers.list<-split(immunecells_markers$gene,immunecells_markers$cluster)
## non-tumor cells
#######转化top gene top 20 
sce.markers.sig_st_sc<-c(sce.cell_type_markers,immunecells_markers.list)
####使用addmodulescore计算空间转录组学分数
DefaultAssay(st_sce)<-"SCT"
st_sce.copy<-AddModuleScore(st_sce,features=sce.markers.sig_st_sc,name=names(sce.markers.sig_st_sc))


###########cor of cell correlation
st_sc_metadata<-st_sce.copy@meta.data
colnames(st_sc_metadata)[c(32:72)]<-names(sce.markers.sig_st_sc)
####
names(st_sc_metadata)
library(corrplot)
library(RColorBrewer)
display.brewer.all(type="seq")
brewer.pal.info[brewer.pal.info$category == "seq",]


correlation<-rcorr(as.matrix(st_sc_metadata[,c(32:72)]))
correlation_r<-correlation$r
correlation_p_values<-correlation$P
write.csv(correlation_r,file.path(Figure3,"corrplot__immune.csv"))
write.csv(correlation_p_values,file.path(Figure3,"corrplot__immune_pvalues.csv"))

cor.martix<-cor(st_sc_metadata[,c(32:72)])
cor.martix<-cor.martix[-c(1:12),c(1:12)]
corrplot(cor.martix)
corrplot(as.matrix(cor.martix), 
         method="circle", 
         #type="upper", 
         #col=colorRampPalette(c("yellow2","goldenrod","green"))(8), 
         col=brewer.pal(n=8, name="PuBuGn"), 
         tl.col="black", 
         tl.srt=45, 
         #diag=FALSE,
         #tl.pos='l',
         #p.mat = as.matrix(correlation_p_values[-c(1:12),c(1:12)]), 
         #sig.level = 0.00000000001, 
         insig = "label_sig")
pdf(file.path(Figure3,"corrplot_immune.pdf"),width =9.82, height = 10.74)
print(corrplot(as.matrix(cor.martix), 
         method="circle", 
         #type="upper", 
         #col=colorRampPalette(c("yellow2","goldenrod","green"))(8), 
         col=brewer.pal(n=8, name="PuBuGn"), 
         tl.col="black", 
         tl.srt=45, 
         #diag=FALSE,
         #tl.pos='l',
         #p.mat = as.matrix(correlation_p_values[-c(1:12),c(1:12)]), 
         #sig.level = 0.00000000001, 
         insig = "label_sig")
)
dev.off()
#######################
cor.martix<-as.matrix(cor.martix)[c(19:23),]
cor.martix<-rbind(cor.martix,cor.martix.non_immune)

pdf(file.path(Figure3,"corrplot_full_immune.pdf"),width =9.82, height = 10.74)
print(corrplot(as.matrix(cor.martix), 
               method="circle", 
               #type="upper", 
               #col=colorRampPalette(c("yellow2","goldenrod","green"))(8), 
               col=brewer.pal(n=8, name="PuBuGn"), 
               tl.col="black", 
               tl.srt=45, 
               #diag=FALSE,
               #tl.pos='l',
               #p.mat = as.matrix(correlation_p_values[-c(1:12),c(1:12)]), 
               #sig.level = 0.00000000001, 
               insig = "label_sig")
)
dev.off()
rownames(cor.martix)
####


cor.martix<-cor.martix[c(1:12),c(1:12)]

pdf(file.path(Figure3,"corrplot_full_st_immune.pdf"),width =7.82, height = 6.74)
print(corrplot(as.matrix(cor.martix), 
               method="circle", 
               type="upper", 
               #col=colorRampPalette(c("yellow2","goldenrod","green"))(8), 
               col=brewer.pal(n=8, name="PuBuGn"), 
               tl.col="black", 
               tl.srt=45, 
               #diag=FALSE,
               #tl.pos='l',
               #p.mat = as.matrix(correlation_p_values[-c(1:12),c(1:12)]), 
               #sig.level = 0.00000000001, 
               insig = "label_sig")
)
dev.off()
