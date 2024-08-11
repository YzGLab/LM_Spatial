
############
stcol<-c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#f5e801', '#08f7f0')

sample<-"LM_2" 
sample<-"LM_1" 
sample<-"C_1" 
sample<-"C_2" 
###########
library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)
##
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
##############
#list.files("/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/stRNA_scMetabolismuniondata")

load(paste0("~/Projects/LM_spatial/ResultsLM/Figure7_sp/",sample,"sample_classification.Rdata"))

sp<-readRDS(paste0("/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/stRNA_scMetabolismuniondata/",sample,"_both.rds"))
##########################
#匹配数据 
########################################################################################################
sp_matraix<-sp$metadata
rownames(sp_matraix)<-sp_matraix$mz
#############################
table(st_sce$orig.ident)
sp_matraix_filter<-sp_matraix[,c(15:ncol(sp_matraix))]
colnames(sp_matraix_filter)<-paste0(colnames(sp_matraix_filter),"-1")
colnames(sp_matraix_filter)<-gsub("\\.","_",colnames(sp_matraix_filter))

###############################################导入空间代谢组学分群数据
sp_cluster<-data.frame(sample_classification)
#sp_cluster<-cbind(sample=sp$samplename$transname,clusters)
sp_cluster$Sample<-paste0(sp_cluster$Sample,"-1")
sp_cluster$Sample<-gsub("\\.","_",sp_cluster$Sample)
rownames(sp_cluster)<-sp_cluster$Sample
colnames(sp_cluster)[c(2:ncol(sp_cluster))]<-c(3,6,8,12)
colnames(sp_cluster)[c(2:ncol(sp_cluster))]<-paste0("cluster",colnames(sp_cluster)[c(2:ncol(sp_cluster))])

#################匹配样本 
sub_sce<-st_sce[,st_sce$orig.ident==sample]
sub_sce<-sub_sce[,colnames(sub_sce) %in% sp_cluster$Sample]

#############建立差异分析的assay 
sp_matraix_filter<-sp_matraix_filter[,colnames(sp_matraix_filter) %in% sp_cluster$Sample]
sp_matraix_assay <- CreateAssayObject(counts =sp_matraix_filter)
sub_sce[["metabolits"]] <- sp_matraix_assay
##########添加信息
sub_sce<-AddMetaData(sub_sce,metadata = sp_cluster)
#####################################################################################
#代谢分群之间相关性（代谢物均值处理，做person 相关性）

####################
#resutl #柱状图
library(janitor)
library(ggplot2)
library(scales)
metadata<-sub_sce@meta.data
colnames(metadata)
##"cluster3"             "cluster6"             "cluster8"             "cluster12" 
##
names(stcol1)<- c('Normal Epi', 'Tumor', 'Fibroblast', 'Hepatocytes', 'Lamina propria', 'B/Plasma', 'Monocyte')
stcol1 <- c('Normal Epi' = '#1f77b4', 
            'Tumor' = '#ff7f0e', 
            'Fibroblast' = '#2ca02c', 
            'Hepatocytes' = '#d62728', 
            'Lamina propria' = '#9467bd', 
            'B/Plasma' = '#f5e801', 
            'Monocyte' = '#08f7f0')

metadata$cluster3<-factor(metadata$cluster3)
p3<-ggplot(data =metadata, aes(x=cluster3, fill = st_celltype)) +
  geom_bar(position="fill")+
  theme_bw()+
  scale_fill_manual(values=stcol1) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="st_celltype")+
  #geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], 
  # label=paste(..count..,paste("(",percent(..count../tapply(..count.., ..x.. ,sum)[..x..]),")",sep=""),sep=" ")),
  # stat="count",color="white", position=position_fill(0.5), vjust=0.5,hjust=0.5,angle = 60)+
  scale_y_continuous(position = "left",labels = percent)+  ####??????????y??????????
  theme(axis.text.y = element_text(size=6, colour = "black"))+
  theme(axis.text.x = element_text(size=6, colour = "black"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 6),
        axis.text.y = element_text(angle = 00, size = 6))+
  theme(legend.position = "left")

p3###

metadata$cluster6<-factor(metadata$cluster6)
p6<-ggplot(data =metadata, aes(x=cluster6, fill = st_celltype)) +
  geom_bar(position="fill")+
  theme_bw()+
  scale_fill_manual(values=stcol1) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="st_celltype")+
  #geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], 
  # label=paste(..count..,paste("(",percent(..count../tapply(..count.., ..x.. ,sum)[..x..]),")",sep=""),sep=" ")),
  # stat="count",color="white", position=position_fill(0.5), vjust=0.5,hjust=0.5,angle = 60)+
  scale_y_continuous(position = "left",labels = percent)+  ####??????????y??????????
  theme(axis.text.y = element_text(size=6, colour = "black"))+
  theme(axis.text.x = element_text(size=6, colour = "black"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 6),
        axis.text.y = element_text(angle = 00, size = 6))+
  theme(legend.position = "left")

p6###
metadata$cluster8<-factor(metadata$cluster8)
p8<-ggplot(data =metadata, aes(x=cluster8, fill = st_celltype)) +
  geom_bar(position="fill")+
  theme_bw()+
  scale_fill_manual(values=stcol1) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="st_celltype")+
  #geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], 
  # label=paste(..count..,paste("(",percent(..count../tapply(..count.., ..x.. ,sum)[..x..]),")",sep=""),sep=" ")),
  # stat="count",color="white", position=position_fill(0.5), vjust=0.5,hjust=0.5,angle = 60)+
  scale_y_continuous(position = "left",labels = percent)+  ####??????????y??????????
  theme(axis.text.y = element_text(size=6, colour = "black"))+
  theme(axis.text.x = element_text(size=6, colour = "black"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 6),
        axis.text.y = element_text(angle = 00, size = 6))+
  theme(legend.position = "left")

p8###

metadata$cluster12<-factor(metadata$cluster12)
p12<-ggplot(data =metadata, aes(x=cluster12, fill = st_celltype)) +
  geom_bar(position="fill")+
  theme_bw()+
  scale_fill_manual(values=stcol1) +
  theme(panel.grid =element_blank()) +
  labs(x="",y="st_celltype")+
  #geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], 
  # label=paste(..count..,paste("(",percent(..count../tapply(..count.., ..x.. ,sum)[..x..]),")",sep=""),sep=" ")),
  # stat="count",color="white", position=position_fill(0.5), vjust=0.5,hjust=0.5,angle = 60)+
  scale_y_continuous(position = "left",labels = percent)+  ####??????????y??????????
  theme(axis.text.y = element_text(size=6, colour = "black"))+
  theme(axis.text.x = element_text(size=6, colour = "black"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 6),
        axis.text.y = element_text(angle = 00, size = 6))+
  theme(legend.position = "left")

p12###
p3+p6+p8+p12

ggsave(file.path(Figure7,paste0(sample,"_sp_cluster_st_boxplot.pdf")),width =10,height = 11 )



