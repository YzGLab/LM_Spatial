######
library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)
library(hexbin)
library(viridisLite)
library(RColorBrewer)
library(cowplot)
stcol<-c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#f5e801', '#08f7f0')
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
dim(st_sce)
sub_sce<-readRDS("data/spatial_combinend_removed_MT_RP.RDS")
##############分别计算不同组织中的差异基因 并合并导出
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
write.csv(sce.markers,file=file.path(Figure3,"sce.markers_LM_vs_CRC_celltypes.csv"))

#############top genes 进行画图处理
sce.cell_type_markers<-split(sce.markers$gene,sce.markers$cluster)###############合并两种差异#
sce.cell_type_markers.sig<-sce.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)# top 10 
sub_sce$st_celltype_loc<-paste0(sub_sce$st_celltype,".",sub_sce$set)



Idents(sub_sce)<-sub_sce$st_celltype_loc
fig4b <-DotPlot(sub_sce,features = unique(sce.cell_type_markers.sig$gene),group.by = "st_celltype_loc") +
  #rotate() +
  scale_color_gradientn(colours=brewer.pal(n=8, name="PuBuGn"), guide = "colourbar") +
  theme_minimal_grid() + theme(axis.text.x=element_text(angle=90, hjust=1))
fig4b
ggsave(file.path(Figure3,"all_LM_CRC_markergenes.pdf"),width =12.99, height =4.23)



DoHeatmap(subset(sub_sce, downsample = 1000),
          size = 2.5,
          features=sce.cell_type_markers.sig$gene,
          slot = 'counts', raster = T,group.colors =stcol) + 
  scale_fill_gradient(
    low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
    #mid = "white",
    high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
    #midpoint = 1,
    guide = "colourbar",
    aesthetics = "fill"
  )
##banes*
#names(sce.cell_type_markers)
#sce.cell_type_markers$LM_Fibroblast[!sce.cell_type_markers$LM_Fibroblast %in% intersect(sce.cell_type_markers$CRC_Fibroblast,sce.cell_type_markers$LM_Fibroblast)]



#####################特异性空间特异性 FIbroblast cell genes 探索 可要可不要。
table(sub_sce$st_celltype)
Idents(sub_sce)<-sub_sce$set
sce.markers3 <- FindAllMarkers(sub_sce[,sub_sce$st_celltype=="Fibroblast"], 
                               only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)

sce.markers3.lm<-sce.markers3[sce.markers3$cluster=="LM" ,] 

write.csv(sce.markers3.lm,file=file.path(Figure3,"sce.markers3.lm_fibroblast.csv"))


sce.markers3.lm<-sce.markers3.lm[sce.markers3.lm$gene %in% sce.markers2[sce.markers2$cluster=="LM_Fibroblast",]$gene,]

sce.markers3.crc<-sce.markers3[sce.markers3$cluster=="CRC" ,] %>% .[.$gene %in% sce.markers1[sce.markers1$cluster=="CRC_Fibroblast",]$gene,]
sce.markers3<-rbind(sce.markers3.crc,sce.markers3.lm)

sce.markers.x.sig<- sce.markers3 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)# top 10 



#DotPlot(immunesce[,immunesce$sub_celltype=="Myeloid"],features = unique(sce.markers.x.sig$gene),group.by = "sub_sub_celltype") + RotatedAxis()+
  #scale_fill_gradientn(colours = viridis(256, option = "D"))# +coord_flip()
fig4a <-DotPlot(sub_sce[,sub_sce$st_celltype=="Fibroblast"],features = unique(sce.markers.x.sig$gene),group.by = "set") +
  #rotate() +
  scale_color_gradientn(colours=brewer.pal(n=8, name="PuBuGn"), guide = "colourbar") +
  theme_minimal_grid() + theme(axis.text.x=element_text(angle=90, hjust=1))
fig4a
ggsave(file.path(Figure4,"mye_markergenes_immune.pdf"),width =6.79, height =3)

subsub_sce<-sub_sce[,sub_sce$st_celltype=="Fibroblast"]
Idents(subsub_sce)<-subsub_sce$set
DoHeatmap(subsub_sce,
          size = 2.5,
          features=sce.markers.x.sig$gene,
          slot = 'counts', raster = T,group.colors =stcol) + 
  scale_fill_gradient(
    low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
    #mid = "white",
    high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
    #midpoint = 1,
    guide = "colourbar",
    aesthetics = "fill"
  )
VolcanoPlot(srt = subsub_sce , group_by = "set")
#####################


library(ggplot2)
library(RColorBrewer)
library(ggrepel)
rm(list = ls())

df <- sce.markers3
## 添加一列logCPM
df$logCPM = runif(179,0,16)
#确定是上调还是下调，用于给图中点上色
df$threshold = factor(ifelse(df$P_value  < 0.05 & abs(df$fd) >= 0.25, ifelse(df$fd >= 0.25 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
df$gene <- row.names(df) #添加一列基因名，以便备注
ggplot(df,aes(x=fd,y= -log10(P_value),size = logCPM,fill = threshold))+
  geom_point(colour = "black", shape = 21, stroke = 0.5)+
  scale_size(limits  = c(2, 16))+ #控制最大气泡和最小气泡，调节气泡相对大小
  scale_fill_manual(values=c("#fe0000","#13fc00","#bdbdbd"))+#确定点的颜色
  geom_text_repel(
    data = df[df$P_value==5.060000e-05,],
    aes(label = gene),
    size = 4.5,
    color = "black",
    segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
  ylab('-log10 (Pvalue)')+#修改y轴名称
  xlab('log2 (FoldChange)')+#修改x轴名称
  geom_vline(xintercept=c(-0.25,0.25),lty=2,col="black",lwd=0.5) +#添加横线|logFoldChange|>0.25
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.5) +#添加竖线padj<0.05
  theme_classic(  # 主题设置，这个是无线条主题
    base_line_size = 1 # 坐标轴的粗细
  )+
  guides(fill=guide_legend(override.aes = list(size=5)))+ #图例圈圈大小
  theme(axis.title.x = element_text(size = 10, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 10,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 7, 
                                   face = "bold")
  )


ggsave("volcano_2.0.pdf", height = 4, width = 5)
