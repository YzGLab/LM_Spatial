
setwd("/home/data/gaoyuzhen/Projects/LM_spatial")
#####match 代谢物
library(readxl)
Allsample_pos <- read.delim2("spatial_metabolism/Metabolism_martix/Allsample-pos-sample_level.xls")
Allsample_neg <- read.delim2("spatial_metabolism/Metabolism_martix/Allsample-neg-sample_level.xls")
#################
colnames(Allsample_pos)[c(14:17)]<-c("C_1.ALL", "C_2.ALL","LM_1.ALL", "LM_2.ALL")
colnames(Allsample_neg)[c(14:17)]<-c("C_1.ALL", "C_2.ALL","LM_1.ALL", "LM_2.ALL")
# 绘制火山图
names(x)
library(ggplot2)
x<-rbind(Allsample_neg,Allsample_pos)
matched_data.x<-x[x$mz %in% combined_LM$gene,]
# 按照 mz 列进行匹配
combined_LM_up<-combined_LM[combined_LM$avg_log2FC>0,]
combined_LM_down<-combined_LM[combined_LM$avg_log2FC<0,]

matched_data.up<-x[x$mz %in% combined_LM_up$gene,]
matched_data.down<-x[x$mz %in% combined_LM_down$gene,]


#####################speclia in LM




final_results$mz<-final_results$metabolits
matched_data.x<-merge(x,final_results,by="mz")
#matched_data.x<-x[x$mz %in% final_results$metabolits,]

write.csv(final_results,file=file.path(Figure8,"full_st_celltype_fib_final_results_metabolites.csv"))

write.csv(matched_data.x,file=file.path(Figure8,"full_st_celltype_fib.csv"))


VlnPlot(sub_sce,features = unique(all_markers_filtered_combined$gene),pt.size = 0,group.by="st_celltype")
VlnPlot(sce[,sce$st_celltype=="Fibroblast"],features =  unique(all_markers_filtered_combined$gene),pt.size = 0,group.by="st_celltype",split.by = "orig.ident")

VlnPlot(sce[,sce$st_celltype=="Fibroblast"],features =  unique(combined_LM$gene),pt.size = 0,group.by="set")

VlnPlot(sub_sce,features =  sce.markers.sig.LM$gene[1:10],pt.size = 0,group.by="st_celltype")






#matched_data.up

load("data/sp_cardnial_dataset_full_four_samples.Rdata")
library(Cardinal) 
combinedImgSet <- combine(msi_data_list[[1]],msi_data_list[[2]])
image( combinedImgSet,mz=sce.markers.sig[sce.markers.sig$cluster=="Fibroblast",]$gene[1:5],asp=1.5,
       colorscale =magma,
       col = discrete.colors,
       #zlim=c(1,20),
       xlim=c(-10,125))




#########接差异分析
sub_sce<-sce[,sce$orig.ident=="LM_1"]

sub_sce<-sce[,sce$set=="LM"]
sub_sce<-sce[,sce$st_celltype=="Fibroblast"]
rt<-sub_sce@assays$metabolits@data %>% data.frame(.,check.names = F)
#rt<-rt[rownames(rt) %in% sce.markers.sig[sce.markers.sig$cluster=="Fibroblast",]$gene,]
#rt<-rt[rownames(rt) %in% combined_LM_down$gene[combined_LM_down$sample_name=="LM_1"],]
rt<-rt[rownames(rt) %in% unique(matched_data.x$mz),]
#rt<-t(scale(t(rt)))
rt <- rt[ rowSums(rt ) > 0.1,]
rt<-t(rt) %>% cbind(spotID=rownames(.),.)
###########
colnames(meta_sub)
#########
meta_sub<-sub_sce@meta.data
meta_sub<- cbind(spotID=rownames(meta_sub),meta_sub)
meta_sub<-merge(meta_sub,rt,by="spotID")
rownames(meta_sub)<-meta_sub$spotID

library(ComplexHeatmap)
rt<-meta_sub
rt<-rt[order(rt$st_celltype,decreasing = T),]
col = list(
  st_celltype= c(
    'Tumor' = '#ff7f0e', 
    'Normal Epi' = '#1f77b4', 
    'Hepatocytes' = '#d62728', 
    'Fibroblast' = '#2ca02c', 
    'B/Plasma' = '#f5e801',
    'Monocyte' = '#08f7f0',
    'Lamina propria' = '#9467bd'
  ),
  orig.ident=c("LM_1"="white","LM_2"="black","C_1"="white","C_2"="black"),
  #cluster3= c('1'="red",'2' = '#9467bd','3'='#1f77b4'),
  #cluster6= c('1'="red",'2' = '#9467bd','3'='#1f77b4',"4"='#2ca02c',"5"='#f5e801',"6"='#08f7f0'),
  #cluster8= c('1'="red",'2' = '#9467bd','3'='#1f77b4',"4"='#2ca02c',"5"='#f5e801',"6"='#08f7f0',"7"= "#C3BC3F","8"="#C7E9C0"),
  #cluster12= c('1'="red",'2' = '#9467bd','3'='#1f77b4'),
  metabolic_activity= circlize::colorRamp2(c(0,1),c('white','#d62728' )),
  B.Plasma= circlize::colorRamp2(c(0,15),c('white','#d62728' )),
  Endothelial= circlize::colorRamp2(c(0,15),c('white','#d62728' )),
  Epithelial= circlize::colorRamp2(c(0,30),c('white','#d62728' )),
  Fibroblasts= circlize::colorRamp2(c(0,30),c('white','#d62728' )),
  Myeloid= circlize::colorRamp2(c(0,20),c('white','#d62728' )),
  T.ILC= circlize::colorRamp2(c(0,10),c('white','#d62728' ))
)



type1<-rt[,c(22,2,26:31,36)]##cluster3
ha1= HeatmapAnnotation(df=data.frame(type1),border = TRUE, 
                       col=col,
                       #height = unit(30, "mm"),
                       gap=unit(0.4, "mm"),
                       #gp = gpar(fontsize=8),
                       #annotation_width =unit(4, "cm"),
                       annotation_name_side ="left"
                       #annotation_height = unit(3, "cm")
                       #width = 1
)

sorted_rt<-rt[,c(45:ncol(rt))]
sorted_rt[]<-apply(sorted_rt,2,as.numeric)
###establish the heatmap
heatcluster<-t(scale(sorted_rt))
ht_opt$message = FALSE
p1<-Heatmap(heatcluster,
            name = "Z-score",
            cluster_columns =T,
            cluster_rows =T,
            #column_order=c("1high","2median","3low"),
            #row_split = rt$st_celltype,
            column_split = rt$set,
            #row_labels = splits,
            gap = unit(1, "mm"),
            col = circlize::colorRamp2(c(-1,0,1),c('#1f77b4',"white",'#d62728' )),
            #col = colorRamp2(c(-4,0,4), c("navyblue", "white", "red")),
            show_heatmap_legend = TRUE,
            row_names_gp = gpar(fontsize =8),
            column_names_gp = gpar(fontsize = 8),
            row_names_side = "left",       
            show_column_names = FALSE,
            show_row_names = TRUE,
            top_annotation = ha1,
            row_dend_side = "right",
            show_row_dend=FALSE,
            show_column_dend = FALSE,
            border = TRUE,
            na_col = "grey",
            column_title_gp = gpar(fill = "WHITE", col = "BLACK", border = "WHITE"),
            row_title_gp = gpar(fill = "WHITE", col = "BLACK", border = "WHITE"),
            row_title = "mz",
            #column_title = "Combinded TME with HIF "
)

draw(p1, newpage = TRUE, 
     column_title = "ST_celltype", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     heatmap_legend_side = "right")

pdf(file = file.path(Figure8,"ST_Celltype_fib_heatmap_LM_diff_metabolites.pdf"),width = 8.5,height =7)
draw(p1, newpage = TRUE, 
     column_title = "ST_celltype", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     heatmap_legend_side = "right")
dev.off()
