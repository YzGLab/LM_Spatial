#####
sce<-readRDS("data/metabolits_scp_st_sp_sce.RDS")
sample<-"LM_2"
##top metabolits in the different 

cluster_diff <- read_csv("ResultsLM/Figure7_sp/LM_2cluster3.csv")
cluster_diff <- read_csv("ResultsLM/Figure7_sp/LM_2cluster6.csv")
cluster_diff <- read_csv("ResultsLM/Figure7_sp/LM_2cluster8.csv")


colnames( cluster_diff)
cluster_diff.sig <- cluster_diff %>% group_by(class) %>% top_n(n =10, wt = statistic)# top 10
cluster_diff.sig<-cluster_diff.sig[order(cluster_diff.sig$class,decreasing = F),]
cluster_diff.sig$mz<-paste0("mz",cluster_diff.sig$mz)



# 定义样本名称
sample_names <- c("LM_1", "LM_2", "C_1", "C_2")

extract_data_meta<-function(sample) {
  # 读取数据
  print(sample)
  load(paste0("~/Projects/LM_spatial/ResultsLM/Figure7_sp/",sample,"sample_classification.Rdata"))
  sp<-readRDS(paste0("/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/stRNA_scMetabolismuniondata/",sample,"_both.rds"))
  sp_matraix<-sp$metadata
  rownames(sp_matraix)<-sp_matraix$mz
  #############################
  #table(st_sce$orig.ident)
  sp_matraix_filter<-sp_matraix[,c(15:ncol(sp_matraix))]
  colnames(sp_matraix_filter)<-paste0(colnames(sp_matraix_filter),"-1")
  colnames(sp_matraix_filter)<-gsub("\\.","_",colnames(sp_matraix_filter))

return(sp_matraix_filter)
}

rownames(sp_matraix_filter)<-paste0("mz",rownames(sp_matraix_filter))


################
sub_sce<-sce[,sce@meta.data$orig.ident==sample]
meta_sub<-sub_sce@meta.data
meta_sub<-cbind(spotID=rownames(meta_sub),meta_sub)
##############################
expression_meta<-extract_data_meta(sample)
rownames(expression_meta)<-paste0("mz",rownames(expression_meta))
# 去除全为0的代谢物
expression_meta <- expression_meta[ rowSums(expression_meta ) > 0.1,]
###########
expression_meta<-t(expression_meta) %>% cbind(spotID=rownames(.),.)
meta_sub<-merge(meta_sub,expression_meta,by="spotID")

###########################
table(meta_sub$st_celltype)

library(ComplexHeatmap)
rt<-meta_sub


rt$st_celltype <- factor(rt$st_celltype, 
                               levels = c("Fibroblast", "Tumor", "Hepatocytes", "Normal Epi", 
                                          "Monocyte", "B/Plasma", "Lamina propria"))


############
rt$cluster3<-as.factor(rt$cluster3)
rt$cluster6<-as.factor(rt$cluster6)
rt$cluster8<-as.factor(rt$cluster8)
rt$cluster12<-as.factor(rt$cluster12)
# 按照 cluster3 进行排序
rt<-rt[order(rt$st_celltype,decreasing = T),]


rt<-rt[order(rt$cluster3,decreasing = F),]
rt<-rt[order(rt$cluster6,decreasing = F),]
rt<-rt[order(rt$cluster8,decreasing = F),]
rt<-rt[order(rt$cluster12,decreasing = T),]

sorted_expression_meta <- rt[,colnames(rt) %in% cluster_diff.sig$mz]
sorted_expression_meta[]<-apply(sorted_expression_meta,2,as.numeric) 



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
  cluster3= c('1'="red",'2' = '#9467bd','3'='#1f77b4'),
  cluster6= c('1'="red",'2' = '#9467bd','3'='#1f77b4',"4"='#2ca02c',"5"='#f5e801',"6"='#08f7f0'),
  cluster8= c('1'="red",'2' = '#9467bd','3'='#1f77b4',"4"='#2ca02c',"5"='#f5e801',"6"='#08f7f0',"7"= "#C3BC3F","8"="#C7E9C0"),
  #cluster12= c('1'="red",'2' = '#9467bd','3'='#1f77b4'),
  B.Plasma= circlize::colorRamp2(c(0,15),c('white','#d62728' )),
  Endothelial= circlize::colorRamp2(c(0,15),c('white','#d62728' )),
  Epithelial= circlize::colorRamp2(c(0,30),c('white','#d62728' )),
  Fibroblasts= circlize::colorRamp2(c(0,30),c('white','#d62728' )),
  Myeloid= circlize::colorRamp2(c(0,20),c('white','#d62728' )),
  T.ILC= circlize::colorRamp2(c(0,10),c('white','#d62728' )),
  #LC= c("0" = "White", "1" = "black"),
  #BCLC= c("A" = "White", "B+C" = "black"),
  AJCC= c("I+II" = "white", "III+IV" = "black"),
  Grade=c("G1-G2" = "white", "G3-G4" = "black"),
  NAFLD= c("NO" = "white", "YES" = "black"),
  #Enapsulation= c("0" = "White", "1" = "black"),
  Metastasis= c("NO" = "White", "YES" = "black"),
  OS= c("0" = "white", "1" = "black"),
  TMB= c("0" = "white", "1" = "black"),
  age60= c("0" = "white", "1" = "black"),
  Vascularinvasion=c("0" = "white", "1" = "black"),
  AFP= c("0" = "White", "1" = "black"),
  xcell_OS_score_Cat= c("high" = "red", "low" = "blue"),
  HIF_OS_score_Cat= c("high" = "red", "low" ="blue")
)

type1<-rt[,c(22,26:31,42)]
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
type2<-rt[,c(26:31)]
ha= HeatmapAnnotation(df=data.frame(type2),border = TRUE, 
                      #col=col,
                      #gp = gpar(fontsize=8),
                      gap = unit(0, "points"),
                      annotation_name_side ="left")
matrix<-c("XXXX","XXXXXs")
ha33 = rowAnnotation(foo = anno_empty(border = TRUE, 
                                      width = max_text_width(unlist(matrix)) + unit(0, "mm")))
###establish the heatmap
heatcluster<-t(scale(sorted_expression_meta))
ht_opt$message = FALSE
p1<-Heatmap(heatcluster,
            name = "Z-score",
            cluster_columns =F,
            cluster_rows =T,
            #column_order=c("1high","2median","3low"),
            #row_split = splits,
            #column_split = rt$HIFxcellcat,
            #row_labels = splits,
            gap = unit(1, "mm"),
            col = circlize::colorRamp2(c(-2,0,2),c('#1f77b4',"white",'#d62728' )),
            #col = colorRamp2(c(-4,0,4), c("navyblue", "white", "red")),
            show_heatmap_legend = TRUE,
            row_names_gp = gpar(fontsize =6),
            column_names_gp = gpar(fontsize = 8),
            row_names_side = "left",
            #row_names_rot = 30,
            #row_names_max_width = unit(3, "cm"),
            show_column_names = FALSE,
            show_row_names = TRUE,
            #width = unit(120, "mm"),
            #height = unit(140, "mm"),
            top_annotation = ha1,
            #bottom_annotation=ha,
            #right_annotation = ha,
            #left_annotation=ha3,
            row_dend_side = "right",
            show_row_dend=FALSE,
            show_column_dend = FALSE,
            #name = "ht2",
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


pdf(file = file.path(Figure7,"ST_Celltype_sp_heatmap_LM2_cluster3.pdf"),width = 8.5,height =5)
pdf(file = file.path(Figure7,"ST_Celltype_sp_heatmap_LM2_cluster6.pdf"),width = 8.5,height =7)
pdf(file = file.path(Figure7,"ST_Celltype_sp_heatmap_LM2_cluster8.pdf"),width = 8.5,height =9)
draw(p1, newpage = TRUE, 
     column_title = "ST_celltype", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     heatmap_legend_side = "right")
dev.off()

for(i in 1:2) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x =0, width = unit(1, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(matrix[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left")
  })
}

