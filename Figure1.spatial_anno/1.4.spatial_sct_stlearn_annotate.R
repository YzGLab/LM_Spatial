###diff for cluster of st sce
library(SCP)
library(Seurat)

load("data/stmetadata.Rdata")
st_sce<-readRDS("data/Colon_HC_spatial.RDS")

DefaultAssay(st_sce)<-"SCT"
#######################################
colnames(stmetadata)
#stcol<-c('#1f77b4','#f87f13','#359c62','#d32929','#69308e','#8c564c','#f33ca9','#f5e801','#08f7f0','#afc7e6')

stcol<-c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#f5e801', '#08f7f0')
DimPlot(st_sce,reduction = "umap",group.by = "leiden",label=T,split.by="orig.ident",cols = stcol)

st_sce$leiden<-factor(st_sce$leiden)
CellDimPlot(
  srt =st_sce, group.by =  "leiden",  stat_plot_size = 0.1,split.by = "orig.ident",ncol = 2,
  reduction = "UMAP", theme_use = "theme_blank",show_stat = T,
  palcolor = stcol,
)
ggsave(file = file.path(Figure2,"st_sce_cluster_CRC_LM_total.pdf"),width = 9.5,height = 9.03)####
######################################
#DimPlot(st_sce, reduction = "umap", label = TRUE) + NoLegend()
#SpatialDimPlot(st_sce, label = TRUE, label.size = 3) + NoLegend()
#cluster10_markers <- FindMarkers(Colon_obj, ident.1 = 1, min.pct = 0.25)
Idents(st_sce)<-st_sce$leiden
sce.markers <- FindAllMarkers(st_sce, group.by="leiden",
                              only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
save(sce.markers,file="data/sce.markers_st.Rdata")
head(sce.markers, n = 5)
st_sce[["percent.mt"]] <- PercentageFeatureSet(st_sce,pattern = "^MT-")
st_sce[["percent.rp"]] <- PercentageFeatureSet(st_sce,pattern = "^RP-")



library(SingleR)
HumanRNA<-HumanPrimaryCellAtlasData()
#HumanRNA<-DatabaseImmuneCellExpressionData()
sce_for_SingleR <- as.SingleCellExperiment(st_sce)
cluster_sce<-st_sce$leiden
#humanImmu <- ImmGenData()
#pred.humanImmu <- SingleR(test = sce_for_SingleR, ref = humanImmu, labels =humanImmu$label.main)
#mouseRNA <- MouseRNAseqData()
pred.sce <- SingleR(test = sce_for_SingleR, ref = HumanRNA, labels =HumanRNA$label.main,clusters=factor(cluster_sce))
#pred.humanRNA <- SingleR(test = sce_for_SingleR, ref = list(BP=Blue.se, HPCA=hpca.se), labels = list(Blue.se$label.main, hpca.se$label.main)) 
table(rownames(pred.sce),pred.sce$labels )

##################
table(st_sce$leiden)
stmetadata$st_celltype<-stmetadata$leiden
stmetadata$st_celltype[stmetadata$leiden=="0"]="Normal Epi"# 
stmetadata$st_celltype[stmetadata$leiden=="1"]="Tumor"# ok
stmetadata$st_celltype[stmetadata$leiden=="2"]="Fibroblast"# OK
stmetadata$st_celltype[stmetadata$leiden=="3"]="Hepatocytes"#OK
stmetadata$st_celltype[stmetadata$leiden=="4"]="Lamina propria"##OK
stmetadata$st_celltype[stmetadata$leiden=="5"]="B/Plasma"
stmetadata$st_celltype[stmetadata$leiden=="6"]="Monocyte"# "ADM"

table(stmetadata$st_celltype)

st_sce<-AddMetaData(st_sce,stmetadata)
table(st_sce$orig.ident)

CellDimPlot(
  srt =st_sce, group.by =  "st_celltype", split.by = "orig.ident", stat_plot_size = 0.1,ncol = 2,
  reduction = "UMAP", theme_use = "theme_blank",show_stat = T,
  palcolor = stcol,
)
st_sce$set<-ifelse(st_sce$orig.ident %like% "C_","CRC","LM")
table(st_sce$st_celltype)

st_sce$st_celltype<-factor(st_sce$st_celltype,levels =c('Normal Epi','Tumor', "Fibroblast","Hepatocytes","Lamina propria","B/Plasma",'Monocyte' ) )
CellDimPlot(
  srt =st_sce, group.by =  "st_celltype", split.by = "set", stat_plot_size = 0.1,ncol = 2,
  reduction = "UMAP", theme_use = "theme_blank",show_stat = T,
  palcolor = stcol,
)
ggsave(file = file.path(Figure2,"st_sce_cluster_CRC_LM_set.pdf"),width = 9.5,height = 9.03)####

##################################

CellStatPlot(st_sce,stat.by ="st_celltype", group.by =c("orig.ident","set"),
             plot_type = "dot",palcolor = stcol)
ggsave(file = file.path(Figure2,"celltype_immune_CRC_LM_Patients_probability.pdf"),width = 9.8,height = 5)

###################
table(st_sce$st_celltype,st_sce$orig.ident)

# 载入必要的库
library(ggplot2)
library(tidyr)
library(dplyr)
# 假设您的原始数据框如下：
data <- data.frame(
  Cell_Type = c('Normal Epi', 'Tumor', 'Fibroblast', 'Hepatocytes', 'Lamina propria', 'B/Plasma', 'Monocyte'),
  C_1 = c(858, 0, 232, 0, 0, 28, 7),
  C_2 = c(1235, 1218, 721, 0, 1249, 105, 25),
  LM_1 = c(1333, 1052, 550, 382, 2, 39, 15),
  LM_2 = c(644, 1041, 517, 1217, 5, 27, 25)
)
# 转换成长格式
long_data <- pivot_longer(data, cols = -Cell_Type, names_to = "Sample", values_to = "Count")

# 计算每个样本的总数
sample_totals <- long_data %>% group_by(Sample) %>% summarise(Total = sum(Count))

# 计算每种细胞类型在每个样本中的比例
long_data <- long_data %>% left_join(sample_totals, by = "Sample") %>% 
  mutate(Proportion = Count / Total)

long_data<-data.frame(long_data)
# 绘制每种细胞类型在四个样本中的比例的小图
long_data$Cell_Type<-factor(long_data$Cell_Type,levels = c('Normal Epi','Tumor', "Fibroblast","Hepatocytes","Lamina propria","B/Plasma",'Monocyte' ))

ggplot(long_data, aes(x = Sample, y = Proportion, group = Cell_Type)) +
  #geom_ribbon(aes(ymin = Proportion - sd, ymax = Proportion + sd, fill = Cell_Type), alpha = 0.2) +
  geom_line(aes(color = Cell_Type), size = 1) +  # 增加线条粗细
  geom_point(aes(color = Cell_Type), size = 3) +  # 增加点的粗细
  scale_fill_manual(values=stcol) +
  scale_color_manual(values=stcol) +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  labs(title = "Cell Type Proportions in Different Samples", x = "Sample", y = "Proportion") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent", colour =stcol,), # 设置面板背景颜色为浅灰色
    plot.background = element_rect(fill ="transparent", colour = stcol), # 设置整个绘图区域的背景颜色为白色
    legend.position = "none"
  )
ggsave(file = file.path(Figure2,"st_sce_cluster_CRC_LM_patients_prop.pdf"),width = 9.5,height = 6.9)####
