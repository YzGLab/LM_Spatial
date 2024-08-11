library(SCP)
library(data.table)
library(stringr)## data from ptyhon SME method for merge four samples
####################################
full_SME<-read.csv("/home/data/gaoyuzhen/Projects/LM_spatial/data/stlearn_results/full_SME_results.csv")
full_SME$rownames<-paste0(str_split(full_SME$X,"-") %>% sapply("[[",1),"-",str_split(full_SME$X,"-") %>% sapply("[[",2))
full_SME$samplenames<-str_split(full_SME$X,"-") %>% sapply("[[",3)
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
stmetadata<-st_sce@meta.data
stmetadata<-cbind(stmetadata,full_SME)
st_sce<-AddMetaData(st_sce,stmetadata)


##########################
DefaultAssay(st_sce)<-"SCT"
#sce.big <- ScaleData(sce.big,verbose = FALSE) %>% RunPCA(pc.gene = sce.big@var.genes,npcs = 30,verbose = FALSE)  %>% RunHarmony(group.by.vars = "orig.ident")
st_sce<- FindNeighbors(st_sce, reduction = "harmony", dims = 1:30)
st_sce <- FindClusters(st_sce, resolution = 0.5)
DimPlot(st_sce,reduction = "umap",label=T,split.by="orig.ident")

#################feature plot for 
colnames(stmetadata)
mycol2<-c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F", "#E41A1C", '#FFD92F', "#C3BC3F", "#FF7F00", "#1B9E77", "#66A61E", "#F1788D")
mycol2<-c( '#A65628', '#984EA3','#FDB462','#BEAED4','#A6D854','#1B9E77','#FFD92F','#E5C494','#8DA0CB', '#E41A1C','#FF7F00')
stcol<-c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#f5e801', '#08f7f0')


DimPlot(st_sce,reduction = "umap",group.by = "leiden",label=T,split.by="orig.ident",cols = stcol)

st_sce$leiden<-factor(st_sce$leiden)
CellDimPlot(
  srt =st_sce, group.by =  "leiden",  stat_plot_size = 0.1,split.by = "orig.ident",ncol = 2,
  reduction = "UMAP", theme_use = "theme_blank",show_stat = T,
  palcolor = stcol,
)
ggsave(file = file.path(Figure1,"st_sce_cluster_CRC_LM_total.pdf"),width = 9.5,height = 9.03)####


