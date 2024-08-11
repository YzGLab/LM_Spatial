library(Seurat)
library(ggplot2)
library(SCP)
##########split cells in scRNA  
##分别做聚类分析 然后在分组 
## 巨噬细胞，CD8细胞，CD4细胞，NK细胞，Mait细胞等， 一个resolution ==1## Fib  endo, epi. 
####load data
sce<-readRDS("data/Colon_HC_combined.sce.rds")
load("data/totalscRNAdata.Rdata")
sce<-AddMetaData(sce,metadata = metadata)
##
colnames(sce@meta.data)
library(scRNAtoolVis)
# make a polar plot
table(sce$sub_celltype)
Idents(sce)<- sce$sub_celltype
sce.markers <- FindAllMarkers(sce, 
                              #only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)

sce$sub_celltype<-factor(sce$sub_celltype,levels = c("CD4","CD8", "Treg", "NK","MAIT","Myeloid","B","Plasma", "Endothelial", "Epithelial" , "Fibroblasts"))

jjVolcano(diffData = sce.markers,
          log2FC.cutoff=0.25,
          #pvalue.cutoff=0.0001,
          adjustP.cutoff=0.001,
          topGeneN=5,
          #col.type = "adjustP",
          #tile.col = corrplot::COL2('RdBu', 15),
          tile.col=  mycol2,
          #tile.col=c("CD4"="#A65628" ,"CD8"="#984EA3" ,"Treg"="#FDB462","NK"="#BEAED4","MAIT"= "#A6D854","Myeloid"= "#1B9E77",
                   # "B"= "#FFD92F","Plasma"= "#E5C494","Endothelial"= "#8DA0CB","Epithelial"= "#E41A1C","Fibroblasts"= "#FF7F00")
          size  =2.5,
          fontface = 'italic',
          polar = T)
####
library(data.table)
non_immmune_metabolism_exp<-fread("/home/data/gaoyuzhen/Projects/LM_spatial/Compass/LM_CRC_nonimmunecells/reactions.tsv")
#non_immmune_metabolism_exp<-data.frame(non_immmune_metabolism_exp)
immmune_metabolism_exp<-fread("/home/data/gaoyuzhen/Projects/LM_spatial/Compass/ImmunecellCompass/reactions.tsv") 
#immmune_metabolism_exp<-data.frame(immmune_metabolism_exp)
metabolism_exp<-merge(non_immmune_metabolism_exp,immmune_metabolism_exp,by="V1")
metabolism_exp[1:5,1:5]
metabolism_exp<-as.data.frame(metabolism_exp)
rownames.metabolism_exp<-metabolism_exp$V1
metabolism_exp<-metabolism_exp[,-1]
rownames(metabolism_exp)<-rownames.metabolism_exp
metabolism_exp[1:5,1:5]
#colnames(metabolism_exp)<-gsub()
sce@assays$RNA@counts[1:5,1:5]
#sce@assays@metabolism<-metabolism_exp
#metabolism_exp<-metabolism_exp[,colnames(metabolism_exp) %in% rownames(metadata)]


#########
sce<-readRDS("data/full_sce_with_metabolism.RDS")
########
##筛选一下
load("~/Projects/LiverCancer/Compass_envData.Rdata")
reaction<-read.csv("/home/data/gaoyuzhen/Projects/LiverCancer/data/reaction_metadata.csv")
compass_reaction<- reaction  %>% 
  left_join(
    select(compass_reaction_partitions, "reaction_id", "reaction_no_direction"),
    by = "reaction_no_direction"
  ) %>%
# Keep only "confident reactions", as defined in our paper.
filter(EC_number!="N/A") %>%
filter(confidence == "0" | confidence == "4") 
#%>%
# mutate(subsystem = case_when(
#  reaction_no_direction == "SPMDOX" ~ "Arginine and Proline Metabolism",
# subsystem == "Citric acid cycle" & !grepl("[m]", formula, fixed = TRUE) ~ "Other",
#  TRUE ~ subsystem
#))
#genes.sig$noreatction<-str_split(genes.sig$genes,"_") %>% sapply("[[",1)
##############
metabolism_exp<-metabolism_exp[rownames(metabolism_exp) %in% compass_reaction$reaction_id,]
metabolism_assay <- CreateAssayObject(counts = metabolism_exp)
sce[["metabolism"]] <- metabolism_assay
##################
sce$sub_celltype<-factor(sce$sub_celltype,levels = c("CD4","CD8", "Treg", "NK","MAIT","Myeloid","B","Plasma", "Endothelial", "Epithelial" , "Fibroblasts"))
DefaultAssay(sce)<-"metabolism"
Idents(sce)<- sce$sub_celltype

sce.markers.metabolits1 <- FindAllMarkers(sce[,sce$immunecells=="immunecells"],  
                            #only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sce.markers.metabolits1,file = file.path(Figure2,"jVolcano_metabolits_sc_immunecells.csv"))

jjVolcano(diffData = sce.markers.metabolits1,
          log2FC.cutoff=0.25,
          #pvalue.cutoff=0.0001,
          adjustP.cutoff=0.001,
          topGeneN=5,
          #col.type = "adjustP",
          #tile.col = corrplot::COL2('RdBu', 15),
          tile.col=  mycol2,
          #tile.col=c("CD4"="#A65628" ,"CD8"="#984EA3" ,"Treg"="#FDB462","NK"="#BEAED4","MAIT"= "#A6D854","Myeloid"= "#1B9E77",
          # "B"= "#FFD92F","Plasma"= "#E5C494","Endothelial"= "#8DA0CB","Epithelial"= "#E41A1C","Fibroblasts"= "#FF7F00")
          size  =1.5,
          fontface = 'italic',
          polar = T)
ggsave(file = file.path(Figure2,"jjVolcano_metabolits_sc_immunecells.pdf"),width = 11.5,height =11)

sce.markers.metabolits2 <- FindAllMarkers(sce[,sce$immunecells=="nonimmunecells"], 
                                          #only.pos = TRUE, 
                                          min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sce.markers.metabolits2,file = file.path(Figure2,"jVolcano_metabolits_sc_immunecells.csv"))
jjVolcano(diffData = sce.markers.metabolits2,
          #log2FC.cutoff=0.25,
          #pvalue.cutoff=0.0001,
          #adjustP.cutoff=0.001,
          topGeneN=10,
          #col.type = "adjustP",
          #tile.col = corrplot::COL2('RdBu', 15),
          tile.col=  mycol2,
          #tile.col=c("CD4"="#A65628" ,"CD8"="#984EA3" ,"Treg"="#FDB462","NK"="#BEAED4","MAIT"= "#A6D854","Myeloid"= "#1B9E77",
          # "B"= "#FFD92F","Plasma"= "#E5C494","Endothelial"= "#8DA0CB","Epithelial"= "#E41A1C","Fibroblasts"= "#FF7F00")
          size  =2.5,
          fontface = 'italic',
          polar = T)
ggsave(file = file.path(Figure2,"jjVolcano_metabolits_sc_non_immunecells.pdf"),width = 11.5,height = 11)
jjVolcano(diffData =sce.markers.meta,log2FC.cutoff=1,
          #pvalue.cutoff=0.0001,
          tile.col= mycol2,
          adjustP.cutoff=0.001)

############scMETAbolism
library(scMetabolism)
countexp.Seurat<-sc.metabolism.Seurat(obj = sce, method = "VISION", imputation = F, ncores = 60, metabolism.type = "KEGG")
dim(countexp.Seurat)
#countexp.Seurat2<-sc.metabolism.Seurat(obj = sce, method = "VISION", imputation = F, ncores = 70, metabolism.type = "REACTOME")
metabolism_assay_pathway <- countexp.Seurat@assays$METABOLISM$score
colnames(metabolism_assay_pathway)<-gsub("\\.","-",colnames(metabolism_assay_pathway))
write.csv(metabolism_assay_pathway,file = file.path(Figure2,"metabolism_assay_pathway_scrNA.csv"))
metabolism_assay_pathway<-CreateAssayObject(counts = metabolism_assay_pathway)

sce[["metabolism_pathway"]] <- metabolism_assay_pathway
DefaultAssay(sce)<-"metabolism_pathway"
Idents(sce)<- sce$sub_celltype
sce$seurat_clusters<-sce$sub_celltype
sce.markers.meta <- FindAllMarkers(sce, 
                                   #only.pos = TRUE, 
                                   min.pct = 0.15, logfc.threshold = 0.15)
#####################################################################################################
write.csv(sce.markers.meta,file = file.path(Figure2,"jjVolcano_metabolismpathway_sc_immunecells.csv"))
sce.markers.meta.sig <- sce.markers.meta %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)# top 10 
sce.markers.meta.sig<-data.frame(sce.markers.meta.sig)
jjDotPlot(object = sce,
          xtree = T,
          ytree = T,
          rescale = T,
          dot.col = c('blue','white','red'),
          rescale.min = -2,
          rescale.max = 2,
          midpoint = 0,
          gene = unique(sce.markers.meta.sig$gene))
ggsave(file = file.path(Figure2,"jjDotPlot_metabolismpathway_sc_immunecells.pdf"),width = 11.5,height = 11)
################################################################
rownames(sce.markers.meta)<-gsub('metabolism',"M",rownames(sce.markers.meta))
sce.markers.meta$gene<-gsub('metabolism',"M",sce.markers.meta$gene)
jjVolcano(diffData = sce.markers.meta,
          log2FC.cutoff=0.25,
          #pvalue.cutoff=0.0001,
          #adjustP.cutoff=0.001,
          topGeneN=5,
          #col.type = "adjustP",
          #tile.col = corrplot::COL2('RdBu', 15),
          tile.col=  mycol2,
          #tile.col=c("CD4"="#A65628" ,"CD8"="#984EA3" ,"Treg"="#FDB462","NK"="#BEAED4","MAIT"= "#A6D854","Myeloid"= "#1B9E77",
          # "B"= "#FFD92F","Plasma"= "#E5C494","Endothelial"= "#8DA0CB","Epithelial"= "#E41A1C","Fibroblasts"= "#FF7F00")
          size  =1.5,
          fontface = 'italic',
          polar = T)
ggsave(file = file.path(Figure2,"jjVolcano_metabolismpathway_sc_immunecells.pdf"),width = 11.5,height = 11)

CellDimPlot(
  srt = sce, group.by =  "sub_celltype",
  palcolor =c('#1B9E77', '#377EB8', '#984EA3', '#FFFFB3', '#FDBF6F', '#E41A1C'),
  # palcolor =mycol,
  #split.by = "site",
  
  reduction = "UMAP", theme_use = "theme_blank"
)

###########

saveRDS(sce,file = "data/full_sce_with_metabolism.RDS")

#### micropooled
expdata1 <- read.delim2("~/Data/Rdata/pan_sc_met/GPCR/Compass/epi/reactions.tsv", header = T)
# rename micropooled data to "linear_gene_expression_matrix.tsv"

metadata1 <- read.delim2("~/Data/Rdata/pan_sc_met/GPCR/Compass/epi/micropools.tsv", header = T)
metadata2 <- sce_repair@meta.data %>% as.data.frame() %>% rownames_to_column("cell_id") %>% dplyr::select("cell_id", "Cell_type_KEGG")
metadata1 <- metadata1 %>% left_join(metadata2, by = c("X" = "cell_id"))

## 取多数/认定micropool之后的细胞属于哪一类
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

metadata3 <- metadata1 %>% group_by(microcluster) %>% summarise(Cell_type = getmode(Cell_type_KEGG)) 
metadata3$cell_id <- paste0("cluster_", metadata3$microcluster)
colnames(metadata3)
write.csv(metadata3, file="/home/data/yushaobo/Data/Rdata/pan_sc_met/GPCR/Compass/epi/cell_metadata.csv", row.names = F)

#