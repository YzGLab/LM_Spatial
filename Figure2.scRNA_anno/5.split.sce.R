#
##########split cells in scRNA  
##分别做聚类分析 然后在分组 
## 巨噬细胞，CD8细胞，CD4细胞，NK细胞，Mait细胞等， 一个resolution ==1## Fib  endo, epi. 

##T cell gene annotation 
#split sce
table(sce$Major_celltype,sce$sub_celltype)
sce$Major_celltype[sce$sub_celltype %in% c("CD4",  "CD8", 'Treg', 'NK', 'MAIT')]="T/ILC"
sce$Major_celltype[sce$sub_celltype %in% c("Myeloid")]="Myeloid"
sce$Major_celltype[sce$sub_celltype %in% c("B", 'Plasma')]="B/Plasma"
sce$immunecells[sce$Major_celltype %in% c("Myeloid","B/Plasma","T/ILC")]="immunecells"
sce$immunecells[!sce$Major_celltype %in% c("Myeloid","B/Plasma","T/ILC")]="nonimmunecells"
table(sce$immunecells)
immunesce<-sce[,sce$immunecells=="immunecells"]
nonimmunesce<-sce[,sce$immunecells=="nonimmunecells"]
saveRDS(sce, file = "data/Colon_HC_combined.sce.rds")
saveRDS(immunesce, file = "data/immunesce.rds")
saveRDS(nonimmunesce, file = "data/nonimmunesce.rds")
##################################
immunescemetadata<-immunesce@meta.data
nonimmunescemetadata<-nonimmunesce@meta.data
colnames(nonimmunescemetadata)
immunescemetadata<-immunescemetadata[,-c(13,14)]
nonimmunescemetadata<-nonimmunescemetadata[,-c(13,14)]
totalmetadata<-rbind(immunescemetadata,nonimmunescemetadata)
library(Seurat)
sce<-AddMetaData(sce,totalmetadata)
metadata<-sce@meta.data
save(metadata,file = "data/totalscRNAdata.Rdata")
##############################################################################
