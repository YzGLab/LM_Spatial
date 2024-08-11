##############
fib_sce<-readRDS("data/fib_sce_scRNA.Rds")
############################
table(fib_sce$sub_celltype_genes)
Idents(fib_sce)<-fib_sce$sub_celltype_genes
#########################################
sce.markers<- FindAllMarkers(fib_sce, 
                             only.pos = TRUE, 
                             min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.x.sig<-sce.markers[sce.markers<0.001,]#########
sce.markers.x.sig<- sce.markers.x.sig %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)# top 10 
sce.markers.x.sig.list<-split(sce.markers.x.sig$gene,sce.markers.x.sig$cluster)

##
Idents(fib_sce)<-fib_sce$set
sce.markers.group<- FindAllMarkers(fib_sce, 
                                   only.pos = TRUE, 
                                   min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.group.sig<-sce.markers.group[sce.markers.group$p_val_adj<0.001,]
sce.markers.group.sig<- sce.markers.group.sig %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)# top 10 


a<-unique(sce.markers.x.sig$gene[sce.markers.x.sig$cluster=="Fib(C4)-CXCL8"])
b<-sce.markers.group.sig[sce.markers.group.sig$cluster=="LM",]$gene
# Remove NAs
markergeness<-intersect(a,b)
#############################epi tumor relation
sce<-readRDS("data/Colon_HC_combined.sce.rds")
load("data/totalscRNAdata.Rdata")
############################
sce<-AddMetaData(sce,metadata)
DefaultAssay(sce)<-'RNA'
# 假设'seurat_obj'是您的Seurat对象
# 首先，找出包含'MT-'或'RP'的基因
mt_genes <- grep("^MT-", rownames(sce), value = TRUE)
rp_genes <- grep("^RP[L|S]", rownames(sce), value = TRUE)
# 合并MT和RP基因列表
genes_to_remove <- union(mt_genes, rp_genes)
# 从Seurat对象中移除这些基因
sc_sce<- subset(sce, features = setdiff(rownames(sce), genes_to_remove))
sce<-NULL
################################################################################
sub_sc_sce<-sc_sce[,sc_sce$sub_celltype=="Epithelial"]

sub_sc_sce$copykay_set<-paste0(sub_sc_sce$set,"_",sub_sc_sce$copykay)
Idents(sub_sc_sce)<-sub_sc_sce$copykay_set
table(sub_sc_sce$copykay_set)
sce.markers.copykay <- FindAllMarkers(sub_sc_sce[,!is.na(sub_sc_sce$copykay)], 
                                      only.pos = TRUE, 
                                      min.pct = 0.25, logfc.threshold = 0.25)
sce.markers.copykay <-sce.markers.copykay[sce.markers.copykay$p_val_adj<0.001,]
sce.markers.copykay.sig<- sce.markers.copykay %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)# top 10 
###################

sce.geneslist<-split(sce.markers.copykay.sig$gene,sce.markers.copykay.sig$cluster)
sce.geneslist<-c(sce.geneslist,sce.markers.x.sig.list)
sce.geneslist[["fib"]]<-markergeness
#################################cor 
#save(sce.cell_type_markers,sce.markers.epi,sce.markers.epi.sig,file='data/spatial_scRNA_differetgenes_clusters.Rdata')
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
#load('data/spatial_scRNA_differetgenes_clusters.Rdata')# from before differen analysis
DefaultAssay(st_sce)<-"SCT"
print(st_sce)

####使用addmodulescore计算空间转录组学分数
st_sce.copy<-AddModuleScore(st_sce,features=sce.geneslist,name=names(sce.geneslist))
###########cor of cell correlation
st_sc_metadata<-st_sce.copy@meta.data
colnames(st_sc_metadata)[c(36:40)]<-names(sce.geneslist)
####
table(st_sc_metadata$st_celltype)
library(corrplot)
st_sc_metadata<-st_sce.copy@meta.data
cor.martix<-cor(st_sc_metadata[st_sc_metadata$st_celltype=="Tumor" &st_sc_metadata$set=="LM" ,c(35:44)])
#cor.martix<-cor.martix[-c(1:12),c(1:12)]
corrplot(cor.martix)


######################################################################################################
library(AUCell)
cells_rankings <- AUCell_buildRankings(st_sce.copy@assays$SCT@data)  # 关键一步

##load gene set, e.g. GSEA lists from BROAD
#c5 <- read.gmt("c5.cc.v7.1.symbols.gmt") ## ALL  GO  tm
cells_AUC <- AUCell_calcAUC(sce.geneslist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
#aucs <- as.numeric(getAUC(cells_AUC[names(sce.markers.epi.list)[4],]))
names(sce.geneslist)
aucs<-getAUC(cells_AUC)
aucs<-t(aucs) %>% data.frame()
st_sce.copy<-AddMetaData(st_sce,aucs)


library(corrplot)
st_sc_metadata<-st_sce.copy@meta.data
#cor.martix<-cor(st_sc_metadata[st_sc_metadata$set=="CRC" ,c(35:44)])
cor.martix<-cor(st_sc_metadata[st_sc_metadata$st_celltype %in% c("Tumor" ,  "Fibroblast","Normal Epi") &st_sc_metadata$set=="CRC" ,c(35:44)])
cor.martix<-cor.martix[-c(1:5),c(1:5)]
corrplot(cor.martix)


##############################
pdf(file.path(Figure4,"corrplot.pdf"),width =7.82, height = 6.74)
print(corrplot(corr=cor.martix,
               #method = "color",
               #order = "hclust",
               tl.col="black",
               addCoef.col = "black",
               number.cex = 0.8,
               col=colorRampPalette(c("blue", "white","#FC4E07"))(50),
)
)
dev.off()
