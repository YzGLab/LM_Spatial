###
fib_sce<-readRDS("data/fib_sce_scRNA.Rds")
##
library(SCENIC)
# https://resources.aertslab.org/cistarget/
db='/home/data/gaoyuzhen/Projects/LiverCancer/data/cisTarget_database'
list.files(db)
#data(list="motifAnnotations_hgnc", package="RcisTarget")
# 保证cisTarget_databases 文件夹下面有下载好 的文件
motifAnnotations_hgnc <- motifAnnotations
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir=db , nCores=16) 
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#saveRDS(cellInfo, file="int/cellInfo.Rds")
data.input<-fib_sce@assays$RNA@data
#data.input<-data.input[,colnames(data.input) %in% fib.metadatas$Cell]
##
exprMat<-as.matrix(data.input) 
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered<-t(exprMat_filtered)
write.csv(exprMat_filtered,file="fib_cells.csv")


rm(list = ls())  
library(SCENIC)
packageVersion("SCENIC")  
library(SCopeLoomR)
scenicLoomPath='/home/data/gaoyuzhen/Projects/LM_spatial/data/cisTarget_databases/sample_SCENIC_fib.loom'
loom <- open_loom(scenicLoomPath)
# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")##
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
###########################


#run Tcell SMC regulation
#sce<-readRDS("fibroblast_sce_annotation.RDS")
sce<-fib_sce
library(pheatmap) 
n=t(scale(t(getAUC(regulonAUC[,] )))) # 'scale'可以对log-ratio数值进行归一化
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
dim(n) 
ac=data.frame(group= as.character(Idents(sce)))
n[1:4,1:4]
n=n[,colnames(n) %in% colnames(sce)]
rownames(ac)=colnames(n) ##
lengths(regulons)
library(ggpubr)
#aucs<-getAUC(regulonAUC)
aucs<-getAUC(regulonAUC)
aucs<-t(aucs)
regulonss<-regulons
for(i in c(1:369)){
  names(regulonss)[i]<-gsub("\\+",paste0(lengths(regulons)[i],"g"),names(regulons)[i])
}
names(regulonss)
colnames(aucs)<-names(regulonss)
#aucs<-t(n)
#aucs<-data.frame(aucs)
Idents(sce)<-sce$sub_celltype_genes
aucs<-cbind(aucs,group= as.character(Idents(sce)))
dim(aucs)
#write.csv(aucs,file="Bcell_TF_idente.csv")
#aucs<-read.csv("Bcell_TF_idente.csv",check.names = F)
head(aucs)
#aucs<-data.frame(aucs)
#group_by(aucs, group) %>% summarize_each(funs(mean), colnames(aucs)[1:81])
aucss<-apply(aucs,2,as.numeric)
rownames(aucss)<-rownames(aucs)
aucsmenas<-aggregate(aucss[,c(1:369)],list(aucs[,370]),mean)#########根据最后确定的TF数量进行构建 heatmap
aucsmenas<-t(aucsmenas)
colnames(aucsmenas)<-aucsmenas[1,]
aucsmenas<-aucsmenas[-c(1),]
aucsmenas[1:5,1:4]
######
data<-aucs
outTab<-c()
for (i in colnames(aucs)[c(1:370)]){
  geneName=i
  data[,i]<-as.numeric(data[,i])
  rt=rbind(expression=data[,i],group=data[,"group"])
  rt=as.matrix(t(rt))
  rt<-data.frame(rt)
  #Means=rowMeans(data[,i])
  #if(is.numeric(Means) & Means>0){
  kTest<-kruskal.test(expression ~ group, data=rt)
  pvalue=kTest$p.value
  outTab=rbind(outTab,cbind(gene=i,pValue=pvalue))
  message(i,"is done now!!!")
  #}
}#
outTab<-data.frame(outTab)
#gene.cox.i[,11] <- p.adjust(gene.cox.i[,9], method = "fdr")
outTab$fdr<- p.adjust(outTab$pValue, method = "fdr")
outTabs<-outTab[outTab$pValue<0.05,]
TF_heatmap<-aucsmenas[outTabs$gene,] 
TF_heatmap<-`rownames<-`(apply(TF_heatmap,2,as.numeric), rownames(TF_heatmap))  

pdf(file.path(Figure4,"TF_scenic_fib.pdf"),height = 7,width=4)
pheatmap(TF_heatmap,
         scale = "row",
         show_colnames =T,
         show_rownames = T,
         cluster_cols = F,
         fontsize = 10,
         angle_col = c('45'),#'','#d62728'
         color=colorRampPalette(c('#1f77b4',"white",'#E41A1C'))(20)
)
dev.off()
#c('#984EA3','#377EB8','#1B9E77','#E41A1C')

