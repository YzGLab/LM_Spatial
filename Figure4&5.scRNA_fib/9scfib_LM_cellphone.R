###CCI
######################
load(file.path(Figure4,"Fibroblasts_markergenes_split_markers_sig.Rdata"))
fib_sce<-readRDS("data/fib_sce_scRNA.Rds")
########################################
scRNA_sce<-readRDS("data/Colon_HC_combined.sce.rds")
data.input<-as.matrix(scRNA_sce@assays$RNA@data)
data.input[1:4,1:5]


metadata<-scRNA_sce@meta.data
#scRNA_sce<-AddMetaData(scRNA_sce,metadata)
#saveRDS(scRNA_sce,file="data/Colon_HC_combined.sce.rds")
#table(metadata$set)
#metadata$set<-ifelse(metadata$orig.ident %in%c("chang1","chang2"),"CRC","LM")
metadata<-metadata[metadata$set=="LM",]



fib_metadata<-fib_sce@meta.data
fib_metadata<-fib_metadata[fib_metadata$set=="LM",]
table(fib_metadata$sub_celltype)
fib_metadata$sub_celltype<-fib_metadata$sub_celltype_genes

metadata<-rbind(metadata[!metadata$sub_celltype=="Fibroblasts",],fib_metadata[,c(1:12)])
table(metadata$sub_celltype)


#write.table(data.input, 'cellphonedb_count.txt', sep='\t', quote=F)
#write.table(data.input, 'cellphonedb_count.txt', sep='\t', quote=F)
#fwrite(data.input, "cellphonedb_count.txt", row.names=F, sep='\t')
########### 找meta.data

meta_fib<-cbind(Cell=rownames(metadata),cellphone_subtype=metadata$sub_celltype)
###

#VlnPlot(scRNA_sce,features = c("SPP1","CD44"),group.by = "sub_celltype",split.by = "set")

#################
colnames(meta_fib)<-c("Cell","cell_type")
meta_fib<-data.frame(meta_fib)
meta_fib<-meta_fib[meta_fib$cell_type!="NA",]
############# 肝转移


data.input_Tumor<-data.input[,colnames(data.input) %in% meta_fib$Cell]
meta_fib<-meta_fib[meta_fib$Cell %in% colnames(data.input_Tumor),]
table(meta_fib$cell_type)

meta_fib<-as.matrix(meta_fib)

write.table(data.input_Tumor, 'cellphonedb_count_Tumor.txt', sep='\t', quote=F)
write.table(meta_fib, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
#data.table::fwrite(data.input_Tumor, "cellphonedb_count_Tumor.txt", row.names=F, sep='\t')
#data.table::fwrite(meta_Tumor, 'cellphonedb_meta.txt', row.names=F, sep='\t')

##########
file_count = "./cellphonedb_count_Tumor.txt"
file_anno = "./cellphonedb_meta.txt"
outdir ="./out"
#projectname="Tumor_mye"
metapath="/home/data/gaoyuzhen/Projects/LM_spatial/cellphonedb_meta.txt"

if (!dir.exists(outdir)){
  dir.create(outdir)
}
system(paste0("docker run --rm --name cellphonedb21 -i -v `pwd`/:/root ydevs/cellphonedb:2.1 cellphonedb method statistical_analysis ", 
              file_anno, " ", file_count, " --counts-data=gene_name", " --iterations=3 --threads=50"))

system(paste0("docker run --rm --name cellphonedb21 -i -v `pwd`/:/root ydevs/cellphonedb:2.1 cellphonedb plot dot_plot --means-path ",
              outdir,'/means.txt',
              " --pvalues-path ",outdir,'/pvalues.txt', " --output-path ",outdir))

system(paste0("docker run --rm --name cellphonedb21 -i -v `pwd`/:/root ydevs/cellphonedb:2.1 cellphonedb plot heatmap_plot", 
              " cellphonedb_meta.txt")) 


######


##cellphone visulization 
#########################
###3.1  cellphone DB for mac metasubtype. 
library(CellChat)
library(tidyr)
library(igraph)
library(pheatmap)
df.net <- read.table("/home/data/gaoyuzhen/Projects/LM_spatial/out/count_network.txt",
                     header = T,sep = "\t",stringsAsFactors = F)
meta.data <- read.table("/home/data/gaoyuzhen/Projects/LM_spatial/cellphonedb_meta.txt",
                        header = T,sep = "\t",stringsAsFactors = F)

groupSize <- as.numeric(table(meta.data$cell_type))

df.net <- spread(df.net, TARGET, count)
rownames(df.net) <- df.net$SOURCE
df.net <- df.net[, -1]
df.net <- as.matrix(df.net)

netVisual_circle(df.net, vertex.weight = groupSize,
                 shape="rectangle",
                 #layout = layout_on_sphere(),
                 weight.scale = T, label.edge= T,
                 title.name = "Number of interactions")
mat <- df.net
#/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure4_fib/fib_LM_cellphone
pdf(file.path(Figure4,paste0("fib_LM_cellphone/Tumor_fib","network_heatmap.pdf")))
pheatmap(mat, cluster_rows = F, 
         cluster_cols = F, color=colorRampPalette(c('#1B9E77',"white",'#E41A1C'))(50))
dev.off()


par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) { 
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)) 
  mat2[i, ] <- mat[i, ]  
  pdf(file.path(Figure4,paste0("fib_LM_cellphone/",rownames(mat)[i],"Tumor_fib_network.pdf")))
  netVisual_circle(mat2, vertex.weight = groupSize, label.edge= T, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()}

mycol2<-c('#1B9E77','#984EA3','#1B9E77','#377EB8','#E41A1C','#B3B3B3','#FB8072','#BEAED4','#A6D854','#E41A1C','#984EA3','#F0027F','#FFD92F','#E5C494','#66C2A5', '#FC8D62','#8DA0CB','#E78AC3','#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999',
          '#8DD3C7', '#FFFFB3','#BEBADA', '#80B1D3','#FDB462','#B3DE69', '#FCCDE5')
mat
pdf(file.path(Figure4,paste0("fib_LM_cellphone/Tumor_fib","network.pdf")))
netVisual_chord_cell_internal(mat,show.legend =T,
                              sources.use =c("Fib(C1)-ABCC9", "Fib(C2)-MYH11", "Fib(C3)-THBS2", "Fib(C4)-CXCL8", "Fib(C5)-TOP2A"),
                              targets.use = c('Myeloid','NK',"MAIT","Plasma","Endothelial",  "Treg", 'B',"CD4","CD8","Epithelial"),
                              #label.edge= T,
                              #color.use =c( '#80B1D3','#FDB462','#B3DE69', '#FCCDE5','#E41A1C','#1B9E77','#377EB8','#984EA3','#E5C494','#FFD92F'), 
                              directional = 1)
dev.off()
#?netVisual_chord_cell_internal
##################
library(data.table)
library(scales)
library(tidyr)
library(tidyverse)
library(ktplots)
library(SingleCellExperiment)
#remotes::install_local("D:/nCov2019-master.zip",upgrade = F,dependencies = T)

means<-read.delim2("/home/data/gaoyuzhen/Projects/LM_spatial/out/means.txt", check.names = FALSE)
pvals <- read.delim2("/home/data/gaoyuzhen/Projects/LM_spatial/out/pvalues.txt", check.names = FALSE)
means2<-apply(means[,c(12:ncol(means))],2,as.numeric)
means2[means2<1]<-0
means<-cbind(means[,c(1:11)],means2)
pvals2<-apply(means[,c(12:ncol(pvals))],2,as.numeric)
pvals<-cbind(means[,c(1:11)],pvals2)
#plot_cpdb("B cell", "CD4T cell", kidneyimmune, 'celltype', means, pvals, split.by = "Experiment", gene.family = 'chemokines')
#viridis::scale_color_viridis(option="H")
#plot_cpdb("13", "14", kidneyimmune, 'cluster', means, pvals,  genes = c("CXCL13", "CD274", "CXCR5"))
##gene.family ===c('chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
##'



rownames(meta.data)<-meta.data$Cell
scRNA_sce<-AddMetaData(scRNA_sce,meta.data)
colnames(scRNA_sce@meta.data)
#?plot_cpdb
# “viridis”,  “magma”, “plasma”, “inferno”, “civids”, “mako”, and “rocket” 
plot_cpdb(cell_type1 ="Fib(C1)-ABCC9|Fib(C2)-MYH11|Fib(C3)-THBS2|Fib(C4)-CXCL8|Fib(C5)-TOP2A",
          cell_type2 ="Myeloid|Epithelial",
          #split.by='cell_type',
          scdata =scRNA_sce[,scRNA_sce$set=="LM"],
          idents ='cell_type', 
          means = means, 
          pvals = pvals, 
          #noir=T,
          col_option = viridis::mako(50),
          keep_significant_only = TRUE,
          #split.by = 'Experiment', 
          default_style = TRUE,
          #cluster_columns = T, 
          gene.family =c('chemokines','costimulatory', 'coinhibitory', 'Treg'),
          cluster_rows = F ) + small_guide() + small_axis(fontsize =6) + small_legend(fontsize = 8,keysize=.5) & viridis::scale_color_viridis(option="A")

ggsave(file.path(Figure4,"metabolism_DimPlot_LM_fib_cellphone.pdf"),width = 7,height =4)

