###network of metabolites and special genes in sp and scRNA
load(file.path(Figure8,"filtered_list_LM_caf_metabolites.Rdata"))
#############################################
st_sce<-readRDS("data/metabolits_scp_st_sp_sce.RDS")
############### fib genes  scRNA fib genes
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
sce.markers.x.sig<- sce.markers.x.sig %>% group_by(cluster) %>% top_n(n =100, wt = avg_log2FC)# top 10 
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
######################################
##sp genes
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
sce.markers.merge<-rbind(sce.markers1,sce.markers2)################################
table(sce.markers$cluster)
sce.markers.merge<-sce.markers.merge[sce.markers.merge$cluster=="LM_Fibroblast" & sce.markers.merge$p_val_adj<0.001,]
#########################################
a
b
c
c<-sce.markers.merge$gene

markergeness<-intersect(a,c)
markergeness
##################





sub_st_sce<-st_sce[,st_sce$st_celltype=="Fibroblast" ]
metabolismexp<-sub_st_sce@assays$metabolits@data[unique(final_results$metabolits),] %>% data.frame()
geneexp<-sub_st_sce@assays$SCT@data[markergeness,] %>% data.frame()

cor_meta_gene<-rbind(metabolismexp,geneexp) %>% data.frame()

cor(t(cor_meta_gene))

corrplot(corr=cor(t(cor_meta_gene)),
         #method = "color",
         #order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white","#FC4E07"))(50),
)

####
library(igraph)
library(graph)
g <- graph.adjacency(
  as.matrix(as.dist(cor(t(cor_meta_gene), method="pearson"))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

#Simplfy the adjacency object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- 'blue'

# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- '#E41A1C'

# Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)

# Change arrow size
# For directed graphs only
#E(g)$arrow.size <- 1.0

# Remove edges below absolute Pearson correlation 0.8
g <- delete_edges(g, E(g)[which(E(g)$weight<0.1)])

# Remove any vertices remaining that have no edges
g <- delete.vertices(g, degree(g)==0)

# Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name

# Change shape of graph vertices
V(g)$shape <- "sphere"

# Change colour of graph vertices
V(g)$color <- "skyblue"

# Change colour of vertex frames
V(g)$vertex.frame.color <- "white"

# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}
vSizes <- (scale01(apply(cor_meta_gene, 1, mean)) + 1.0) * 10

# Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 2.0

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="prim")

# Plot the tree object
plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph")


mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1

coords <- layout_(g, as_star())


par(mfrow=c(1,2))

plot(
  mst.clustering, mst,
  #layout=layout.fruchterman.reingold,
  layout= layout_with_graphopt,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.8,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=1,
  edge.width=edgeweights,
  edge.arrow.mode=1,
  main="top gene interaction with top metabolism")

plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="My first graph")




