###分群热图
library(ggSCvis)
library(scales)
library(tidyr)
library(tidyverse)
library(patchwork)
DefaultAssay(sub_sce)<- "metabolits"

sub_sce<-sce[,sce$st_celltype=="Fibroblast"]
##
#######

Idents(sub_sce)<-sub_sce$set
sce.markers.x <- FindAllMarkers(sub_sce, 
                                only.pos = TRUE,  test.use = "wilcox",
                                min.pct = 0.0, logfc.threshold = 0.0)
sce.markers.sig <- sce.markers.x %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10

table(sce.markers.x$cluster)

##########
library(ggplot2)
######special volcano plot
x<-sce.markers.x
x$label<- rownames(x)
head(x)
x$logFC<-log2(x$avg_log2FC)
x$P.Value<-x$p_val_adj


#plot_mode <- "classic" #
plot_mode <- "advanced" #
#
logFCcut <- 1.5 #log2-foldchange
pvalCut <- 0.05 #P.value
adjPcut <- 0.05 #adj.P.value

#for advanced mode
logFCcut2 <- 2.5
logFCcut3 <- 5
pvalCut2 <- 0.0001
pvalCut3 <- 0.00001

#
xmin <- (range(x$logFC)[1]- (range(x$logFC)[1]+ 10))
xmax <- (range(x$logFC)[1]+ (10-range(x$logFC)[1]))
ymin <- 0
ymax <- -log10(x$P.Value)[3] * 1.1


mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13",
           "#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D",
           "#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
if (plot_mode == "classic"){
  # setting for color
  x$color_transparent <- ifelse((x$P.Value < pvalCut & x$logFC > logFCcut), "red", ifelse((x$P.Value < pvalCut & x$logFC < -logFCcut), "blue","grey30"))
  # ??????setting for size
  size <- ifelse((x$P.Value < pvalCut & abs(x$logFC) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  # setting for color
  n1 <- length(x[, 1])
  cols <- rep("grey30", n1)
  names(cols)<- rownames(x)
  
  cols[x$P.Value < pvalCut & x$logFC >logFCcut]<- "#FB9A99"
  cols[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- "#ED4F4F"
  cols[x$P.Value < pvalCut & x$logFC < -logFCcut]<- "#B2DF8A"
  cols[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- "#329E3F"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  # setting for size
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  
  #?????????????
  size[x$P.Value < pvalCut & x$logFC > logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC > logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC > logFCcut3]<- 6
  size[x$P.Value < pvalCut & x$logFC < -logFCcut]<- 2
  size[x$P.Value < pvalCut2 & x$logFC < -logFCcut2]<- 4
  size[x$P.Value < pvalCut3 & x$logFC < -logFCcut3]<- 6
  
} else {
  stop("Unsupport mode")
}

# Construct the plot object
p1 <- ggplot(data=x, aes(logFC, -log10(P.Value), label = label)) +
  geom_point(alpha = 0.8, size = size, colour = x$color_transparent) +
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  #ylim(c(ymin,ymax)) + 
  #scale_x_continuous(
  # breaks = c(-10, -5, -logFCcut, 0, logFCcut, 5, 10), 
  # labels = c(-10, -5, -logFCcut, 0, logFCcut, 5, 10),
  #limits = c(-11, 11) 
  #) +
  #xlim(c(xmin, xmax)) + 
  #geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
  #    linetype="longdash", lwd = 0.5) + 
  # geom_hline(yintercept = -log10(pvalCut), color="grey40", 
  #linetype="longdash", lwd = 0.5) +

theme_bw(base_size = 12#, base_family = "Times" #???????
) +
  theme(panel.grid=element_blank())
p1

if (plot_mode == "advanced") {
  p1 <- p1 + 
    geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey40", 
               linetype="longdash", lwd = 0.5) +
    geom_hline(yintercept = -log10(pvalCut2), color="grey40", 
               linetype="longdash", lwd = 0.5)
}
p1
ggsave("GSE28425_volcano1.pdf",width = 6,height=7)
##########interation of genes in the BCC tumor and adjucent tissue
library(ggplot2)
library(ggrepel)
library(ggthemes)
n = 5.5
p1 + geom_text_repel(data=x,aes(x = logFC, y = -log10(P.Value), 
                                label = ifelse(abs(logFC) > n, rownames(x),"")),
                     colour="darkred", size = 5, box.padding = unit(0.35, "lines"), 
                     point.padding = unit(0.3, "lines"))

####
ggsave("GSE28425_volcano.pdf",width = 6,height=7)
selectgenes<-data.frame(rownames(x))

###
targetgenes<<-x[x$logFC>1.5,]

targetgenes<-c("GJB6","NTRK3","CBX2","STMN1","MDK", "TOP2A","CENPF","KIF20A")


p2 <- p1 + 
  geom_point(data = selectgenes, alpha = 1, size = 4.6, shape = 1, 
             stroke = 1, 
             color = "black") +
  scale_color_manual(values = mycol) + 
  geom_text_repel(data = selectgenes, 
                  show.legend = FALSE, 
                  size = 5, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines")) +
  guides(color=guide_legend(title = NULL)) 

p2


