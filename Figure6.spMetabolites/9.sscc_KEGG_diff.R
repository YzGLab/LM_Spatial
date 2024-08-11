library(FELLA)
data("FELLA.sample")


graph <- buildGraphFromKEGGREST(
  organism = "hsa", 
  filter.path = c("01100", "01200", "01210", "01212", "01230"))

tmpdir <- paste0(tempdir(), "/my_database")
# Mke sure the database does not exist from a former vignette build
# Otherwise the vignette will rise an error 
# because FELLA will not overwrite an existing database
unlink(tmpdir, recursive = TRUE)  
buildDataFromGraph(
  keggdata.graph = graph, 
  databaseDir = tmpdir, 
  internalDir = FALSE, 
  matrices = "diffusion", 
  normality = "diffusion", 
  niter = 50)

###################################################
### code chunk number 4: 01_loadkegg
###################################################
fella.data <- loadKEGGdata(
  databaseDir = tmpdir, 
  internalDir = FALSE, 
  loadMatrix = "diffusion"
)
save(fella.data,file="data/fella_data_metabolits_KEGG.Rdata")


###
sample_names <- c("LM_1", "LM_2", "C_1", "C_2")

#sample<-

cluster_names<-c("cluster3","cluster6","cluster8")


sample<-sample_names[1]
clusters<-cluster_names[1]
resutls_cluster_metabolits<-paste0("matched_",sample,clusters,".csv")


cluster_diff<-read.csv(paste0("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure7_sp/",resutls_cluster_metabolits))
cluster_diff.sig <- cluster_diff %>% group_by(class) %>% top_n(n =50, wt = statistic)# top 10

classnames<-unique(cluster_diff.sig$class)
class_clusters<-classnames[3]
kegg<-cluster_diff.sig[cluster_diff.sig$class==class_clusters,]$KEGG
# 去除单独的分号和双分号
filtered_kegg <- kegg[!kegg %in% c(";", ";;")]
# 将含有分号的字符串分裂开来
split_kegg <- unlist(strsplit(filtered_kegg, ";"))

myAnalysis <- enrich(
  compounds = split_kegg, 
  method = "diffusion", 
  approx = "normality", 
  data = fella.data)

plot_names<-paste0(sample," ", clusters,"_",class_clusters,"_FELLA")
plot(
  x = myAnalysis, 
  method = "diffusion", 
  main = plot_names, 
  threshold = 0.01, 
  nlimit = 60,
  data = fella.data)

pdf(file.path(Figure7,paste0( plot_names, ".pdf")),width = 10,height = 10)

dev.off()

myTable <- generateResultsTable(
  object = myAnalysis, 
  method = "diffusion", 
  threshold = 0.01, 
  data = fella.data)

write.csv(myTable,file=paste0(sample,"_", clusters,"_",class_clusters,"_FELLA_KEGG.csv"))




##############
# 加载必要的库
library(dplyr)

# 定义函数
process_sample_cluster <- function(sample, cluster) {
  # 生成文件名
  resutls_cluster_metabolits <- paste0("matched_", sample, cluster, ".csv")
  print(resutls_cluster_metabolits)
  
  # 读取数据
  cluster_diff <- read.csv(paste0("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure7_sp/", resutls_cluster_metabolits))
  cluster_diff.sig <- cluster_diff %>% group_by(class) %>% top_n(n = 30, wt = statistic)
  cluster_diff.sig<-cluster_diff.sig[cluster_diff.sig$statistic>10,]
  # 处理每个class中的KEGG
  classnames <- unique(cluster_diff.sig$class)
  for (class_clusters in classnames) {
    kegg <- cluster_diff.sig[cluster_diff.sig$class == class_clusters,]$KEGG
    
    # 去除单独的分号和双分号
    filtered_kegg <- kegg[!kegg %in% c(";", ";;")]
    # 将含有分号的字符串分裂开来
    split_kegg <- unlist(strsplit(filtered_kegg, ";"))
    
    # 进行富集分析，添加错误处理
    tryCatch({
      myAnalysis <- enrich(
        compounds = split_kegg, 
        method = "diffusion", 
        approx = "normality", 
        data = fella.data
      )
      
      plot_names <- paste0("FELLA_", sample, "_", cluster, "_", class_clusters)
      # 绘图
      pdf(file.path(Figure7, paste0(plot_names, ".pdf")), width = 10, height = 10)
      plot(
        x = myAnalysis, 
        method = "diffusion", 
        main = "Diffusion analysis in FELLA", 
        threshold =0.001, 
        limit =10,
        data = fella.data,
        vertex.label.cex = 0.5
      )
      dev.off()
      
      # 生成结果表
      myTable <- generateResultsTable(
        object = myAnalysis, 
        method = "diffusion", 
        threshold =0.001, 
        data = fella.data
      )
      
      # 写入CSV文件
      write.csv(myTable, file = file.path(Figure7, paste0("FELLA_KEGG_", sample, "_", cluster, "_", class_clusters, ".csv")))
    }, error = function(e) {
      # 打印错误信息，并继续下一个组合
      message(paste("Error in processing:", sample, cluster, class_clusters, ":", e$message))
    })
  }
}

# 样本和群集名称
sample_names <- c("LM_1", "LM_2", "C_1", "C_2")
cluster_names <- c("cluster3", "cluster6", "cluster8")

# 循环处理
for (sample in sample_names) {
  for (cluster in cluster_names) {
    process_sample_cluster(sample, cluster)
  }
}

# 样本和群集名称
sample_names <- c("LM_1", "LM_2", "C_1", "C_2")
cluster_names <- c("cluster3", "cluster6", "cluster8")

# 计算Entry.type中pathway的含量
pathway_counts <- list()
myTable <- read.csv("ResultsLM/Figure7_sp/FELLA_0.01/FELLA_KEGG_C_2_cluster8_1.csv")
for (sample in sample_names) {
  for (cluster in cluster_names) {
    csv_filename <- paste0("FELLA_KEGG_", sample, "_", cluster, "_", 8, ".csv")
    print(csv_filename)
    
    tryCatch({
      myTable <- read.csv(paste0("ResultsLM/Figure7_sp/FELLA_0.01/", csv_filename))
      pathway_count <- sum(myTable$Entry.type == "pathway")
      pathway_counts[[csv_filename]] <- pathway_count
    }, error = function(e) {
      message(paste("Error in reading or processing:", csv_filename, ":", e$message))
    })
  }
}



# 对数据框按p.score进行排序
myTable <- myTable %>% arrange(p.score)
ggplot(myTable, aes(x = -log10(p.score), y = reorder(KEGG.name, -p.score))) +
  geom_bar(stat = "identity",position=position_dodge(width=0.6),width=0.8) +
  labs(x = "Log10(p.score)", y = "KEGG.name") +
  #theme_minimal() +
  scale_color_manual(values="red") +
  theme(axis.text.y = element_text(size = 10, hjust = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 30, size = 8),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")
#stat_compare_means(label = "p.signif")
if(length(twotype)>0){
  ggsave(file.path(fig.path6,paste0("pro_diff/",genes,"_protein_Inpancancer.pdf")),width =3*length(twotype),height = 6 )
} else {
  ggsave(file.path(fig.path6,paste0("pro_diff/",genes,"_protein_Inpancancer.pdf")),width =2,height = 6 )
}

#######################导出数据
myExp_csv <- "cluster_table.csv"
exportResults(
  format = "csv", 
  file = myExp_csv, 
  method = "diffusion", 
  threshold = 0.1, 
  object = myAnalysis, 
  data = fella.data)


myExp_graph <- paste0(myTempDir, "/graph.RData")
exportResults(
  format = "igraph", 
  file = myExp_graph, 
  method = "pagerank", 
  threshold = 0.1, 
  object = myAnalysis, 
  data = fella.data)