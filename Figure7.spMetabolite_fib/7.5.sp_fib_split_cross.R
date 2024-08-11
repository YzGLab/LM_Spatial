####

load(file.path(Figure8,"significant_features.Rdata"))


# 返回结果列表
names(feat.rank.mat.list)
# 定义筛选条件的函数
filter_criteria <- function(data) {
  data[data[, "orthoVipVn"] > 1 & 
         log2(data[, "foldchange"]) > 0.2 & 
         data[, "pvalue"] < 0.05, ]
}




# 筛选各个列表中的代谢物
filtered_list <- lapply(feat.rank.mat.list, filter_criteria)



# 找出每个样本的三个比较结果中的共同代谢物
get_common_metabolites <- function(list) {
  metabolites_names_list <- lapply(list, rownames)
  common_metabolites <- Reduce(intersect, metabolites_names_list)
  
  # 合并共同代谢物的统计值
  combined_data <- do.call(rbind, lapply(names(list), function(name) {
    data <- list[[name]]
    #data$metabolits<-rownames(data)
    common_data <- data[rownames(data) %in% common_metabolites, ]
    common_data$Comparison <- name
    common_data
  }))
  
  combined_data
}

filtered_list<-filtered_list[c(1,2,3)]
final_results<-get_common_metabolites(filtered_list)
final_results

names(filtered_list)
save(filtered_list,file=file.path(Figure8,"filtered_list_LM_caf_metabolites.Rdata"))


load(file.path(Figure8,"filtered_list_LM_caf_metabolites.Rdata"))
# Export each list as a CSV file
for (name in names(filtered_list)) {
  tem<-filtered_list[[name]]
  names(tem)
  
  write.csv(filtered_list[[name]], file=file.path(Figure8,paste0(name, "diff_metabolites.csv")), row.names = FALSE)
}

# Assuming 'x' is your first data frame and 'filtered_list' contains your second data frames

# Function to merge and export each data frame in filtered_list with 'x'
export_merged_csv <- function(x, filtered_list) {
  for (name in names(filtered_list)) {
    tem <- filtered_list[[name]]
    
    # Merge based on 'mz' in 'x' and 'metabolits' in 'tem'
    merged_df <- merge(tem, x, by.x = "metabolits", by.y = "mz")
    
    # Export merged data frame to CSV
    write.csv(merged_df, file.path(Figure8,paste0(name, "diff_metabolites.csv")), row.names = FALSE)
  }
}

# Example usage
export_merged_csv(x, filtered_list)



library(ggvenn)
library(ggVennDiagram)
# 提取每个列表中的代谢物名称
metabolites_list <- lapply(filtered_list, rownames)
names(metabolites_list) <- names(filtered_list)

VennDiagram::venn.diagram(filtered_list)
# 创建维恩图
library(VennDiagram)
# 创建维恩图
venn.plot <- venn.diagram(
  x = metabolites_list,
  category.names = names(metabolites_list),
  filename = NULL,
  output = TRUE,
  imagetype = "png",
  height = 480, 
  width = 480, 
  resolution = 300,
  lwd = 2,
  fill = c("#999999", "#E69F00", "#56B4E9"),
  alpha = 0.5,
  label.col = "black",
  cex = 1.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "outer",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  cat.dist = c(0.05, 0.05, 0.05),
  cat.pos = 0
)

# 保存图像
pdf("VennDiagram.pdf")
grid.draw(venn.plot)
dev.off()

# 输出最终结果
filtered_list1<-filtered_list[c(1,3)]
final_results1<-get_common_metabolites(filtered_list1)
final_results1

filtered_list2<-filtered_list[c(2,3)]
final_results2<-get_common_metabolites(filtered_list2)
final_results2
final_results <- rbind(final_results1, final_results2)

# 输出最终结果
final_results
###
final_results<-final_results[final_results$metabolits %in% rownames(feat.rank.mat.list[["LM_vs_CRC"]]),]
save(filtered_list,final_results,file=file.path(Figure8,"filtered_list_LM_caf_metabolites.Rdata"))
##################################################################
a<-unique(final_results$metabolits)
#####################################################################
dev.off()
sub_sce<-sce[,sce$st_celltype %in% c("Fibroblast", "Hepatocytes","Tumor","Normal Epi")]
sub_sce$groups<-ifelse(sub_sce$st_celltype=="Fibroblast",1,0)
VlnPlot(sub_sce[,sub_sce$set=="LM"],features =unique(final_results$metabolits)[1:20],pt.size = 0 ,group.by = "st_celltype",cols =stcol)
ggsave("VlnPlot_target_metabolites.pdf")

VlnPlot(sub_sce[,sub_sce$set=="LM"],features =unique(final_results$metabolits)[11:20],pt.size = 0 ,group.by = "groups")
stcol<-c('#2ca02c', '#d62728', '#1f77b4', '#ff7f0e', '#9467bd', '#f5e801', '#08f7f0')


VlnPlot(sub_sce[,sub_sce$orig.ident=="LM_2"],features =unique(rownames(final_results))[1:4],pt.size = 0 ,group.by = "st_celltype")
