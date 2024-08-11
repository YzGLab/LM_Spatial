#
##
tissue_sp<-read.delim2("/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/tissue_metabolism/GVSCC_diff_all_features.txt")
###################
head(tissue_sp)
############################
tissue_sp$HMDB

#"HMDB0011154"                             NA                                        "HMDB10404"                               "HMDB10404"                              
# "HMDB10404"                               "HMDB11482"                               "HMDB11482"                               "HMDB11487;HMDB11488;HMDB11517;HMDB11518" "HMDB11487;HMDB11488;HMDB11517;HMDB11518" "HMDB11494;HMDB11495;HMDB11524;HMDB11525"

matched_data.x$Compound.ID
# 假设match.data.x和tissue_sp已经被加载为数据框

colnames(matched_data.x)
colnames(tissue_sp)
# 先将match.data.x的Compound.ID按分号分割成多个代谢物ID
library(dplyr)
library(tidyr)

# 分割Compound.ID并将每个代谢物单独列出来
match.data.x <- matched_data.x %>%
  separate_rows(Metabolites, sep = ";")

# 对tissue_sp的HMDB列也进行同样的处理
tissued_sp <- tissue_sp %>%
  separate_rows(MS2.name, sep = ";")

# 对tissue_sp的MS1.name列进行处理，按逗号分割
tissued_sp <- tissued_sp %>%
  separate_rows(MS1.name, sep = ",")

# 用inner_join进行匹配，保留一致的代谢物（HMDB匹配）
matched_hmdb <- match.data.x %>%
  inner_join(tissued_sp, by = c("Metabolites" =  "MS2.name"))

# 用inner_join进行匹配，保留一致的代谢物（MS1.name匹配）
matched_ms1 <- match.data.x %>%
  inner_join(tissued_sp, by = c("Metabolites" = "MS1.name"))

# 合并匹配结果
matched_data <- bind_rows(matched_hmdb, matched_ms1) %>%
  distinct()


write.csv(matched_data,file=file.path(Figure8,"matched_data_sp_tissue.csv"))

# 查看匹配结果
print(matched_data)


###################ROC analysis



#ROC

# Install and load necessary libraries
#install.packages("pROC")
#install.packages("dplyr")
#install.packages("ggplot2")
getwd()

roc_results<-read.delim2("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure8_feature_fib_sp/tissue_results.txt")
# 假设 roc_results 是你的数据框
# 使用 dplyr 包进行操作
library(pROC)
library(dplyr)
library(ggplot2)

# 假设 roc_results 是你的数据框
# 创建新的 GroupType 列
roc_results <- roc_results %>%
  mutate(GroupType = ifelse(grepl("^G", Group), 1, 0))

# 检查新的 GroupType 列是否创建正确
head(roc_results)

roc_results$GroupType<-as.numeric(roc_results$GroupType)

# 要分析的代谢物列表
metabolites <- c("Suberic.acid", "Aspartyl.Aspartate", "Tetraethylene.glycol", "Triacetin", 
                 "N1.Acetylspermine", "Fulvestrant", "Muricatenol")

# 为每个代谢物生成 ROC 曲线的函数
plot_roc <- function(data, metabolite) {
  data[[metabolite]]<-as.numeric(data[[metabolite]]                            )
  roc_obj <- roc(data$GroupType, data[[metabolite]],plot=TRUE, print.thres=FALSE, print.auc=TRUE,ci=TRUE,
                 xlim=c(1,0),ylim=c(0,1), smooth=TRUE,col="#6BAED6")
  #plot(roc_obj, main = paste("ROC for", metabolite))
  return(auc(roc_obj))  # 返回 AUC 作为参考
}



# 为 LM 组中的每个代谢物生成 ROC 曲线
cat("LM Group ROC Curves:\n")
for (metabolite in metabolites) {
  cat("Plotting ROC for LM group, metabolite:", metabolite, "\n")
  data[[metabolite]]<-as.numeric(data[[metabolite]]                            )
  pdf(file.path(Figure8,paste0("pro_dif_",metabolite,"_ROC.pdf")),width =4,height = 4)
  roc_obj <- roc(data$GroupType, data[[metabolite]],plot=TRUE, print.thres=FALSE, print.auc=TRUE,ci=TRUE,
                 xlim=c(1,0),ylim=c(0,1), smooth=TRUE,col="#6BAED6")
  dev.off()
  #print(plot_roc(roc_results, metabolite))
}



# 使用 GLM 合并代谢物，获取总的代谢物评分 ROC
glm_model <- glm(GroupType ~ ., data = roc_results[, c("GroupType", metabolites[-1])], family = binomial)

# 预测概率
roc_results$predicted_score <- predict(glm_model, type = "response")

# 生成总的代谢物评分的 ROC 曲线
roc_combined <- roc(roc_results$GroupType, roc_results$predicted_score,plot=TRUE, print.thres=FALSE, print.auc=TRUE,ci=TRUE,
                    xlim=c(1,0),ylim=c(0,1),col="#6BAED6")

combined_auc <- auc(roc_combined)

cat("Combined AUC:", combined_auc, "\n")
