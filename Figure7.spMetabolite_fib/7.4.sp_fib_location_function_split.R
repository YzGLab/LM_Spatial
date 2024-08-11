sce<-readRDS("data/metabolits_scp_st_sp_sce.RDS")
metadata<-sce@meta.data
metabolites_obj<-sce@assays$metabolits@data %>% data.frame(.,check.names = FALSE)
write.csv(metabolites_obj, file=file.path(Figure8, "full_metabolites_obj_four_samples_metabolites.csv"), row.names = TRUE)
# Load necessary library
library(MetaboAnalystR)

# 定义样本名称和比较组
dataset_names <- c("LM_1", "LM_2")
comparison_groups <- list(
  LM_1 = list(c("Fibroblast", "Hepatocytes"), c("Fibroblast", "Tumor"), c("Fibroblast", "Normal Epi")),
  LM_2 = list(c("Fibroblast", "Hepatocytes"), c("Fibroblast", "Tumor"), c("Fibroblast", "Normal Epi"))
)

# 初始化结果存储列表
feat.rank.mat.list <- list()

# 遍历每个样本
for (sample in dataset_names) {
  
  print(sample)
  
  # 遍历当前样本的每个比较组
  for (comparison in comparison_groups[[sample]]) {
    
    celltype1 <- comparison[1]
    celltype2 <- comparison[2]
    
    # 过滤当前样本和比较组的元数据
    metadata_sample <- metadata[metadata$orig.ident == sample & metadata$st_celltype %in% c(celltype1, celltype2), ]
    
    # 创建标签
    Lable <- data.frame(Sample = rownames(metadata_sample), 
                        Lable = ifelse(metadata_sample$st_celltype == "Fibroblast", 1, 0))
    
    # 子集化代谢物数据并转置
    Metabolites_sample <- t(metabolites_obj[, colnames(metabolites_obj) %in% rownames(metadata_sample)])
    Metabolites_sample <- data.frame(Sample = rownames(Metabolites_sample), Metabolites_sample)
    
    # 合并标签和代谢物数据
    Metadat <- merge(Lable, Metabolites_sample, by = "Sample", all = FALSE)
    
    # 创建目录
    sample_dir <- paste0("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure8_feature_fib_sp/results_diff_meta_", sample, "_", "Fibroblast", "_vs_", celltype2)
    if (!dir.exists(sample_dir)) {
      dir.create(sample_dir)
    }
    setwd(sample_dir)
    
    # 写入合并数据到CSV文件
    write.csv(Metadat, file = paste0(sample_dir, "/metadata.csv"), row.names = FALSE)
    rm(mSet)
    
    # 初始化 MetaboAnalystR
    mSet <- InitDataObjects("conc", "stat", FALSE)
    mSet <- Read.TextData(mSet, paste0(sample_dir, "/metadata.csv"), "rowu", "disc")
    mSet <- SanityCheckData(mSet)
    mSet <- ReplaceMin(mSet)
    mSet <- RemoveMissingPercent(mSet, percent = 0.2)
    mSet <- ImputeMissingVar(mSet, method = "min")
    mSet <- SanityCheckData(mSet)
    mSet <- FilterVariable(mSet, "F", 25, "none", -1, "mean", 0)
    mSet <- PreparePrenormData(mSet)
    mSet <- Normalization(mSet, "NULL", "NULL", "NULL", ratio = FALSE, ratioNum = 20)
    mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width = NA)
    mSet <- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width = NA)
    
    # 统计分析
    mSet <- FC.Anal(mSet, 2.0, 1, FALSE)
    mSet <- PlotFC(mSet, "fc_1_", "pdf", 72, width = NA)
    mSet <- Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
    mSet <- PlotTT(mSet, "tt_1_", "PDF", 72, width = NA)
    mSet <- Volcano.Anal(mSet, FALSE, 2.0, 1, F, 0.1, TRUE, "raw")
    mSet <- PlotVolcano(mSet, "volcano_1_", 1, 0, "PDF", 72, width = NA, -1)
    mSet <- OPLSR.Anal(mSet, reg = TRUE)
    mSet <- PlotOPLS2DScore(mSet, "opls_score2d_0_", "PDF", 72, width = NA, 1, 2, 0.95, 0, 0, "na")
    mSet <- PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "PDF", 72, width = NA)
    mSet <- PlotOPLS.Imp(mSet, "opls_imp_0_", "PDF", 72, width = NA, "vip", "tscore", 20, FALSE)
    mSet <- PlotOPLS.MDL(mSet, "opls_mdl_0_", "PDF", 72, width = NA)
    mSet <- PlotOPLS2DScore(mSet, "opls_score2d_1_", "PDF", 72, width = NA, 1, 2, 0.95, 0, 1, "na")
    mSet <- UpdateOPLS.Splot(mSet, "all")
    mSet <- PlotOPLS.Splot(mSet, "opls_splot_1_", "all", "PDF", 72, width = NA)
    mSet <- PlotOPLS.Imp(mSet, "opls_imp_1_", "PDF", 72, width = NA, "vip", "tscore", 20, FALSE)
    mSet <- PCA.Anal(mSet)
    mSet <- PlotPCAPairSummary(mSet, "pca_pair_1_", "PDF", 72, width = NA, 5)
    mSet <- PlotPCAScree(mSet, "pca_scree_1_", "PDF", 72, width = NA, 5)
    mSet <- PlotPCA2DScore(mSet, "pca_score2d_1_", "PDF", 72, width = NA, 1, 2, 0.95, 0, 0)
    
    # 提取显著特征
    feat.rank.mat <- cbind(
      orthoVipVn = mSet$analSet$oplsda$orthoVipVn,
      foldchange = mSet$analSet$fc$fc.all,
      pvalue = mSet$analSet$tt$p.value
    )
   
    feat.rank.mat <- data.frame(feat.rank.mat)
    feat.rank.mat$metabolits<-gsub("X","",rownames(feat.rank.mat))
    feat.rank <- feat.rank.mat[
      feat.rank.mat[, 1] > 1 & 
        abs(log2(feat.rank.mat[, 2])) > 0.5 & 
        feat.rank.mat[, 3] < 0.001, 
    ]
    
    # 存储结果到列表
    feat.rank.mat.list[[paste0(sample, "_", celltype1, "_vs_", celltype2)]] <- feat.rank.mat
    
    # 保存显著特征到CSV文件
    write.csv(feat.rank.mat, file = paste0(sample_dir, "/all_features.csv"), row.names = TRUE)
    write.csv(feat.rank, file = paste0(sample_dir, "/significant_features.csv"), row.names = TRUE)
    rm(feat.rank.mat)
    rm(feat.rank)
  }
}

# 返回结果列表
names(feat.rank.mat.list)


# 定义筛选条件的函数
filter_criteria <- function(data) {
  data[data[, "orthoVipVn"] > 1 & 
         log2(data[, "foldchange"]) > 0 & 
         data[, "pvalue"] < 0.05, ]
}

# 筛选各个列表中的代谢物
filtered_list <- lapply(feat.rank.mat.list, filter_criteria)

# 根据样本名称提取子列表
LM_1_list <- filtered_list[grepl("LM_1", names(filtered_list))][c(1,3)]
LM_2_list <- filtered_list[grepl("LM_2", names(filtered_list))][c(1,3)]

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

# 获取 LM_1 和 LM_2 的共同代谢物
LM_1_common <- get_common_metabolites(LM_1_list)
LM_2_common <- get_common_metabolites(LM_2_list)

# 合并 LM_1 和 LM_2 的结果
final_results <- rbind(LM_1_common, LM_2_common)

# 输出最终结果
final_results



