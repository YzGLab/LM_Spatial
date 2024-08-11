#####sp-fib_metabolits
library(Cardinal) 
library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)
#save(sce.cell_type_markers,sce.markers.epi,sce.markers.epi.sig,file='data/spatial_scRNA_differetgenes_clusters.Rdata')
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
meta_data_st<-st_sce@meta.data
rownames(meta_data_st)<-gsub("-1","",rownames(meta_data_st))
#rownames(meta_data_st)<-gsub("-1","",rownames(meta_data_st))
rownames(meta_data_st) <- gsub("(.+?)_(.+?)_", "\\1_\\2.", rownames(meta_data_st))
# 定义样本名称
sample_names <- c("LM_1", "LM_2", "C_1", "C_2")

# 初始化一个空列表来存储MSImagingExperiment对象
msi_data_list <- list()

for (sample in sample_names) {
  # 读取数据
  print(sample)
  samples_both <- readRDS(paste0("/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/stRNA_scMetabolismuniondata/", sample, "_both.rds"))
  
  # 根据样本名称选择列
  if (sample == "LM_1") {
    data_matrix <- samples_both$metadata[, c(2, 15:3158)]
  } else if (sample == "LM_2") {
    data_matrix <- samples_both$metadata[, c(2, 15:3399)]
  } else if (sample == "C_1") {
    data_matrix <- samples_both$metadata[, c(2, 15:213)]
  } else if (sample == "C_2") {
    data_matrix <- samples_both$metadata[, c(2, 15:3957)]
  }
  
  rownames(data_matrix) <- data_matrix[, 1]
  data_matrix <- data_matrix[, -1]
  data_matrix[] <- apply(data_matrix, 2, as.numeric)
  data_matrix <- data_matrix[, colnames(data_matrix) %in% samples_both$samplename$transname]
  
  # 提取坐标信息
  positions <- samples_both$samplename[samples_both$samplename$transname %in% colnames(data_matrix), ]
  coord <- data.frame(
    x = as.numeric(sub(".*-(\\d+)-.*", "\\1", positions$pos)),
    y = as.numeric(sub(".*-(\\d+)$", "\\1", positions$pos))
  )
  rownames(coord) <- positions$transname
  
  # 确保 mz 按升序排列并与 data_matrix 的行数匹配
  mz <- sort(as.numeric(rownames(data_matrix)))
  
  if (length(mz) != nrow(data_matrix)) {
    stop("The length of mz does not match the number of rows in data_matrix.")
  }
  
  # 确保 coord 的行数与 data_matrix 的列数匹配
  if (nrow(coord) != ncol(data_matrix)) {
    stop("The number of rows in coord does not match the number of columns in data_matrix.")
  }
  
  # 检查 data_matrix 的类型和维度
  if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }
  # 创建 ImageArrayList 对象
  idata <- ImageArrayList(data_matrix)
  # 创建 MassDataFrame 对象
  fdata <- MassDataFrame(mz = mz)
  
  # 创建 PositionDataFrame 对象
  run <- factor(rep(sample, nrow(coord)))
  pdata <- PositionDataFrame(run = run, coord = coord)
  st_data<- meta_data_st[meta_data_st$orig.ident==sample ,]########添加空间转录组学信息
  st_data<-st_data[rownames(st_data) %in% colnames(data_matrix),]
  #metadata(pdata)<-st_data
  # 创建 MSImagingExperiment 对象
  msi_data <- MSImagingExperiment(
    imageData=idata,
    featureData=fdata,
    pixelData=pdata,
    metadata=st_data
  )
  
  # 将对象添加到列表中
  msi_data_list[[sample]] <- msi_data
  
  # 绘制图像
  image(msi_data)
}


# 合并所有样本的 MSImagingExperiment 对象
#combinedImgSet <- combine(msi_data_list[[1]],msi_data_list[[2]],)
combinedImgSet <- combine(msi_data_list[[3]],msi_data_list[[4]],msi_data_list[[1]],msi_data_list[[2]])
#combinedImgSet <- combine(msi_data_list[[4]])
# 保存对象
save(msi_data_list, combinedImgSet, file = "data/sp_cardnial_dataset_full_four_samples.Rdata")

image(combinedImgSet, 
      #col = stcol,
      #normalize.image =  "linear",
      #contrast.enhance =  "suppression",
      colorscale = viridis,
      #smooth.image =  "gaussian",
      asp=1.5,
      alpha.power = 1.5,
      xlim=c(-10,125)
)


image(msi_data_list[[4]], asp=1.5,feature.groups="st_celltype",superpose = F)

###############################################
set.seed(2020)
###
msi_data<-msi_data_list[[1]]
# 基于 celltype 过滤数据
celltypes <- unique(msi_data@elementMetadata@metadata$st_celltype)

##############################
set.seed(1)
mse <- msi_data
cls <- msi_data@metadata$st_celltype
# fit a PLS model with 3 components
pls <- OPLS(mse,cls)
plot(pls, type="scores")
# visualize predictions
## S4 method for signature 'SpatialPLS'
topFeatures(pls, n = Inf, sort.by ="vip")
## S4 method for signature 'SpatialPLS,missing'

samples <-pixelData(msi_data)$run
msi_data@elementMetadata@metadata$run<-ifelse(msi_data@elementMetadata@metadata$st_celltype=="Fibroblast",1,0)
pixelData(msi_data)$run<-msi_data@metadata$st_celltype
pixelData(msi_data)$run<-as.numeric(factor(pixelData(msi_data)$run))
fit <- meansTest(
  x = msi_data, 
  #data=st_data,
  ~ run, 
  #random = NULL,  # 如果没有随机效应，可以设为 NULL
  #samples = samples, 
  #response = "intensity"
)
topFeatures(fit, n = Inf, sort.by = "statistic")
msi_data@metadata$design$pixelData@run
table(msi_data@elementMetadata@metadata$run)
# 查看结果
str(msi_data)
str(x)
# 假设 msi_data 是你的 MSImagingExperiment 对象，并且 celltype 是你的样本元数据
# 筛选出感兴趣的细胞类型
selected_cells <- msi_data@elementMetadata@metadata$st_celltype %in% c("Tumor", "Fibroblast")

###
samples <- replace(msi_data@elementMetadata@metadata$st_celltype, !selected_cells, NA)
samples<-as.factor(samples)


# 检查 selected_cells 的长度是否与 msi_data 的行数相同
print(length(selected_cells))
print(ncol(msi_data))
# 更新 MSImagingExperiment 对象，仅保留所选样本
filtered_msi_data <- msi_data[,selected_cells ]
filtered_msi_data$st_celltype
filtered_msi_data$st_celltype<-factor(filtered_msi_data$st_celltype)
metadata(filtered_msi_data@metadata)
# 执行差异分析
fit <- meansTest(
  x = msi_data, 
  ~ st_celltype, 
  #random = NULL,  # 如果没有随机效应，可以设为 NULL
  samples = samples, 
  response = "intensity"
)


summary(fit)

