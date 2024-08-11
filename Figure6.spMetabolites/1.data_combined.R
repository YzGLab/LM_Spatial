##
library(data.table)
setwd("/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/stRNA_scMetabolismuniondata")
path <- "/home/data/gaoyuzhen/Projects/LM_spatial/spatial_metabolism/stRNA_scMetabolismuniondata"
sample <- list.files(path)

for (samp in sample) {
  print(samp)
  samplename <- list()
  metadata <- list()
  transdata <- list()
  for (mode in c("neg","pos")) {
    samplename[[mode]] <- read.table(paste0(path,"/",samp,"/",mode,"/samplename.txt"),header = T,check.names=FALSE)
    metadata[[mode]] <- read.table(paste0(path,"/",samp,"/",mode,"/metadata.txt"),header = T,check.names=FALSE)
    transdata[[mode]] <- read.table(paste0(path,"/",samp,"/",mode,"/transdata.txt"),header = T,check.names=FALSE)
  }
  ## 正负模式空转barcodeid交集
  transname <- intersect(samplename[["neg"]]$transname,samplename[["pos"]]$transname)
  ## samplename文件处理
  dat <- data.frame(transname = transname,
                    neg = samplename$neg$metaname[match(transname,samplename$neg$transname)],
                    pos = samplename$pos$metaname[match(transname,samplename$pos$transname)])
  colneg <- c(colnames(metadata$neg)[1:13],dat$neg)
  neg <- metadata$neg %>% dplyr::select(all_of(colneg))
  colnames(neg)[14:ncol(neg)] <- dat$transname
  neg$mode <- "neg"
  colpos <- c(colnames(metadata$pos)[1:13],dat$pos)
  pos <- metadata$pos %>% dplyr::select(all_of(colpos))
  colnames(pos)[14:ncol(pos)] <- dat$transname
  pos$mode <- "pos"
  data <- rbind(neg,pos)
  ## 处理后的代谢数据,以barcodeid注释
  data <- data %>% dplyr::select("ID","mz","mode",everything())
  ## 处理后的转录数据
  coltrans <- c(colnames(transdata$neg)[1:2],dat$transname)
  transdatanew <- transdata$neg %>% dplyr::select(all_of(coltrans))
  ## 写出最后矩阵
  data0 <- list()
  data0[["metadata"]] <- data
  data0[["samplename"]] <- dat
  data0[["transdata"]] <- transdatanew
  # write.xlsx(data0,paste0(samp,"_both.xlsx"))
  saveRDS(data0,paste0(path,"/",samp,"_both.rds"))
}
