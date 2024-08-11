library(SPATA2)
library(SPATAData)
library(tidyverse)
###
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
###
Colon_data_gan1 <-
  initiateSpataObject_10X(
    directory_10X = '/home/data/gaoyuzhen/Projects/LM_spatial/spatial/LM_1', # the directory from which to load the data
    sample_name = "gan1"
  )
getDirectoryInstructions(object = Colon_data_gan1)
# set/change the current default directory
Colon_data_gan1 <- setSpataDir(Colon_data_gan1, dir = "data/spata_objects/Colon_spatia_gan1.RDS")
# quickly save the `spata2` object under the default directory 
saveSpataObject(object = Colon_data_gan1)
##############################
Colon_data_gan1 <- setActiveMatrix(object = Colon_data_gan1, mtr_name = "scaled")
##机器学习去除噪声 
Colon_data_gan1 <-
  runAutoencoderDenoising(
    object = Colon_data_gan1, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    #layers = c(128, 64, 32), 
    dropout = 0.1
  )
getExpressionMatrixNames(object =Colon_data_gan1)
# active expression matrix after denoising
getActiveMatrixName(object = Colon_data_gan1)
###############################################################
# discard the one added above, to create new one
#Colon_data_gan1 <- discardFeatures(object = Colon_data_gan1, feature_names = "histology")
##################################
Colon_data_gan1 <- createSpatialSegmentation(Colon_data_gan1)
#################################################
Colon_data_gan1@fdata$gan1$histology<-ifelse(Colon_data_gan1@fdata$gan1$histology=="tumor","tumor","non-tumor")
#install.packages("Cairo",force=T)"#f3a546", "#cf3e53",
plotSurface(object =Colon_data_gan1, color_by = "histology", pt_clrp= "turbo",clrp_adjust = c("non-tumor" ="#00a2b3", "tumor" =  "#cf3e53"))
################################
ggsave(file = file.path(Figure3,"LM1_plotCnvHeatmap_histology.pdf"),width =11.97,height = 6.5)####
# names of grouping variables
getGroupingOptions(object = Colon_data_gan1)

# only histology
plotSurface(object =Colon_data_gan1, pt_alpha = 0,bcsp_rm=T)
#save(stmetadata,file = "data/stmetadata.Rdata")
#load("data/stmetadata.Rdata")
########################################################

# 添加其他变量stmetadata
#####################
stmetadata<-st_sce@meta.data;stmetadata$barcodes<-stmetadata$rownames
LM1<-stmetadata[stmetadata$orig.ident=="LM_1",]
LM1$leiden<-factor(LM1$leiden)
#LM1<- lapply(LM1, function(x) x)
# 尝试再次运行 addFeatures 函数
#####################
Colon_data_gan1 <- 
  addFeatures(
    object = Colon_data_gan1, 
    feature_df =LM1, 
    overwrite = TRUE
  )

###########CNV analysis 
Colon_data_gan1 <-
  runCnvAnalysis(
    object = Colon_data_gan1,
    directory_cnv_folder = "data/cnv-results", # example directory
    cnv_prefix = "Chr"
  )
cnv_results <- getCnvResults(Colon_data_gan1)

names(cnv_results)
##
Colon_data_gan1<-readRDS("data/spata_objects/Colon_spatia_gan1.RDS")
plotCnvHeatmap(object = Colon_data_gan1, across ="histology",clrp_adjust = c("non-tumor" ="#00a2b3", "tumor" =  "#cf3e53"))
ggsave(file = file.path(Figure3,"LM1_plotCnvHeatmap.pdf"),width =12.97,height = 6.5)####

################
plotCnvHeatmap(object = Colon_data_gan1, across ="st_celltype",
               clrp_adjust = c("Fibroblast" ="#54990F", "Tumor" =  "#D95F02",
                               "Hepatocytes"='#E41A1C',"Normal Epi"= "#3182BD",
                               "B/Plasma"='#FDB462',"Monocyte"="#00a2b3","Lamina propria"='#984EA3')
)
ggsave(file = file.path(Figure3,"LM1_plotCnvHeatmap_st_celltyp.pdf"),width =12.97,height = 6.5)####

###############################展示单基因
Colon_data_gan1 <- setActiveMatrix(object = Colon_data_gan1, mtr_name = "denoised")
############################

###############################展示多基因
plotSurfaceComparison(
  object = Colon_data_gan1, 
  color_by = c("COL1A1", "SPP1", "CD8A"), 
  pt_clrsp = "Greens 3", 
  display_image = TRUE, 
  smooth = TRUE, 
  alpha_by = TRUE
) 
ggsave(file = file.path(Figure3,"LM1_SurfaceComparison_genes.pdf"),width =8.97,height = 6.5)####



