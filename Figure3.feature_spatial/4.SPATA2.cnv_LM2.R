library(SPATA2)
library(SPATAData)
library(tidyverse)
###
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
####

launchSpataData()
SPATA2::sourceDataFrame()
###
Colon_data_gan2 <-
  initiateSpataObject_10X(
    directory_10X = '/home/data/gaoyuzhen/Projects/LM_spatial/spatial/LM_2', # the directory from which to load the data
    sample_name = "gan2"
  )
getDirectoryInstructions(object = Colon_data_gan2)
# set/change the current default directory
Colon_data_gan2 <- setSpataDir(Colon_data_gan2, dir = "data/spata_objects/Colon_spatia_gan2.RDS")
# quickly save the `spata2` object under the default directory 
saveSpataObject(object = Colon_data_gan2)
##############################
Colon_data_gan2 <- setActiveMatrix(object = Colon_data_gan2, mtr_name = "scaled")
##机器学习去除噪声 
Colon_data_gan2 <-
  runAutoencoderDenoising(
    object = Colon_data_gan2, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    #layers = c(128, 64, 32), 
    dropout = 0.1
  )
getExpressionMatrixNames(object =Colon_data_gan2)
# active expression matrix after denoising
getActiveMatrixName(object = Colon_data_gan2)
###############################################################
# discard the one added above, to create new one
#Colon_data_gan2 <- discardFeatures(object = Colon_data_gan2, feature_names = "histology")
##################################
Colon_data_gan2 <- createSpatialSegmentation(Colon_data_gan2)
#################################################
Colon_data_gan2@fdata$gan2$histology<-ifelse(Colon_data_gan2@fdata$gan2$histology=="tumor","tumor","non-tumor")
#install.packages("Cairo",force=T)"#f3a546", "#cf3e53",
plotSurface(object =Colon_data_gan2, color_by = "histology", pt_clrp= "turbo",clrp_adjust = c("non-tumor" ="#00a2b3", "tumor" =  "#cf3e53"))
################################
ggsave(file = file.path(Figure3,"LM2_plotCnvHeatmap_histology.pdf"),width =11.97,height = 6.5)####
# names of grouping variables
getGroupingOptions(object = Colon_data_gan2)
###################


# only histology
plotSurface(object =Colon_data_gan2, pt_alpha = 0,bcsp_rm=T)
#save(stmetadata,file = "data/stmetadata.Rdata")
#load("data/stmetadata.Rdata")
########################################################

# 添加其他变量stmetadata
#####################
stmetadata<-st_sce@meta.data;stmetadata$barcodes<-stmetadata$rownames
LM2<-stmetadata[stmetadata$orig.ident=="LM_2",]
LM2$leiden<-factor(LM2$leiden)
colnames(LM2)
#LM2<- lapply(LM2, function(x) x)
# 尝试再次运行 addFeatures 函数
#####################
Colon_data_gan2 <- 
  addFeatures(
    object = Colon_data_gan2, 
    feature_df =LM2, 
    overwrite = TRUE
  )

###########CNV analysis 
Colon_data_gan2 <-
  runCnvAnalysis(
    object = Colon_data_gan2,
    directory_cnv_folder = "data/SPATA2_cnv-results", # example directory
    cnv_prefix = "Chr"
  )
cnv_results <- getCnvResults(Colon_data_gan2)

names(cnv_results)
##
Colon_data_gan2<-readRDS("data/spata_objects/Colon_spatia_gan2.RDS")
plotCnvHeatmap(object = Colon_data_gan2, across ="histology",clrp_adjust = c("non-tumor" ="#00a2b3", "tumor" =  "#cf3e53"))
ggsave(file = file.path(Figure3,"LM2_plotCnvHeatmap.pdf"),width =12.97,height = 6.5)####


c("#54990F", "grey", '#E41A1C', '#984EA3','#FDB462', "#3182BD")

table(LM2$"st_celltype")


plotCnvHeatmap(object = Colon_data_gan2, across ="st_celltype",
               clrp_adjust = c("Fibroblast" ="#54990F", "Tumor" =  "#D95F02",
                               "Hepatocytes"='#E41A1C',"Normal Epi"= "#3182BD",
                               "B/Plasma"='#FDB462',"Monocyte"="#00a2b3","Lamina propria"='#984EA3')
               )
ggsave(file = file.path(Figure3,"LM2_plotCnvHeatmap_st_celltyp.pdf"),width =12.97,height = 6.5)####


###############################展示单基因
Colon_data_gan2 <- setActiveMatrix(object = Colon_data_gan2, mtr_name = "denoised")
############################
plotSurface(
  object = Colon_data_gan2, 
  color_by ="COL1A1",smooth = TRUE, 
  pt_size = 1.8,na_rm=F
) +  ggpLayerThemeCoords(unit = "px")+
  labs(title = "Denoised")


###############################展示多基因
plotSurfaceComparison(
  object = Colon_data_gan2, 
  color_by = c("COL1A1", "SPP1", "CD8A"), 
  pt_clrsp = "Greens 3", 
  display_image = TRUE, 
  smooth = TRUE, 
  alpha_by = TRUE
) 
ggsave(file = file.path(Figure3,"LM2_SurfaceComparison_genes.pdf"),width =8.97,height = 6.5)####

Colon_data_gan2 <- setActiveMatrix(object = Colon_data_gan2, mtr_name = "scaled")

a<-plotSurface(
  object = Colon_data_gan2, 
  color_by = "GPR182", 
  pt_size = 1.8
) + 
  labs(title = "Scaled (not denoised)") + 
  legendNone()

Colon_data_gan2 <- setActiveMatrix(object = Colon_data_gan2, mtr_name = "denoised")

b<-plotSurface(
  object = Colon_data_gan2, 
  color_by = "GPR182",
  pt_size = 1.8
) + 
  labs(title = "Denoised")
a+b
