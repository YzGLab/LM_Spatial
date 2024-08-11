library(SPATA2)
library(SPATAData)
library(tidyverse)
###

launchSpataData()
SPATA2::sourceDataFrame()
###
Colon_data_C2 <-
  initiateSpataObject_10X(
    directory_10X = '/home/data/gaoyuzhen/Projects/LM_spatial/spatial/C_2', # the directory from which to load the data
    sample_name = "chang2"
  )

getDirectoryInstructions(object = Colon_data_C2)
# set/change the current default directory
Colon_data_C2 <- setSpataDir(Colon_data_C2, dir = "data/spata_objects/Colon_spatia_chang2.RDS")
# quickly save the `spata2` object under the default directory 
saveSpataObject(object = Colon_data_C2)
##############################
Colon_data_C2 <- setActiveMatrix(object = Colon_data_C2, mtr_name = "scaled")
##机器学习去除噪声 
Colon_data_C2 <-
  runAutoencoderDenoising(
    object = Colon_data_C2, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    #layers = c(128, 64, 32), 
    dropout = 0.1
  )
getExpressionMatrixNames(object =Colon_data_C2)
# active expression matrix after denoising
getActiveMatrixName(object = Colon_data_C2)
###############################################################
# discard the one added above, to create new one
#Colon_data_C2 <- discardFeatures(object = Colon_data_C2, feature_names = "histology")
##################################
Colon_data_C2 <- createSpatialSegmentation(Colon_data_C2)
#################################################
Colon_data_C2@fdata$chang2$histology<-ifelse(Colon_data_C2@fdata$chang2$histology=="tumor","tumor","non-tumor")
#install.packages("Cairo",force=T)


# 添加其他变量stmetadata
#####################
stmetadata<-st_sce@meta.data;stmetadata$barcodes<-stmetadata$rownames
C2<-stmetadata[stmetadata$orig.ident=="C_2",]
C2$leiden<-factor(C2$leiden)
colnames(C2)
# 尝试再次运行 addFeatures 函数
#####################
Colon_data_C2 <- 
  addFeatures(
    object = Colon_data_C2, 
    feature_df =C2, 
    overwrite = TRUE
  )


Colon_data_C2<-readRDS("data/spata_objects/Colon_spatia_chang2.RDS")
plotSurface(object =Colon_data_C2, color_by = "histology", pt_clrp= "turbo",clrp_adjust = c("non-tumor" ="#00a2b3", "tumor" =  "#cf3e53"))
ggsave(file = file.path(Figure3,"C2_plotCnvHeatmap_histology.pdf"),width =11.97,height = 6.5)####
# names of grouping variables
getGroupingOptions(object = Colon_data_C2)

################
plotCnvHeatmap(object = Colon_data_C2, across ="st_celltype",
               clrp_adjust = c("Fibroblast" ="#54990F", "Tumor" =  "#D95F02",
                               "Hepatocytes"='#E41A1C',"Normal Epi"= "#3182BD",
                               "B/Plasma"='#FDB462',"Monocyte"="#00a2b3","Lamina propria"='#984EA3')
)
ggsave(file = file.path(Figure3,"C2_plotCnvHeatmap_st_celltyp.pdf"),width =12.97,height = 6.5)####

# only histology
plotSurface(object =Colon_data_C2, pt_alpha = 0,bcsp_rm=T)
#save(stmetadata,file = "data/stmetadata.Rdata")
#load("data/stmetadata.Rdata")
########################################################

# 添加其他变量stmetadata
#####################
stmetadata<-st_sce@meta.data;stmetadata$barcodes<-stmetadata$rownames
C2<-stmetadata[stmetadata$orig.ident=="C_2",]
C2$leiden<-factor(C2$leiden)
#LM2<- lapply(LM2, function(x) x)
# 尝试再次运行 addFeatures 函数
#####################
Colon_data_C2 <- 
  addFeatures(
    object = Colon_data_C2, 
    feature_df =C2, 
    overwrite = TRUE
  )

###########CNV analysis 
Colon_data_C2 <-
  runCnvAnalysis(
    object = Colon_data_C2,
    directory_cnv_folder = "data/cnv-results/C2", # example directory
    cnv_prefix = "Chr"
  )
cnv_results <- getCnvResults(Colon_data_C2)

names(cnv_results)
##
plotCnvHeatmap(object = Colon_data_C2, across ="histology",clrp_adjust = c("non-tumor" ="#00a2b3", "tumor" =  "#cf3e53"))
ggsave(file = file.path(Figure3,"C2_plotCnvHeatmap.pdf"),width =11.97,height = 6.5)####


###############################展示单基因
Colon_data_C2 <- setActiveMatrix(object = Colon_data_C2, mtr_name = "denoised")
############################
a<-plotSurface(
  object = Colon_data_C2, 
  color_by ="GPR171",smooth = TRUE, 
  pt_size = 1.8,na_rm=F
) +  ggpLayerThemeCoords(unit = "px")+
  labs(title = "Denoised")

b<-plotSurface(
  object = Colon_data_gan2, 
  color_by ="GPR171",smooth = TRUE, 
  pt_size = 1.8,na_rm=F
) +  ggpLayerThemeCoords(unit = "px")+
  labs(title = "Denoised")
a+b
ggsave(file = file.path(Figure3,"GPR171_expression_C2_lM2.pdf"),width =11.97,height = 6.5)####



###############################展示多基因
plotSurfaceComparison(
  object = Colon_data_C2, 
  color_by = c("COL1A1", "SPP1", "CD8A"), 
  pt_clrsp = "Greens 3", 
  display_image = TRUE, 
  smooth = TRUE, 
  alpha_by = TRUE
) 

ggsave(file = file.path(Figure3,"C2_SurfaceComparison_genes.pdf"),width =8.97,height = 6.5)####

Colon_data_C2 <- setActiveMatrix(object = Colon_data_C2, mtr_name = "scaled")

a<-plotSurface(
  object = Colon_data_C2, 
  color_by = "GPR182", 
  pt_size = 1.8
) + 
  labs(title = "Scaled (not denoised)") + 
  legendNone()

Colon_data_C2 <- setActiveMatrix(object = Colon_data_C2, mtr_name = "denoised")

b<-plotSurface(
  object = Colon_data_C2, 
  color_by = "GPR182",
  pt_size = 1.8
) + 
  labs(title = "Denoised")
a+b
