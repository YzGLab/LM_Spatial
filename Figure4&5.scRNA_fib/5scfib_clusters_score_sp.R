####
library(SPATA2)
library(SPATAData)
library(tidyverse)
###
Colon_data_C2<-readRDS("data/spata_objects/Colon_spatia_chang2.RDS")
Colon_data_gan2<-readRDS("data/spata_objects/Colon_spatia_gan2.RDS")

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
##############
###
st_metadata<-st_sce.copy@meta.data

st_metadata_c2<-st_metadata[st_metadata$orig.ident=="C_2",]
rownames(st_metadata_c2)<-gsub("C_2_","",rownames(st_metadata_c2))
#rownames(st_metadata_c2)<-str_split(rownames(st_metadata_c2),"-") %>% sapply("[[",1)
colnames(st_metadata_c2)
st_c2<-Colon_data_C2@fdata$chang2[,1]
st_metadata_c2[st_metadata_c2$st_celltype!="Fibroblast",c(32:36)]<-0
Colon_data_C2@fdata$chang2<-cbind(Colon_data_C2@fdata$chang2,st_metadata_c2[,c(32:36)])
Colon_data_C2@fdata$chang2
plotSurface(
  object = Colon_data_C2, 
  color_by ="Fib(C4)-CXCL8",
  #smooth = TRUE, 
  pt_size = 1.8,na_rm=F
) +  ggpLayerThemeCoords(unit = "px")+
  labs(title = "Denoised")
ggsave(file = file.path(Figure4,"C2_SurfaceComparison_Fib(C4)_CXCL8_special.pdf"),width =9.3,height = 11.5)####

table(st_metadata$orig.ident)
st_metadata_gan2<-st_metadata[st_metadata$orig.ident=="LM_2",]
rownames(st_metadata_gan2)<-gsub("LM_2_","",rownames(st_metadata_gan2))
#rownames(st_metadata_c2)<-str_split(rownames(st_metadata_c2),"-") %>% sapply("[[",1)
colnames(st_metadata_gan2)
head(st_metadata_gan2)
st_metadata_gan2[st_metadata_gan2$st_celltype!="Fibroblast",c(32:36)]<-0
Colon_data_gan2@fdata$gan2[1,1]
Colon_data_gan2@fdata$gan2<-cbind(Colon_data_gan2@fdata$gan2,st_metadata_gan2[,c(32:36)])
Colon_data_gan2@fdata$gan2

#####################
#stmetadata<-st_sce@meta.data;stmetadata$barcodes<-stmetadata$rownames
#LM2<-stmetadata[stmetadata$orig.ident=="LM_2",]
#LM2$leiden<-factor(LM2$leiden)
#colnames(LM2)
##LM2<- lapply(LM2, function(x) x)
# 尝试再次运行 addFeatures 函数
######################
#Colon_data_gan2 <- 
#  addFeatures(
 #   object = Colon_data_gan2, 
  #  feature_df =LM2, 
  #  overwrite = TRUE
 # )




plotSurface(
  object = Colon_data_gan2, 
  color_by ="Fib(C4)-CXCL8",
  #smooth = TRUE, 
  pt_size = 1.8,na_rm=F
) +  ggpLayerThemeCoords(unit = "px")+
  labs(title = "Denoised")
ggsave(file = file.path(Figure4,"LM2_SurfaceComparison_Fib(C4)_CXCL8_special.pdf"),width =9.3,height = 11.5)####
