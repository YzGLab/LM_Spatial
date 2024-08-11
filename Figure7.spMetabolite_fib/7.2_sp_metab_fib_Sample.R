sce<-readRDS("data/metabolits_scp_st_sp_sce.RDS")

metadata<-sce@meta.data
metabolites_obj<-sce@assays$metabolits@counts %>% data.frame(.,check.names = FALSE)
#########################

sample_names <- c("LM_1", "LM_2", "C_1", "C_2")

sample<-sample_names[[1]]

metadata_sample<-metadata[metadata$orig.ident==sample,]

table(metadata_sample$st_celltype)



Lable<-cbind(Sample=rownames(metadata_sample),Lable=ifelse(metadata_sample$st_celltype=="Fibroblast",1,0))

Metabolites_sample<-metabolites_obj[,colnames(metabolites_obj) %in% rownames(metadata_sample)] %>% t()

Metadat<-cbind(Lable,Metabolites_sample)
  
write.csv(Metadat,file="metadata.csv",row.names = FALSE)
  
  
  


#???? ????ǰ??????
#rm(list=ls())
#remotes::install_github("xia-lab/MetaboAnalystR",force = TRUE)
getwd()
###
library(MetaboAnalystR)


mSet<-NULL

mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, "metadata.csv", "rowu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-SanityCheckData(mSet)
mSet<-RemoveMissingPercent(mSet, percent=0.2)
mSet<-ImputeMissingVar(mSet, method="min")
mSet<-SanityCheckData(mSet)
mSet<-FilterVariable(mSet, "none", "F", 25)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
#mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
#mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
#################################
######ֻ??Ҫ??????????
#####??FC????
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
#??tttest
mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
###??-OPLS ????
mSet<-OPLSR.Anal(mSet, reg=TRUE)
#mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
#mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width=NA, "vip", "tscore", 20,FALSE)
#mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", "png", 72, width=NA)
mSet<-PlotOPLS.Imp(mSet, "opls_imp_1_", "png", 72, width=NA, "vip", "tscore", 20,FALSE)
###????ɸѡ
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.1, TRUE, "raw")
mSet<-PlotVolcano(mSet, "volcano_0_",1, "pdf", width=NA)
######
feat.rank.mat<-c()
feat.rank.mat<-cbind(orthoVipVn=mSet$analSet$oplsda$orthoVipVn,
                     #logfoldchange=mSet$analSet$fc$fc.log,
                     foldchange=mSet$analSet$fc$fc.all,
                     pvalue=mSet$analSet$tt$p.value)
feat.rank.mat<-data.frame(feat.rank.mat)
feat.rank<-feat.rank.mat[feat.rank.mat[,1]>1 & abs(log2(feat.rank.mat[,2]))>0 & feat.rank.mat[,3]<0.001, ]



###########################################################????ѡ????????
#####??FC????
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
mSet<-PlotTT(mSet, "tt_0_", "png", 72, width=NA)
mSet<-Volcano.Anal(mSet, FALSE, 2.0, 0, F, 0.1, TRUE, "raw")
mSet<-PlotVolcano(mSet, "volcano_0_",1, "pdf", width=NA)
mSet<-PlotCorrHeatMap(mSet, "corr_0_", "png", 72, width=NA, "col", "pearson", "bwm", "overview", F, F, "0")

##??PCA ????
mSet<-PCA.Anal(mSet)
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "png", 72, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "png", 72, width=NA, 1,2);
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", "png", 72, width=NA, 1,2)
mSet<-PlotPCA3DLoading(mSet, "pca_loading3d_0_", "json", 1,2,3)

##?? PLS ????
mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLSPairSummary(mSet, "pls_pair_0_", "png", 72, width=NA, 5)
mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
mSet<-PlotPLSLoading(mSet, "pls_loading_0_", "png", 72, width=NA, 1, 2);
mSet<-PlotPLS3DLoading(mSet, "pls_loading3d_0_", "json", 1,2,3)

###??-PLSDA ????
#install.packages("pls")
library(pls)
mSet<-PLSDA.CV(mSet, "T",5, "Q2")
mSet<-PlotPLS.Classification(mSet, "pls_cv_0_", "png", 72, width=NA)
mSet<-PlotPLS.Imp(mSet, "pls_imp_0_", "png", 72, width=NA, "vip", "Comp. 1", 15,FALSE)

###??-SPLSR ????
mSet<-SPLSR.Anal(mSet, 5, 10, "same", "Mfold")
mSet<-PlotSPLSPairSummary(mSet, "spls_pair_0_", "png", 72, width=NA, 5)
mSet<-PlotSPLS2DScore(mSet, "spls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotSPLS3DScoreImg(mSet, "spls_score3d_0_", "png", 72, width=NA, 1,2,3, 40)
mSet<-PlotSPLSLoading(mSet, "spls_loading_0_", "png", 72, width=NA, 1,"overview");
mSet<-PlotSPLSDA.Classification(mSet, "spls_cv_0_", "png", 72, width=NA)
mSet<-PlotSPLS3DLoading(mSet, "spls_loading3d_0_", "json", 1,2,3)


###??-OPLS ????
mSet<-OPLSR.Anal(mSet, reg=TRUE)
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width=NA, "vip", "tscore", 15,FALSE)
mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", "png", 72, width=NA)
mSet<-PlotOPLS.Imp(mSet, "opls_imp_1_", "png", 72, width=NA, "vip", "tscore", 15,FALSE)


##Advanced Significance Analysis
###??-SAM ????
mSet<-SAM.Anal(mSet, "d.stat", FALSE, TRUE, 0.0, "sam_imp_0_")
mSet<-PlotSAM.Cmpd(mSet, "sam_imp_0_", "png", 72, width=NA)
mSet<-PlotSAM.FDR(mSet, "sam_view_0_", "png", 72, width=NA)

###??EBAM
mSet<-EBAM.Init(mSet, FALSE, TRUE, FALSE, -99.0, 0.9, "ebam_view_0_", "ebam_imp_0_")
mSet<-PlotEBAM.Cmpd(mSet, "ebam_imp_0_", "png", 72, width=NA)


##Cluster Analysis
###????ͼ????ͼ
mSet<-PlotHCTree(mSet, "tree_0_", "png", 72, width=NA, "euclidean", "ward.D")
mSet<-PlotHeatMap(mSet, "heatmap_0_", "png", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview", T, T, NULL, T, F)
###??Kmeans????
mSet<-Kmeans.Anal(mSet, 3)
mSet<-PlotKmeans(mSet, "km_0_", "png", 72, width=NA, "default", "F")
mSet<-PlotClustPCA(mSet, "km_pca_0_", "png", 72, width=NA, "default", "km", "F")

###??SOM????
mSet<-SOM.Anal(mSet, 1,3,"linear","gaussian")
mSet<-PlotSOM(mSet, "som_0_", "png", 72, width=NA, "default", "T")
mSet<-PlotClustPCA(mSet, "som_pca_0_", "png", 72, width=NA, "default", "som", "F")

####Classification & Feature Selection
###??RF????
mSet<-RF.Anal(mSet, 500,7,1)
mSet<-PlotRF.Classify(mSet, "rf_cls_0_", "png", 72, width=NA)
mSet<-PlotRF.VIP(mSet, "rf_imp_0_", "png", 72, width=NA)
mSet<-PlotRF.Outlier(mSet, "rf_outlier_0_", "png", 72, width=NA)

###??RSVM????
mSet<-RSVM.Anal(mSet, 10)
mSet<-PlotRSVM.Classification(mSet, "svm_cls_0_", "png", 72, width=NA)
mSet<-PlotRSVM.Cmpd(mSet, "svm_imp_0_", "png", 72, width=NA)

#####
#mSet<-SaveTransformedData(mSet)



###??????????ɸѡ????
write.table(feat.rank.mat,file="results/onefactor_results_total.xls",row.names = TRUE,sep = "\t")
write.table(feat.rank,file="results/onefactor_results.xls",row.names = TRUE,sep = "\t")

###???????? ??ת?ڶ???Ȼ??
one_feat.rank.mat<-feat.rank.mat
one_feat.rank.mat<-cbind(names=rownames(one_feat.rank.mat),one_feat.rank.mat)
save(one_feat.rank.mat,file="undeleteFile/one_feat.rank.mat.Rdata")
###############


