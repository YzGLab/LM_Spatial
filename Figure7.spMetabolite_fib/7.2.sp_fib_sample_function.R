
setwd("/home/data/gaoyuzhen/Projects/LM_spatial")


sce<-readRDS("data/metabolits_scp_st_sp_sce.RDS")
metadata<-sce@meta.data
metabolites_obj<-sce@assays$metabolits@counts %>% data.frame(.,check.names = FALSE)
#########################
# Clear environment
#rm(list=ls())
# Load necessary library
library(MetaboAnalystR)

# Define sample names
sample_names <- c("LM_1", "LM_2", "C_1", "C_2")
feat.rank.mat.list<-list()
# Loop through each sample
for (sample in sample_names[1:2]) {

  print(sample)
  # Filter metadata for the current sample
  metadata_sample <- metadata[metadata$orig.ident == sample, ]
  
  # Create a label for Fibroblast cells
  Lable <- cbind(Sample=rownames(metadata_sample), 
                 Lable=ifelse(metadata_sample$st_celltype == "Fibroblast", 1, 0))
  
  # Subset metabolites for the current sample and transpose
  Metabolites_sample <- t(metabolites_obj[, colnames(metabolites_obj) %in% rownames(metadata_sample)])
  Metabolites_sample<-cbind(Sample=rownames(Metabolites_sample),Metabolites_sample)
  # Combine labels and metabolite data
  Metadat <- merge(Lable, Metabolites_sample,by="Sample",all = FALSE)
  
  # Create a directory for the current sample if it doesn't exist
  sample_dir <- paste0("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure8_feature_fib_sp/results_diff_meta_", sample)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir)
  }
  #setwd_paths<-paste0("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure8_feature_fib_sp/",sample_dir)
  setwd(sample_dir)
  # Write the combined data to a CSV file
  write.csv(Metadat, file=paste0(sample_dir, "/metadata.csv"), row.names = FALSE)
  rm(mSet)
  # Initialize MetaboAnalystR
  mSet <- InitDataObjects("conc", "stat", FALSE)
  mSet <- Read.TextData(mSet, paste0(sample_dir, "/metadata.csv"), "rowu", "disc")
  mSet <- SanityCheckData(mSet)
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet);
  mSet <- RemoveMissingPercent(mSet, percent = 0.2)
  mSet <- ImputeMissingVar(mSet, method = "min")
  mSet<-SanityCheckData(mSet)
  mSet<-FilterVariable(mSet, "F", 25, "none", -1, "mean", 0)
  mSet <- PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
  mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
  mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
  

  # Perform statistical analyses
  ##########
  mSet<-FC.Anal(mSet, 1.5, 1, FALSE)##1代表阳性组
  mSet<-PlotFC(mSet, "fc_1_", "pdf", 72, width=NA)
  #tttest
  mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
  #mSet<-PlotTT(mSet, "tt_1_", "PDF", 72, width=NA)
  #Volcano
  mSet<-Volcano.Anal(mSet, FALSE, 1.5, 1, F, 0.1, TRUE, "raw")#1代表阳性组
  mSet<-PlotVolcano(mSet, "volcano_1_",1, 0, "PDF", 72, width=NA, -1)
  ###OPLS 
  mSet<-OPLSR.Anal(mSet, reg=TRUE)
  mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", "PDF", 72, width=NA, 1,2,0.95,0,0, "na")
  mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "PDF", 72, width=NA);
  mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "PDF", 72, width=NA, "vip", "tscore", 20,FALSE)
  mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", "PDF", 72, width=NA)
  mSet<-PlotOPLS2DScore(mSet, "opls_score2d_1_", "PDF", 72, width=NA, 1,2,0.95,0,1, "na")
  mSet<-UpdateOPLS.Splot(mSet, "all");
  mSet<-PlotOPLS.Splot(mSet, "opls_splot_1_", "all", "PDF", 72, width=NA);
  mSet<-PlotOPLS.Imp(mSet, "opls_imp_1_", "PDF", 72, width=NA, "vip", "tscore", 20,FALSE)
  ##PCA
  mSet<-PCA.Anal(mSet)
  mSet<-PlotPCAPairSummary(mSet, "pca_pair_1_", "PDF", 72, width=NA, 5)
  mSet<-PlotPCAScree(mSet, "pca_scree_1_", "PDF",72, width=NA, 5)
  mSet<-PlotPCA2DScore(mSet, "pca_score2d_1_", "PDF",72, width=NA, 1,2,0.95,0,0)
  
  # Extract significant features
  feat.rank.mat <- cbind(
    orthoVipVn = mSet$analSet$oplsda$orthoVipVn,
    foldchange = mSet$analSet$fc$fc.all,
    pvalue = mSet$analSet$tt$p.value
  )
  feat.rank.mat <- data.frame(feat.rank.mat)
  feat.rank.mat$metabolits<-gsub("X","",rownames(feat.rank.mat))
  feat.rank <- feat.rank.mat[
    feat.rank.mat[, 1] > 1 & 
      abs(log2(feat.rank.mat[, 2])) > 0 & 
      feat.rank.mat[, 3] < 0.05, 
  ]
  feat.rank.mat.list[[sample]]<-feat.rank.mat
  # Save significant features to CSV
  write.csv(feat.rank.mat, file=paste0(sample_dir, "/all_features.csv"), row.names = TRUE)
  write.csv(feat.rank, file=paste0(sample_dir, "/significant_features.csv"), row.names = TRUE)
}

# Print the current working directory
print(getwd())
###################


# Load necessary library
library(MetaboAnalystR)
# Define sample names
dataset_names <- c("CRC", "LM")

# Loop through each sample
for (sample in dataset_names) {
  
  print(sample)
  # Filter metadata for the current sample
  metadata_sample <- metadata[metadata$set == sample, ]
  
  # Create a label for Fibroblast cells
  Lable <- cbind(Sample=rownames(metadata_sample), 
                 Lable=ifelse(metadata_sample$st_celltype == "Fibroblast", 1, 0))
  
  # Subset metabolites for the current sample and transpose
  Metabolites_sample <- t(metabolites_obj[, colnames(metabolites_obj) %in% rownames(metadata_sample)])
  
  # Combine labels and metabolite data
  Metadat <- cbind(Lable, Metabolites_sample)
  
  # Create a directory for the current sample if it doesn't exist
  sample_dir <- paste0("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure8_feature_fib_sp/results_diff_meta_", sample)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir)
  }
  #setwd_paths<-paste0("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure8_feature_fib_sp/",sample_dir)
  setwd(sample_dir)
  # Write the combined data to a CSV file
  write.csv(Metadat, file=paste0(sample_dir, "/metadata.csv"), row.names = FALSE)
  rm(mSet)
  # Initialize MetaboAnalystR
  mSet <- InitDataObjects("conc", "stat", FALSE)
  mSet <- Read.TextData(mSet, paste0(sample_dir, "/metadata.csv"), "rowu", "disc")
  mSet <- SanityCheckData(mSet)
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet);
  mSet <- RemoveMissingPercent(mSet, percent = 0.2)
  mSet <- ImputeMissingVar(mSet, method = "min")
  mSet<-SanityCheckData(mSet)
  mSet<-FilterVariable(mSet, "F", 25, "none", -1, "mean", 0)
  mSet <- PreparePrenormData(mSet)
  mSet <- Normalization(mSet, "NULL", "NULL", "NULL", ratio = FALSE, ratioNum = 20)
  mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
  mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
  
  
  # Perform statistical analyses
  ##########
  mSet<-FC.Anal(mSet, 1.5, 1, FALSE)
  mSet<-PlotFC(mSet, "fc_1_", "pdf", 72, width=NA)
  #tttest
  mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "fdr", FALSE)
  mSet<-PlotTT(mSet, "tt_1_", "PDF", 72, width=NA)
  #Volcano
  mSet<-Volcano.Anal(mSet, FALSE, 1.5, 1, F, 0.1, TRUE, "raw")
  mSet<-PlotVolcano(mSet, "volcano_1_",1, 0, "PDF", 72, width=NA, -1)
  ###OPLS 
  mSet<-OPLSR.Anal(mSet, reg=TRUE)
  mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", "PDF", 72, width=NA, 1,2,0.95,0,0, "na")
  mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "PDF", 72, width=NA);
  mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "PDF", 72, width=NA, "vip", "tscore", 20,FALSE)
  mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", "PDF", 72, width=NA)
  mSet<-PlotOPLS2DScore(mSet, "opls_score2d_1_", "PDF", 72, width=NA, 1,2,0.95,0,1, "na")
  mSet<-UpdateOPLS.Splot(mSet, "all");
  mSet<-PlotOPLS.Splot(mSet, "opls_splot_1_", "all", "PDF", 72, width=NA);
  mSet<-PlotOPLS.Imp(mSet, "opls_imp_1_", "PDF", 72, width=NA, "vip", "tscore", 20,FALSE)
  ##PCA
  mSet<-PCA.Anal(mSet)
  mSet<-PlotPCAPairSummary(mSet, "pca_pair_1_", "PDF", 72, width=NA, 5)
  mSet<-PlotPCAScree(mSet, "pca_scree_1_", "PDF",72, width=NA, 5)
  mSet<-PlotPCA2DScore(mSet, "pca_score2d_1_", "PDF",72, width=NA, 1,2,0.95,0,0)
  
  # Extract significant features
  feat.rank.mat <- cbind(
    orthoVipVn = mSet$analSet$oplsda$orthoVipVn,
    foldchange = mSet$analSet$fc$fc.all,
    pvalue = mSet$analSet$tt$p.value
  )
  feat.rank.mat <- data.frame(feat.rank.mat)
  feat.rank.mat$metabolits<-gsub("X","",rownames(feat.rank.mat))
  feat.rank <- feat.rank.mat[
    feat.rank.mat[, 1] > 1 & 
      abs(log2(feat.rank.mat[, 2])) > 0 & 
      feat.rank.mat[, 3] < 0.001, 
  ]
  feat.rank.mat.list[[sample]]<-feat.rank.mat
  # Save significant features to CSV
  write.csv(feat.rank.mat, file=paste0(sample_dir, "/all_features.csv"), row.names = TRUE)
  write.csv(feat.rank, file=paste0(sample_dir, "/significant_features.csv"), row.names = TRUE)
  #save(feat.rank.mat.list,feat.rank,file=paste0(sample_dir, "/significant_features.Rdata"))
}



# Print the current working directory
print(getwd())
############

