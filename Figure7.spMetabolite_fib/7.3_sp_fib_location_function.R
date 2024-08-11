sce<-readRDS("data/metabolits_scp_st_sp_sce.RDS")
metadata<-sce@meta.data
metabolites_obj<-sce@assays$metabolits@counts %>% data.frame(.,check.names = FALSE)
#########################
# Clear environment
rm(list=ls())
setwd("/home/data/gaoyuzhen/Projects/LM_spatial")
#for LM vs CRC
# Load necessary library
# Load necessary library
library(MetaboAnalystR)

sample<-"LM_vs_CRC"
# Filter metadata for the current sample
metadata_sample <- metadata[metadata$st_celltype =="Fibroblast", ]

# Create a label for Fibroblast cells
Lable <- cbind(Sample=rownames(metadata_sample), 
               Lable=ifelse(metadata_sample$set == "LM", 1, 0))

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
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet)
mSet<-FilterVariable(mSet, "F", 25, "none", -1, "mean", 0)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)


# Perform statistical analyses
##########
mSet<-FC.Anal(mSet, 1.5, 1, FALSE)
mSet<-PlotFC(mSet, "fc_1_", "PDF", 72, width=NA)

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

feat.rank.mat.LM_fib<- data.frame(feat.rank.mat)
feat.rank.mat.LM_fib$metabolits<-gsub("X","",rownames(feat.rank.mat.LM_fib))
feat.rank <-feat.rank.mat.LM_fib[
  feat.rank.mat.LM_fib[, 1] > 1 & 
    log2(feat.rank.mat[, 2]) > 0.58 & 
    feat.rank.mat.LM_fib[, 3] < 0.001, 
]
feat.rank.mat.list[["LM_vs_CRC"]]<-feat.rank.mat.LM_fib

# Save significant features to CSV
write.csv(feat.rank.mat.LM_fib, file=paste0(sample_dir, "/all_features.csv"), row.names = TRUE)
write.csv(feat.rank, file=paste0(sample_dir, "/significant_features.csv"), row.names = TRUE)


save(feat.rank.mat.list,file=file.path(Figure8,"significant_features.Rdata"))



