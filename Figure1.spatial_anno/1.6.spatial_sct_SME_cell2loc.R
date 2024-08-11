##cell2location and stlearn SME clusters
####################################
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
stmetadata<-st_sce@meta.data
##################
C1<-stmetadata[stmetadata$orig.ident=='C_1',]
C2<-stmetadata[stmetadata$orig.ident=='C_2',]
LM1<-stmetadata[stmetadata$orig.ident=='LM_1',]
LM2<-stmetadata[stmetadata$orig.ident=='LM_2',]
#####################################################组合cell2location 信息四个样本。
###
C2_subcelltype_six_celltypes <- read_csv("data/stlearn_results/C2_subcelltype_six_celltypes.csv")
C1_subcelltype_six_celltypes <- read_csv("data/stlearn_results/C1_subcelltype_six_celltypes.csv")
LM1_subcelltype_six_celltypes <- read_csv("data/stlearn_results/LM1_subcelltype_six_celltypes.csv")
LM2_subcelltype_six_celltypes <- read_csv("data/stlearn_results/LM2_subcelltype_six_celltypes.csv")
colnames(C1)
C1<-cbind(C1_subcelltype_six_celltypes[,c(9:15)],C1) 
table(C1$leiden,C1$cell2loc)
C2<-cbind(C2_subcelltype_six_celltypes[,c(9:15)] ,C2) 
table(C2$leiden,C2$cell2loc)
LM1<-cbind(LM1_subcelltype_six_celltypes[,c(9:15)] ,LM1) 
table(LM1$leiden,LM1$cell2loc)
LM2<-cbind(LM2_subcelltype_six_celltypes[,c(9:15)] ,LM2) 
table(LM2$leiden,LM2$cell2loc)

stmetadata<-rbind(C1,C2,LM1,LM2)


st_sce<-AddMetaData(st_sce,stmetadata)

colnames(stmetadata)
table(stmetadata$orig.ident)
table(stmetadata$st_celltype)
######

write.csv(stmetadata,file = "data/stmetadata_full.csv")

save(stmetadata,file="data/stmetadata.Rdata")
saveRDS(st_sce,file="data/Colon_HC_spatial.RDS")


####
library(ggplot2)
library(reshape2)
colnames(LM2)
# Assuming your data is in a data frame called 'cell_data'
# Replace 'cell_data' with the name of your actual data frame
# Reshaping the data for ggplot
melted_data <- melt(LM2[,c(1:6,23)], id.vars = "leiden", variable.name = "Cell_Type", value.name = "Amount")
melted_data$Cell_Type<-factor(melted_data$Cell_Type)
melted_data$leiden<-factor(melted_data$leiden)
ggplot(melted_data, aes(x =leiden, y = Amount, fill =Cell_Type,colour= Cell_Type)) +
  scale_fill_discrete(name = "Cell_Type")+
  #geom_boxplot(alpha = 0.2,outlier.size=0)+
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values=c("#54990F", "grey", '#E41A1C', '#984EA3','#FDB462', "#3182BD")) +
  scale_color_manual(values=c("#54990F", "grey", '#E41A1C', '#984EA3','#FDB462', "#3182BD")) +
  theme_minimal() +
  # theme_classic2()+
  theme_bw()+
  labs(title = "Cell Type Amounts by Leiden Clustering",
       x = "Cluster",
       y = "Average abundence from Cell2location") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #coord_flip()
####
ggsave(file = file.path(Figure1,"st_sce_stlearn_immunecells_LM2.pdf"),width = 9.5,height = 4.03)
##
# Reshaping the data for ggplot
melted_data <- melt(LM1[,c(1:6,23)], id.vars = "leiden", variable.name = "Cell_Type", value.name = "Amount")
melted_data$Cell_Type<-factor(melted_data$Cell_Type)
melted_data$leiden<-factor(melted_data$leiden)
ggplot(melted_data, aes(x =leiden, y = Amount, fill =Cell_Type,colour= Cell_Type)) +
  scale_fill_discrete(name = "Cell_Type")+
  #geom_boxplot(alpha = 0.2,outlier.size=0)+
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values=c("#54990F", "grey", '#E41A1C', '#984EA3','#FDB462', "#3182BD")) +
  scale_color_manual(values=c("#54990F", "grey", '#E41A1C', '#984EA3','#FDB462', "#3182BD")) +
  theme_minimal() +
  # theme_classic2()+
  theme_bw()+
  labs(title = "Cell Type Amounts by Leiden Clustering",
       x = "Cluster",
       y = "Average abundence from Cell2location") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #coord_flip()
####
ggsave(file = file.path(Figure1,"st_sce_stlearn_immunecells_LM1.pdf"),width = 9.5,height = 4.03)


# Reshaping the data for ggplot
melted_data <- melt(C2[,c(1:6,23)], id.vars = "leiden", variable.name = "Cell_Type", value.name = "Amount")
melted_data$Cell_Type<-factor(melted_data$Cell_Type)
melted_data$leiden<-factor(melted_data$leiden)
ggplot(melted_data, aes(x =leiden, y = Amount, fill =Cell_Type,colour= Cell_Type)) +
  scale_fill_discrete(name = "Cell_Type")+
  #geom_boxplot(alpha = 0.2,outlier.size=0)+
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values=c("#54990F", "grey", '#E41A1C', '#984EA3','#FDB462', "#3182BD")) +
  scale_color_manual(values=c("#54990F", "grey", '#E41A1C', '#984EA3','#FDB462', "#3182BD")) +
  theme_minimal() +
  # theme_classic2()+
  ylim(c(0,40))+
  theme_bw()+
  labs(title = "Cell Type Amounts by Leiden Clustering",
       x = "Cluster",
       y = "Average abundence from Cell2location") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #coord_flip()
####
ggsave(file = file.path(Figure1,"st_sce_stlearn_immunecells_C2.pdf"),width = 9.5,height = 4.03)

# Reshaping the data for ggplot
melted_data <- melt(C1[,c(1:6,23)], id.vars = "leiden", variable.name = "Cell_Type", value.name = "Amount")
melted_data$Cell_Type<-factor(melted_data$Cell_Type)
melted_data$leiden<-factor(melted_data$leiden)
ggplot(melted_data, aes(x =leiden, y = Amount, fill =Cell_Type,colour= Cell_Type)) +
  scale_fill_discrete(name = "Cell_Type")+
  #geom_boxplot(alpha = 0.2,outlier.size=0)+
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values=c("#54990F", "grey", '#E41A1C', '#984EA3','#FDB462', "#3182BD")) +
  scale_color_manual(values=c("#54990F", "grey", '#E41A1C', '#984EA3','#FDB462', "#3182BD")) +
  theme_minimal() +
  # theme_classic2()+
  ylim(c(0,50))+
  theme_bw()+
  labs(title = "Cell Type Amounts by Leiden Clustering",
       x = "Cluster",
       y = "Average abundence from Cell2location") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #coord_flip()
####
ggsave(file = file.path(Figure1,"st_sce_stlearn_immunecells_C1.pdf"),width = 9.5,height = 4.03)




