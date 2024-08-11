###full data with metabolits
library(readxl)
Allsample_pos <- read.delim2("spatial_metabolism/Metabolism_martix/Allsample-pos-sample_level.xls")
Allsample_neg <- read.delim2("spatial_metabolism/Metabolism_martix/Allsample-neg-sample_level.xls")
##"C_1.ALL.pos"  "C_2.ALL.pos"  "LM_1.ALL.pos" "LM_2.ALL.pos"  前两个样本 vs 后两个样本，计算fold change，筛选c组1.2,0.5，和LM组前一百代谢
# 安装并加载必要的包
#if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)
Allsample_pos[,c("C_1.ALL.pos" , "C_2.ALL.pos",  
                 "LM_1.ALL.pos" ,"LM_2.ALL.pos")]<-apply(Allsample_pos[,c("C_1.ALL.pos" , "C_2.ALL.pos", 
                                                                                                          "LM_1.ALL.pos" ,"LM_2.ALL.pos")],
                                                                                         2,
                                                                                         as.numeric)

# 假设你的数据存储在 Allsample_pos 数据框中
# 确保列名没有重复
colnames(Allsample_pos) <- make.names(colnames(Allsample_pos), unique = TRUE)

# 提取数值列的名称
C_columns <- c("C_1.ALL.pos", "C_2.ALL.pos")
LM_columns <- c("LM_1.ALL.pos", "LM_2.ALL.pos")

# 计算 fold change
# 计算 fold change 并过滤均值小于 5 的代谢物
Allsample_pos_fc <- Allsample_pos %>%
  rowwise() %>%
  mutate(
    C_mean = mean(c_across(all_of(C_columns)), na.rm = TRUE),
    LM_mean = mean(c_across(all_of(LM_columns)), na.rm = TRUE),
    fold_change = LM_mean/C_mean 
  ) %>%
  filter(C_mean >= 1 & LM_mean >= 1) %>%
  ungroup()

# 根据 fold change 对 C 组排序并筛选前 100 个代谢物
C_top100 <- Allsample_pos_fc %>%
  arrange(desc(fold_change)) %>%
  filter(fold_change <1) 

# 根据 fold change 对 LM 组排序并筛选前 100 个代谢物
LM_top100 <- Allsample_pos_fc %>%
  arrange(fold_change) %>%
  filter(fold_change >3) 

# 创建一个新的数据框，合并 C_top100 和 LM_top100
merged_top100 <- bind_rows(
  C_top100 %>% mutate(Group = "C_top100"),
  LM_top100 %>% mutate(Group = "LM_top100")
)

# 输出合并后的表格
print(merged_top100)

# 如果需要将结果写入 CSV 文件
write.csv(merged_top100,file.path(Figure7, "merged_metabolits_patient_different.csv"), row.names = FALSE)

##########################

Allsample_neg[,c("C_1.ALL.neg" , "C_2.ALL.neg",  
                 "LM_1.ALL.neg" ,"LM_2.ALL.neg")]<-apply(Allsample_neg[,c("C_1.ALL.neg" , "C_2.ALL.neg", 
                                                                          "LM_1.ALL.neg" ,"LM_2.ALL.neg")],
                                                         2,
                                                         as.numeric)


# 提取数值列的名称
C_columns <- c("C_1.ALL.neg", "C_2.ALL.neg")
LM_columns <- c("LM_1.ALL.neg", "LM_2.ALL.neg")

# 计算 fold change
# 计算 fold change 并过滤均值小于 5 的代谢物
Allsample_neg_fc <- Allsample_neg %>%
  rowwise() %>%
  mutate(
    C_mean = mean(c_across(all_of(C_columns)), na.rm = TRUE),
    LM_mean = mean(c_across(all_of(LM_columns)), na.rm = TRUE),
    fold_change = LM_mean/C_mean 
  ) %>%
  filter(C_mean >= 1 & LM_mean >= 1) %>%
  ungroup()

# 根据 fold change 对 C 组排序并筛选前 100 个代谢物
C_top100 <- Allsample_neg_fc %>%
  arrange(desc(fold_change)) %>%
  filter(fold_change <0.8) 

# 根据 fold change 对 LM 组排序并筛选前 100 个代谢物
LM_top100 <- Allsample_neg_fc %>%
  arrange(fold_change) %>%
  filter(fold_change >2) 

# 创建一个新的数据框，合并 C_top100 和 LM_top100
merged_top100_neg <- bind_rows(
  C_top100 %>% mutate(Group = "C_top100"),
  LM_top100 %>% mutate(Group = "LM_top100")
)

# 输出合并后的表格
print(merged_top100_neg)
# 如果需要将结果写入 CSV 文件
write.csv(merged_top100_neg,file.path(Figure7, "merged_metabolits_patient_different_neg.csv"), row.names = FALSE)############导入在线网站进行 通路富集分析


###根据fdr 画火山图
# 绘制火山图
library(ggplot2)
Allsample_neg_fc$sig<-ifelse(Allsample_neg_fc$fold_change>2,"LM_high",ifelse(Allsample_neg_fc$fold_change<0.8,"CRC_high","Not"))
Allsample_neg_fc$sig<-factor(Allsample_neg_fc$sig,levels=c("LM_high","Not","CRC_high"))
volcano_plot <- ggplot(Allsample_neg_fc, aes(x = log2(fold_change), y = log10(C_mean+LM_mean),col=sig)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_hline(yintercept =log2(2), linetype="dashed", color = "#999999")+  #添加y轴线
  geom_vline(xintercept= c(log2(0.8),log2(2)), linetype="dashed", color = "#999999")+
  scale_color_manual(values = c( "#E41A1C", "#54990F","#3182BD")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "Log2 Mean") +
  theme_bw() + #下面就是ggplot主题修饰
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1, colour = "black")) +
  theme(axis.title =element_text(size = 14),axis.text =element_text(size = 10, color = "black"),
        plot.title =element_text(hjust = 0.5, size = 16))

# 显示火山图
print(volcano_plot)
ggsave(file = file.path(Figure7,"sp_volcano_plot_topmetabolits_patientGroup.pdf"),width = 5.3,height =4.4)

Allsample_pos_fc$sig<-ifelse(Allsample_pos_fc$fold_change>2,"LM_high",ifelse(Allsample_pos_fc$fold_change<0.5,"CRC_high","Not"))
Allsample_pos_fc$sig<-factor(Allsample_pos_fc$sig,levels=c("LM_high","Not","CRC_high"))
volcano_plot <- ggplot(Allsample_pos_fc, aes(x = log2(fold_change), y = log10(C_mean+LM_mean),col=sig)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_hline(yintercept =log2(2), linetype="dashed", color = "#999999")+  #添加y轴线
  geom_vline(xintercept= c(-1,log2(2)), linetype="dashed", color = "#999999")+
  scale_color_manual(values = c( "#E41A1C", "#54990F","#3182BD")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "Log10 Mean") +
  theme_bw() + #下面就是ggplot主题修饰
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1, colour = "black")) +
  theme(axis.title =element_text(size = 14),axis.text =element_text(size = 10, color = "black"),
        plot.title =element_text(hjust = 0.5, size = 16))

# 显示火山图
print(volcano_plot)
ggsave(file = file.path(Figure7,"sp_volcano_plot_pos_fc_topmetabolits_patientGroup.pdf"),width = 5.3,height =4.4)
#
colnames(Allsample_pos_fc)[c(14:17)]<-c("C_1.ALL", "C_2.ALL","LM_1.ALL", "LM_2.ALL")
colnames(Allsample_neg_fc)[c(14:17)]<-c("C_1.ALL", "C_2.ALL","LM_1.ALL", "LM_2.ALL")
# 绘制火山图
library(ggplot2)
x<-rbind(Allsample_neg_fc,Allsample_pos_fc)
x$label<-x$Metabolites
x$sig<-factor(x$sig,levels=c("LM_high","Not","CRC_high"))
p1<-volcano_plot <- ggplot(x, aes(x = log2(fold_change), y = log10(C_mean+LM_mean),col=sig, label = label)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_hline(yintercept =log2(2), linetype="dashed", color = "#999999")+  #添加y轴线
  geom_vline(xintercept= c(-1,log2(2)), linetype="dashed", color = "#999999")+
  scale_color_manual(values = c( "#E41A1C", "#54990F","#3182BD")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "Log2 Mean") +
  theme_bw() + #下面就是ggplot主题修饰
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1, colour = "black")) +
  theme(axis.title =element_text(size = 14),axis.text =element_text(size = 10, color = "black"),
        plot.title =element_text(hjust = 0.5, size = 16))


selectgenes <- x[x$fold_change>10 | x$fold_change<0.2,]
selectgenes$label<-selectgenes$Metabolites
head(selectgenes)
selectgenes<-data.frame(selectgenes)
##
library(ggpubr)
library(ggrepel)
p2 <- p1 + 
  # 在感兴趣的基因外面画个黑色圈
  geom_point(data = selectgenes, alpha = 1, size = 4.6, shape = 1, 
             stroke = 1, #圈粗细
             color = "black") +
  
  # 显示感兴趣的基因的基因名
  geom_text_repel(data = selectgenes, 
                  colour="darkred", size = 5, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines"))
p2
# 显示火山图
print(volcano_plot)
ggsave(02,file = file.path(Figure7,"sp_volcano_plot_sample_ful_topmetabolits_patientGroup.pdf"),width = 8,height =5)
