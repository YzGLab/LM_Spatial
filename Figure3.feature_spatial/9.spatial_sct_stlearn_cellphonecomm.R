#####
C1_lr_summary <- read_csv("data/stlearn_spaital_comm/C1/C1_lr_summary.csv")
C2_lr_summary <- read_csv("data/stlearn_spaital_comm/C2/C2_lr_summary.csv")
LM2_lr_summary <- read_csv("data/stlearn_spaital_comm/LM2/LM2_lr_summary.csv")
LM1_lr_summary <- read_csv("data/stlearn_spaital_comm/LM1/LM1_lr_summary.csv")



data %>% 
  #filter(pvalue <0.05) %>% # 如果不想把p值大于0.05的放在图上，去掉最前面的#号
  ggplot(aes(correlation,forcats::fct_reorder(symbol,correlation))) +
  geom_segment(aes(xend=0,yend=symbol)) +
  geom_point(aes(col=pvalue,size=abs(correlation))) +
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3")) +
  #scale_color_viridis_c(begin = 0.5, end = 1) +
  scale_size_continuous(range =c(2,8))  +
  theme_minimal() +
  ylab(NULL)

df<-C2_lr_summary
colnames(df)[1]<-"CC"
df$ratio<-df$n_spots_sig/df$n_spots
ggplot(df[1:30,],aes(x=CC, y=ratio))+
  #=geom_hline(yintercept = seq(0, 10, 2.5),linetype = 2, color = "lightgray",size=1)+
  geom_line()+
  geom_segment(aes(x=CC,xend=CC,y=0,yend=ratio),color="lightgray",size = 1.5)+
  geom_point(size=3,aes(color=ratio))+
  #scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#41ab5d")) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  theme(axis.text.x = element_text(size=8,  hjust = 1,angle=45, colour = "black"))+
  labs(x="",y="n_spots_sig/n_spots")
  ylim(c(190,320))


  #####
  C1_lr_summary <- read_csv("data/stlearn_spaital_comm/C1/C1_lr_summary.csv")
  C2_lr_summary <- read_csv("data/stlearn_spaital_comm/C2/C2_lr_summary.csv")
  LM2_lr_summary <- read_csv("data/stlearn_spaital_comm/LM2/LM2_lr_summary.csv")
  LM1_lr_summary <- read_csv("data/stlearn_spaital_comm/LM1/LM1_lr_summary.csv")
  
  
  library(ggpubr)
  df<-C2_lr_summary
  colnames(df)[1]<-"CC"
  df$ratio<-df$n_spots_sig/df$n_spots
###########棒状图 可视化
  #########################
  ggdotchart(df[1:50,], x = "CC", y = "ratio",
             color = "ratio",                                # Color by groups
             #palette = c("#00AFBB",  "#FC4E07"), # Custom color palette"#E7B800",
             sorting = "descending",                        # Sort value in descending order
             #sorting = "ascending",
             add = "segments",                             # Add segments from y = 0 to dots
             dot.size ="ratio",   
             #group = "sig",  
             ggtheme = theme_pubr()                        # ggplot2 theme
  )
  
  
  library(ggpubr)
  library(readr)
  
  # Define a function to generate and save the plot
  generate_and_save_plot <- function(df, output_file) {
    colnames(df)[1] <- "CC"
    df$ratio <- df$n_spots_sig / df$n_spots
    
    p <- ggdotchart(df[1:50,], x = "CC", y = "ratio",
                    color = "ratio",                                # Color by groups
                    sorting = "descending",                        # Sort value in descending order
                    add = "segments",                              # Add segments from y = 0 to dots
                    dot.size = "ratio",   
                    ggtheme = theme_pubr()                         # ggplot2 theme
    )
    
    ggsave(output_file, plot = p)
  }
  
  # Load datasets
  C1_lr_summary <- read_csv("data/stlearn_spaital_comm/C1/C1_lr_summary.csv")
  C2_lr_summary <- read_csv("data/stlearn_spaital_comm/C2/C2_lr_summary.csv")
  LM2_lr_summary <- read_csv("data/stlearn_spaital_comm/LM2/LM2_lr_summary.csv")
  LM1_lr_summary <- read_csv("data/stlearn_spaital_comm/LM1/LM1_lr_summary.csv")
  
  # Apply the function to each dataset and save the plots
  generate_and_save_plot(C1_lr_summary, "C1_lr_summary_plot.pdf")
  generate_and_save_plot(C2_lr_summary, "C2_lr_summary_plot.pdf")
  generate_and_save_plot(LM2_lr_summary, "LM2_lr_summary_plot.pdf")
  generate_and_save_plot(LM1_lr_summary, "LM1_lr_summary_plot.pdf")
  
  
  