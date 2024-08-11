###prognosis of cell of tumor 
#计算不同细胞的预后 ，并以meta-结果输出

load("~/Projects/CRCprojects/AllBulkdata_CRC/Colon_rectal_2021.Rdata")
load("~/Projects/PanCancerdata/Pancancerexp_split.Rdata")
#############
library(limma)
library(GSVA)
library(data.table)
library(survminer)
library(survival)
library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)
#save(sce.cell_type_markers,sce.markers.epi,sce.markers.epi.sig,file='data/spatial_scRNA_differetgenes_clusters.Rdata')
st_sce<-readRDS("data/Colon_HC_spatial.RDS")
##############
Idents(sub_sce)<-sub_sce$st_celltype
sce.markers1 <- FindAllMarkers(sub_sce[,sub_sce$set=="CRC"], 
                               only.pos = TRUE, 
                               min.pct = 0.5, logfc.threshold = 0.5)
sce.markers2 <- FindAllMarkers(sub_sce[,sub_sce$set=="LM"], 
                               only.pos = TRUE, 
                               min.pct = 0.5, logfc.threshold = 0.5)
#sce.markers.sig <- sce.markers %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)# top 10

sce.markers1$cluster<-paste0("CRC_",sce.markers1$cluster)
sce.markers2$cluster<-paste0("LM_",sce.markers2$cluster)
sce.markers<-rbind(sce.markers1,sce.markers2)################################
sce.cell_type_markers<-split(sce.markers$gene,sce.markers$cluster)###############合并两种差异


#sce.cell_type_markers
#sce.markers.epi.sig_st_sc<-sce.markers.sig_st_sc[1:12]
names(exp)
names(phe)
gmt.list<-sce.cell_type_markers
#for RFS
##
Scorelist<-list()
for(set in names(exp)){
  myexp<-exp[[set]]
  myphe<-phe[[set]]
  print(set)
  #myexp<- normalizeBetweenArrays(myexp)
  norm.expMat<-as.matrix(myexp)
  gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Gaussian",parallel.sz=60)###Array data,ssgsea.norm=TRUE
  #gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Poisson",parallel.sz=10)#RNA-Seq 
  gsva_es<-t(gsva_es)
  gsva_es<-cbind(rownames(gsva_es),gsva_es)
  Scorelist[[set]]<-gsva_es
  message(set,"GSVA is done!")
}
###for tcga cohort 
for(set in names(exp)[c(7,8)]){
  myexp<-exp[[set]]
  myphe<-phe[[set]]
  print(set)
  #myexp<- normalizeBetweenArrays(myexp)
  norm.expMat<-as.matrix(myexp)
  #gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Gaussian",parallel.sz=20)###Array data,ssgsea.norm=TRUE
  gsva_es <- gsva(norm.expMat,gmt.list,method="ssgsea",abs.ranking=F,kcdf="Poisson",parallel.sz=10)#RNA-Seq 
  gsva_es<-t(gsva_es)
  gsva_es<-cbind(rownames(gsva_es),gsva_es)
  Scorelist[[set]]<-gsva_es
  message(set,"GSVA is done!")
}

names(Scorelist)
save(Scorelist,file = "data/Scorelist_allCRC_st_sc.Rdata")
load("data/Scorelist_allCRC_st_sc.Rdata")
#for RFS
b<-c()
for (set in names(phe)){
  myphe<-phe[[set]]
  a<-ifelse("RFS" %in% colnames(myphe),"RFS","NO")
  a<-cbind(a,set)
  b<-rbind(b,a)
}##
b
set<-"TCGA_COAD"
#genes<-c("GPR171","CD8A")
colorectal_cox.sig<-c()
for(set in names(exp)[c(1,3:14,18:19)]){
  print(set)
  gsva_es<-Scorelist[[set]]
  colnames(gsva_es)<-c('geo_accession',names(gmt.list))
  myphe<-phe[[set]]
  #colnames(PANCAN.survival)
  #myexp<-exp[[set]]
  #mygenexp<-myexp[rownames(myexp) %in% genes,]
  #mygenexp<-cbind(geo_accession=colnames(mygenexp),t(mygenexp))
  #myphe<-phe[[set]]
  #myphe<-merge(mygenexp,myphe,by="geo_accession")
  #colnames(PANCAN.survival)
  rt<-merge(gsva_es,myphe,by="geo_accession")
  rt$RFS.time<-as.numeric(as.character(rt$RFS_time))
  rt$RFS<-as.numeric(as.character(rt$RFS))
  #rt[,"CD8A"]<-as.numeric(rt[,"CD8A"])
  rt<-rt[which(rt$RFS!='NA'),]
  rt<-rt[which(rt$RFS.time!='NA'),]
  if(nrow(rt)>10){
    for (sig in names(gmt.list)){
      rt[,sig]<-as.numeric(as.character(rt[,sig]))
      #rt[,sig]<-rt[,sig]/rt[,"CD8A"]
      sur.cut<-surv_cutpoint(rt,time = "RFS.time", event = "RFS", variables = sig)
      rt$score_group=ifelse(rt[,sig]>(sur.cut$cutpoint[,1]),1,0)##bug 
      #rt$score_group=rt[,sig]>median(rt[,sig])
      #rt$score_group=rt[,sig]>summary(rt[,sig])[5]
      Gcox1<-coxph(Surv(RFS.time,RFS)~score_group,data=rt)
      GSum<-summary(Gcox1)
      HR<-round(GSum$coefficients[,2],3)
      Pvalue<-round(GSum$coefficients[,5],3)
      CI<-paste0(round(GSum$conf.int[,3:4],3),collapse='-')
      coeff<-round(GSum$coefficients[1],3)
      se<-GSum$coefficients[1,3]
      low<-round(GSum$conf.int[,3],3)
      up<-round(GSum$conf.int[,4],3)
      cox.p<-data.frame('Characteristics'= sig,
                        'Hazard Ratio'= HR,
                        'CI95'=CI,
                        "coeff"=coeff,
                        "se"=se,
                        "low"=low,
                        "up"=up,
                        'P-value'=Pvalue,
                        "set"=set,
                        "numberofpatients"=nrow(rt))
      colorectal_cox.sig=rbind(colorectal_cox.sig,cox.p)
      message(set,"cox regression is Done!")
      ##plot the km curves
      if (mean(rt$RFS.time) > 50) {
        rt$RFS.time <- rt$RFS.time / 30
      } else {
        rt$RFS.time <- rt$RFS.time
      }
      kmfit<- survfit(Surv(RFS.time,RFS)~score_group,data=rt)
      #data.survdiff <- survdiff(Surv(RFS.time,RFS) ~ score_group,data=rt)
      #p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
      p.val<-Pvalue
      #HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
      #up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
      #low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
      HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
      CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
      print(paste0(set,sig," p= ",round(p.val,2)," ",HR))
      # pdf(paste("plot", i, ".pdf", sep = ""))
      if (!dir.exists("sc_sp_celltypes_rfs_COHORT")){
        dir.create("./sc_sp_celltypes_rfs_COHORT")
      }
      if(p.val<0.05){
        pdf(paste0("./sc_sp_celltypes_rfs_COHORT/",set,"_",gsub("/","_",sig),"_",HR,"RFS",".pdf"),width = 6.49, height = 6.58,onefile = F )
        print(ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
                         # Change legends: title & labels
                         main = "Recurrencefree-Survival",
                         legend.title = paste0(set,"-",sig),
                         legend.labs = c("Low","High"),
                         xlim = c(0, 60),
                         xlab="Time(months)",
                         ylab="Recurrence-free-Survival",
                         size = 1,
                         #fun="cumhaz",
                         #fun='event',
                         # Add p-value and tervals
                         #pval = TRUE,
                         #test.for.trend = TRUE,###group more than 2groups
                         break.time.by = 10,
                         #conf.int = TRUE,
                         #group.by=,
                         # Add risk table
                         risk.table = TRUE,
                         tables.height = 0.185,
                         tables.theme = theme_cleantable(),
                         palette = c("#1F78B4", "#E31A1C"),
                         #ggtheme = theme_bw(), # Change ggplot2 theme
                         #font.title="OS",
                         font.main =15,
                         font.x =  15,
                         font.y = 15,
                         font.tickslab =25,
                         #在左下???标出pvalue、HR???95% CI
                         #???小的p value标为p < 0.001
                         pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                    paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n")))
        invisible(dev.off())
      }
    }
  }
}
#colorectal_cox.sig<-colorectal_cox.sig[colorectal_cox.sig$P.value<0.05,]
write.csv(colorectal_cox.sig,file="/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure3_stRNA_features/Prognosis/cox.sig_RFS.csv")###

library(DT)
datatable(colorectal_cox.sig, 
          extensions = 'Buttons',
          options = list(dom = 'Bfrtip', 
                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
#### meta 结果 汇总~
#colorectal_cox.sig<-colorectal_cox.sig[c(1,9),]
library(metafor)
colorectal_cox.sig<-colorectal_cox.sig[colorectal_cox.sig$numberofpatients>100,]
NMFcox<-split(colorectal_cox.sig,colorectal_cox.sig$Characteristics)
#names(NMFcox)<-c("NMF1", "NMF2" ,"NMF3" , "NMF4")
setnames<-names(gmt.list)
NMFcoxmeta<-c()
NMFcoxlist<-list()
for (set in setnames) {NMFcox.i<-NMFcox[[set]]
res.rma_each<- metafor::rma(yi = NMFcox.i[,4], sei = NMFcox.i[,5],method="FE")
results<-round(c(res.rma_each$pval,exp(res.rma_each$b),exp(res.rma_each$ci.lb),exp(res.rma_each$ci.ub)),3)
NMFcoxlist[[set]]<-results;
#results<-cbind(set,results)
NMFcoxmeta<-rbind(NMFcoxmeta,results)}

NMFcoxmeta<-data.frame(NMFcoxmeta)
colnames(NMFcoxmeta)<-c("pvalue","HR","95CI_L","95CI_U")
rownames(NMFcoxmeta)<-setnames
head(NMFcoxmeta)
library(DT)
datatable(NMFcoxmeta,
          options = list(pageLength = 5))
datatable(NMFcoxmeta, 
          extensions = 'Buttons',
          options = list(dom = 'Bfrtip', 
                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
write.csv(NMFcoxmeta,file="/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure3_stRNA_features/Prognosis/colorectal_cox.sig_RFS_meta_results.csv")###
save(NMFcoxmeta,colorectal_cox.sig,file="/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure3_stRNA_features/Prognosis/colorectal_cox.sig_RFS_meta.Rdata")
##
######forest
library("forestplot")
library("magrittr")
library("checkmate")
library(grid)
##########
rt<-NMFcoxmeta[-(13:17),]
colnames(rt)
rt$Hazard.Ratio<-round(rt$HR,2)
rt$up<-round(rt$`95CI_U`,2)
rt$low<-round(rt$`95CI_L`,2)
#colnames(rt)
#rt<-rt[order(rt$Hazard.Ratio,decreasing = T),]
######
tabletext <- cbind(c("\nCell names",NA,rownames(rt),NA),
                   c("Hazard Ratio\n(95% CI)",NA, 
                     paste0(format(rt$HR,nsmall=2),
                            " (",format(rt$low,nsmall = 2),"-",format(rt$up,nsmall = 2),")",sep=""),NA),
                   c("p-value",NA,rt$pvalue,NA))
pdf(paste0("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure3_stRNA_features/Prognosis/","forestplot_RFS.pdf"),width = 7,height = 5,onefile=FALSE)
#pdf(paste0(cellnames[3],"forestplot_RFS.pdf"),width = 8,height = 5,onefile=FALSE)
forestplot(labeltext=tabletext, #
           mean=c(NA,NA,rt$HR,NA),#HR
           lower=c(NA,NA,rt$low,NA), #95%
           upper=c(NA,NA,rt$up,NA),#95%
           #title="Hazard Ratio",
           graph.pos=2,#�֦�
           graphwidth = unit(.3,"npc"),#
           #fn.ci_norm="fpDrawDiamondCI",#
           col=fpColors(box="steelblue", lines="black", zero = "black"),#
           #col=fpColors(box="black", lines="black", zero = "black"),#
           #boxsize=c(NA,NA,NA,rt$numberofpatients,NA)/200,#
           boxsize=0.3,
           #lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#
           zero=1,#zero
           #xlog = TRUE,
           lwd.zero=1,#zero
           #grid = structure(c(rt[1,]$Hazard.Ratio), gp = gpar(col = "black", lty=2,lwd=2)),#()
           xticks = c(0.5,0.75, 1,1.25,1.5),#
           #clip = c(0.1,2.5), 
           #lwd.xaxis=2,#X
           #xlab="     <-Favour Combination  Therapy       Favour Target  Therapy->",#X
           hrzl_lines=list("3" = gpar(lwd=2, col="black"),#
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#
                           "16" = gpar(lwd=2, col="black")),#"nrow(rt)+5
           txt_gp=fpTxtGp(label=gpar(cex=1.25),#
                          ticks=gpar(cex=1.25),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.25)),
           #is.summary = c(T,rep(F,27)),#
           #lineheight = unit(.75,"cm"),#
           align=c("l","c","c"),#
           #cex=10,
           colgap = unit(0.1,"cm"),#
           #mar=unit(rep(0.25, times = 4), "cm"),#
           new_page = T#
)
dev.off()
#########
#FOR OS
b<-c()
for (set in names(phe)){
  myphe<-phe[[set]]
  a<-ifelse("OS" %in% colnames(myphe),"OS","NO")
  a<-cbind(a,set)
  b<-rbind(b,a)
}
b
#genes<-c("GPR171","CD8A")
colorectal_cox.sig<-c()
for(set in names(exp)[c(1,4:9,11:14,18,19)]){
  gsva_es<-Scorelist[[set]]
  colnames(gsva_es)<-c('geo_accession',names(gmt.list))
  myphe<-phe[[set]]
  #myexp<-exp[[set]]
  #mygenexp<-myexp[rownames(myexp) %in% genes,]
  #mygenexp<-cbind(geo_accession=colnames(mygenexp),t(mygenexp))
  #myphe<-phe[[set]]
  #myphe<-merge(mygenexp,myphe,by="geo_accession")
  #colnames(PANCAN.survival)
  rt<-merge(gsva_es,myphe,by="geo_accession")
  #rt[,"CD8A"]<-as.numeric(rt[,"CD8A"])
  rt$OS.time<-as.numeric(as.character(rt$OS_time))
  rt$OS<-as.numeric(as.character(rt$OS))
  rt<-rt[which(rt$OS!='NA'),]
  if(nrow(rt)>10){
    for (sig in names(gmt.list)){
      rt[,sig]<-as.numeric(as.character(rt[,sig]))
      #rt[,sig]<-rt[,sig]/rt[,"CD8A"]
      sur.cut<-surv_cutpoint(rt,time = "OS.time", event = "OS", variables = sig)
      rt$score_group=ifelse(rt[,sig]>(sur.cut$cutpoint[,1]),1,0)##bug 
      #rt$score_group=rt[,sig]>median(rt[,sig])
      Gcox1<-coxph(Surv(OS.time,OS)~score_group,data=rt)
      GSum<-summary(Gcox1)
      HR<-round(GSum$coefficients[,2],3)
      Pvalue<-round(GSum$coefficients[,5],4)
      CI<-paste0(round(GSum$conf.int[,3:4],3),collapse='-')
      coeff<-round(GSum$coefficients[1],3)
      se<-GSum$coefficients[1,3]
      low<-round(GSum$conf.int[,3],3)
      up<-round(GSum$conf.int[,4],3)
      cox.p<-data.frame('Characteristics'= sig,
                        'Hazard Ratio'= HR,
                        'CI95'=CI,
                        "coeff"=coeff,
                        "se"=se,
                        "low"=low,
                        "up"=up,
                        'P-value'=Pvalue,
                        "set"=set,
                        "numberofpatients"=nrow(rt))
      colorectal_cox.sig=rbind(colorectal_cox.sig,cox.p)
      message(set,"cox regression is Done!")
      kmfit<- survfit(Surv(OS.time,OS)~score_group,data=rt)
      p.val<-Pvalue
      HR <- paste("Hazard Ratio= ", round(HR,2), sep = "")
      CI <- paste("95% CI: ", paste(round(low,2), round(up,2), sep = " - "), sep = "")
      print(paste0(set,sig," p= ",round(p.val,2)," ",HR))
      # pdf(paste("plot", i, ".pdf", sep = ""))
      if (!dir.exists("sc_sp_celltypes_os_COHORT")){
        dir.create("./sc_sp_celltypes_os_COHORT")
      }
      if(p.val<0.05){
        pdf(paste0("./sc_sp_celltypes_os_COHORT/",set,"_",gsub("/","_",sig),"_",HR,"OS",".pdf"),width = 6.49, height = 6.58,onefile = F )
        print(ggsurvplot(kmfit,#surv.median.line = "hv", # Add medians survival
                         # Change legends: title & labels
                         main = "Survival curve",
                         legend.title = paste0(set,"-",sig),
                         legend.labs = c("Low","High"),
                         #xlim = c(0, 200),
                         xlab="Time(months)",
                         ylab="Overall Survival probability",
                         size = 0.5,
                         #fun="cumhaz",
                         #fun='event',
                         # Add p-value and tervals
                         #pval = TRUE,
                         #test.for.trend = TRUE,###group more than 2groups
                         break.time.by = 30,
                         #conf.int = TRUE,
                         #group.by=,
                         # Add risk table
                         risk.table = TRUE,
                         tables.height = 0.185,
                         tables.theme = theme_cleantable(),
                         palette = c("#1F78B4", "#E31A1C"),
                         #ggtheme = theme_bw(), # Change ggplot2 theme
                         #font.title="OS",
                         font.main =15,
                         font.x =  15,
                         font.y = 15,
                         font.tickslab =25,
                         #在左下???标出pvalue、HR???95% CI
                         #???小的p value标为p < 0.001
                         pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                    paste("p = ",round(p.val,3), sep = "")), HR, CI, sep = "\n")))
        invisible(dev.off())
      }
    }
  }
}
#colorectal_cox.sig<-colorectal_cox.sig[colorectal_cox.sig$P.value<0.05,]
######################################################################
write.csv(colorectal_cox.sig,file="/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure3_stRNA_features/Prognosis/cox.sig_OS.csv")###
#### meta 结果 汇总~ 
library(metafor)
NMFcox<-split(colorectal_cox.sig,colorectal_cox.sig$Characteristics)
#names(NMFcox)<-c("C1", "C2" ,"C3","C4")
setnames<-names(gmt.list)
NMFcoxmeta<-c()
NMFcoxlist<-list()
for (set in setnames) {NMFcox.i<-NMFcox[[set]]
res.rma_each<- metafor::rma(yi = NMFcox.i[,4], sei = NMFcox.i[,5],method="FE")
results<-round(c(res.rma_each$pval,exp(res.rma_each$b),exp(res.rma_each$ci.lb),exp(res.rma_each$ci.ub)),3)
NMFcoxlist[[set]]<-results;
#results<-cbind(set,results)
NMFcoxmeta<-rbind(NMFcoxmeta,results)}

NMFcoxmeta<-data.frame(NMFcoxmeta)
colnames(NMFcoxmeta)<-c("pvalue","HR","95CI_L","95CI_U")
rownames(NMFcoxmeta)<-names(gmt.list)

head(NMFcoxmeta)
library(DT)
datatable(NMFcoxmeta, filter = "top", 
          options = list(pageLength = 5))
datatable(NMFcoxmeta, 
          extensions = 'Buttons',
          options = list(dom = 'Bfrtip', 
                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
write.csv(NMFcoxmeta,file="/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure3_stRNA_features/Prognosis/colorectal_cox.sig_OS_meta_results.csv")####
######
######forest
library("forestplot")
library("magrittr")
library("checkmate")
library(grid)
##########
rt<-NMFcoxmeta[-(13:17),]
colnames(rt)
rt$Hazard.Ratio<-round(rt$HR,2)
rt$up<-round(rt$`95CI_U`,2)
rt$low<-round(rt$`95CI_L`,2)
#colnames(rt)
#rt<-rt[order(rt$Hazard.Ratio,decreasing = T),]
######
tabletext <- cbind(c("\nCell names",NA,rownames(rt),NA),
                   c("Hazard Ratio\n(95% CI)",NA, 
                     paste0(format(rt$HR,nsmall=2),
                            " (",format(rt$low,nsmall = 2),"-",format(rt$up,nsmall = 2),")",sep=""),NA),
                   c("p-value",NA,rt$pvalue,NA))

#pdf(paste0(cellnames[1],"forestplot_OS.pdf"),width = 8,height = 5,onefile=FALSE)
pdf(paste0("/home/data/gaoyuzhen/Projects/LM_spatial/ResultsLM/Figure3_stRNA_features/Prognosis/","forestplot_OS.pdf"),width = 7,height = 5,onefile=FALSE)
forestplot(labeltext=tabletext, #
           mean=c(NA,NA,rt$HR,NA),#HR
           lower=c(NA,NA,rt$low,NA), #95%
           upper=c(NA,NA,rt$up,NA),#95%
           #title="Hazard Ratio",
           graph.pos=2,#�֦�
           graphwidth = unit(.3,"npc"),#
           #fn.ci_norm="fpDrawDiamondCI",#
           col=fpColors(box="steelblue", lines="black", zero = "black"),#
           #col=fpColors(box="black", lines="black", zero = "black"),#
           #boxsize=c(NA,NA,NA,rt$numberofpatients,NA)/200,#
           boxsize=0.3,
           #lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#
           zero=1,#zero
           #xlog = TRUE,
           lwd.zero=1,#zero
           #grid = structure(c(rt[1,]$Hazard.Ratio), gp = gpar(col = "black", lty=2,lwd=2)),#()
           xticks = c(0.5,0.75, 1,1.25,1.5),#
           #clip = c(0.1,2.5), 
           #lwd.xaxis=2,#X
           #xlab="     <-Favour Combination  Therapy       Favour Target  Therapy->",#X
           hrzl_lines=list("3" = gpar(lwd=2, col="black"),#
                           #"4" = gpar(lwd=60,lineend="butt", columns=c(1:4), col="#99999922"),#
                           "15" = gpar(lwd=2, col="black")),#"nrow(rt)+5
           txt_gp=fpTxtGp(label=gpar(cex=1.25),#
                          ticks=gpar(cex=1.25),
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.25)),
           #is.summary = c(T,rep(F,27)),#
           #lineheight = unit(.75,"cm"),#
           align=c("l","c","c"),#
           #cex=10,
           colgap = unit(0.1,"cm"),#
           #mar=unit(rep(0.25, times = 4), "cm"),#
           new_page = T#
)
dev.off()
##
