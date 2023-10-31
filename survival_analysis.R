############### overall survival analysis
library("survival")
library("survminer")
library(forestplot)
library(dplyr)
library(stringr)
library(ggpubr)
library(reshape2)
library(ggplot2)
cli <- read.csv("~/TCGA-Clinical-Data-Resource.csv") ###clinical information
colnames(cli)[1] <- "sample_ID" 

os <- cli[,c(1,3,4,6,25,26)]
os <- os[which(os$OS == 1 |os$OS == 0),]
os$OS <- as.numeric(as.character(os$OS))
os$OS.time <- as.numeric(as.character(os$OS.time))

sp$predict <- gsub("\\.", "-", sp$predict)
colnames(hyper_hypo_coef) <- gsub("\\.", "-", colnames(hyper_hypo_coef))
names <- c("Hyper-MS1", "Hyper-MS2","Hyper-MS3","Hypo-MS1", "Hypo-MS2", "Hypo-MS3", "Hypo-MS4", "Hypo-MS5", "Hypo-MS6", "Hypo-MS7")

###data prepared
sp_hyper <- subset(sp, predict %in% c("Hyper.MS1", "Hyper.MS2")) %>% dplyr::filter(!duplicated(sample_ID)) %>% dplyr::rename(hyper_predict  = "predict")
sp_hypo <- subset(sp, !predict %in% c("Hyper.MS1", "Hyper-MS2", "Hyper.MS3")) %>% dplyr::filter(!duplicated(sample_ID)) %>% dplyr::rename(hypo_predict  = "predict")

cof_os <- merge(hyper_hypo_coef, os) %>% filter(!duplicated(sample_ID))
cof_os$OS.Time <- round(cof_os$OS.time/30, 0)
cof_os_predict <- merge(cof_os, merge(sp_hyper[c(1,3)], sp_hypo[c(1,3)]))  ## 2386   18

###calculation of the hazard ratios for each signature
one.cancer <- c("Hyper.MS3", "Hypo.MS1", "Hypo.MS2", "Hypo.MS3", "Hypo.MS5", "Hypo.MS7")
cox_sig <- NULL
for (i in colnames(cof_os_predict)[2:11]){
      ms <- cof_os_predict[c("sample_ID","cancertype", i, "OS", "OS.Time")]
      ms_can <- subset(ms, sample_ID %in% subset(sp, predict == i)$sample_ID)
      can.num <- melt(table(ms_can$cancertype))
      if(i %in% c(one.cancer)){
      cancer.pos <- subset(can.num, value>50)$Var1
      }else{
      cancer.pos <- subset(can.num, value>0)$Var1
      }
      for (j in cancer.pos){
          cat(paste(i, j, "\n")) 
          ms_new <- subset(ms, cancertype == j)              
          ms_new$sig_type <- ifelse(ms_new[,i] < mean(ms_new[,i]), "Low", "High") 
          ms_new$sig_type <- factor(ms_new$sig_type, levels = c("Low", "High"))
          cox <- coxph(formula = Surv(OS.Time, OS) ~ sig_type, data = ms_new)
          pv <- round(summary(cox)$coefficients[1, 5], 4)
          hr <- round(summary(cox)$coefficients[1, 2], 2)
          re <- data.frame(sig = i,
                           cancertype = j,
                           HR = hr,
                           pvalue = pv)
          cox_sig <- rbind(cox_sig, re)
          }
        }


###Plots of significant hazard ratios of the signatures
cox_sig_signi <- subset(cox_sig, pvalue <= 0.05)

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

sig.signi.can_sur <-  lapply(unique(cox_sig_signi$sig), function(i){
                      can.ty <- unique(subset(cox_sig_signi, sig == i)$cancertype)
                      can <- lapply(can.ty, function(j) {
                              cat(paste(i, j, "\n"))               
                              ms <- subset(cof_os_predict[c("sample_ID","cancertype", i, "OS", "OS.Time")], cancertype == j)
                              ms$sig_type <- ifelse(ms[,i]< mean(ms[,i]), "Low", "High") 
                              ms$sig_type <- factor(ms$sig_type, levels = c("Low", "High"))
                              lsur <- survfit(Surv(OS.Time, OS) ~ sig_type, data = ms)
                              cox <- coxph(formula = Surv(OS.Time, OS) ~ sig_type, data = ms)
                              pv <- round(summary(cox)$coefficients[1, 5], 4)
                              hr <- round(summary(cox)$coefficients[1, 2], 2)
                              if(pv < 0.0001){
                                pv1 <- "p < 0.0001"
                              }else{pv1 <- paste0("p = ", pv)}
                              p <- ggsurvplot(lsur, 
                                    title = j, 
                                    xlab = "Time (months)", 
                                    ylab = "\n Overall survival",
                                    #surv.median.line = "hv", 
                                    font.title = c(15),
                                    font.legend = c(15),
                                    #legend = c(0.6, 0.8),
                                    legend.title = gsub("\\.","-", i), size = 1, 
                                    font.x = c(15,"black"),
                                    font.y = c(15,"black"),
                                    font.tickslab = c(15,"black"),
                                    #linetype = "sig", 
                                    #conf.int = TRUE,
                                    #add.all = TRUE,
                                    #pval = TRUE, 
                                    #pval.method = FALSE,
                                    #pval.coord = c(0.2*(max(ms$OS.Time)), 0.6),
                                    pval.size = 8, 
                                    legend.labs = c("Low", "High"),
                                    palette = c("#899DA4","#DC863B"),
                                    ggtheme = custom_theme()) +
                                    guides(colour = guide_legend(nrow = 2))
                              p$plot <- p$plot +theme(legend.background=element_blank(),
                                              plot.margin = unit(c(0.8, 0.5, 0.3, 0.2),"cm"),
                                              axis.text.x = element_text(size = 15)) +
                                        ggplot2::annotate("text",x = max(ms$OS.Time, na.rm =TRUE)* 0.45, y = 0.3,
                                              label = paste("HR :",hr)) + 
                                        ggplot2::annotate("text",x = max(ms$OS.Time, na.rm =TRUE)* 0.45, y = 0.15,
                                              label = paste("(","95%CI:",round(summary(cox)$conf.int[1,3],2),"-",round(summary(cox)$conf.int[1,4],2),")",sep = ""))+
                                        ggplot2::annotate("text",x = max(ms$OS.Time, na.rm =TRUE)* 0.45, y = 0.05,
                                                        label = pv1)
                        return(p)
                        })    
}) 
sig_res1 <- arrange_ggsurvplots(sig.signi.can_sur[[1]], nrow =1, ncol = 3)
sig_res2 <- arrange_ggsurvplots(sig.signi.can_sur[[2]], nrow=1, ncol = 2)
sig_res3 <- arrange_ggsurvplots(sig.signi.can_sur[[3]], nrow =1, ncol = 3)
sig_res4 <- arrange_ggsurvplots(sig.signi.can_sur[[4]], nrow =1, ncol = 1)
p1 <- ggarrange(plotlist=c(sig_res1, sig_res2), nrow =1, ncol =2, widths = c(3,2)) 
p2 <- ggarrange(plotlist=c(sig_res3, sig_res4), nrow =1, ncol =3, widths = c(3,1,1) )

sig_res <- arrange_ggsurvplots(
c(sig.signi.can_sur[[1]],
sig.signi.can_sur[[2]], 
sig.signi.can_sur[[3]],
sig.signi.can_sur[[4]]), ncol = 5) 
p1 <- ggarrange(plotlist=sig_res, nrow =2, ncol =1) 

pdf("~/signature.pos.cancer_survival_analysis_cancertype.median_single.cox_significant.new.pdf", height = 7, width = 16)
ggarrange(p1, nrow = 1, ncol = 1)
dev.off()



###Plots of non-significant hazard ratios of the signatures
cox_sig_not.signi <- subset(cox_sig, pvalue > 0.05)

sig.not.signi.can_sur <-  lapply(unique(cox_sig_not.signi$sig), function(i){
                      can.ty <- unique(subset(cox_sig_not.signi, sig == i)$cancertype)
                      can <- lapply(can.ty, function(j) {
                              cat(paste(i, j, "\n"))               
                              ms <- subset(cof_os_predict[c("sample_ID","cancertype", i, "OS", "OS.Time")], cancertype == j)
                              ms$sig_type <- ifelse(ms[,i]< mean(ms[,i]), "Low", "High") 
                              ms$sig_type <- factor(ms$sig_type, levels = c("Low", "High"))
                              lsur <- survfit(Surv(OS.Time, OS) ~ sig_type, data = ms)
                              cox <- coxph(formula = Surv(OS.Time, OS) ~ sig_type, data = ms)
                              pv <- round(summary(cox)$coefficients[1, 5], 4)
                              hr <- round(summary(cox)$coefficients[1, 2], 2)
                              if(pv < 0.0001){
                                pv1 <- "p < 0.0001"
                              }else{pv1 <- paste0("p = ", pv)}
                              p <- ggsurvplot(lsur, 
                                    title = j, 
                                    xlab = "Time (months)", 
                                    ylab = "\n Overall survival",
                                    #surv.median.line = "hv", 
                                    font.title = c(15),
                                    font.legend = c(15),
                                    #legend = c(0.6, 0.8),
                                    legend.title = gsub("\\.","-", i), size = 1, 
                                    font.x = c(15,"black"),
                                    font.y = c(15,"black"),
                                    font.tickslab = c(15,"black"),
                                    #linetype = "sig", 
                                    #conf.int = TRUE,
                                    #add.all = TRUE,
                                    #pval = TRUE, 
                                    #pval.method = FALSE,
                                    #pval.coord = c(0.2*(max(ms$OS.Time)), 0.6),
                                    pval.size = 6, 
                                    legend.labs = c("Low", "High"),
                                    palette = c("#899DA4","#DC863B"),
                                    ggtheme = custom_theme()) +
                                    guides(colour = guide_legend(nrow = 2))
                              p$plot <- p$plot +theme(legend.background=element_blank(),
                                              plot.margin = unit(c(0.8, 1, 0.3, 0.2),"cm"),
                                              axis.text.x = element_text(size = 15)) +
                                        ggplot2::annotate("text",x = max(ms$OS.Time, na.rm =TRUE)* 0.45, y = 0.3,
                                              label = paste("HR :",hr)) + 
                                        ggplot2::annotate("text",x = max(ms$OS.Time, na.rm =TRUE)* 0.45, y = 0.15,
                                              label = paste("(","95%CI:",round(summary(cox)$conf.int[1,3],2),"-",round(summary(cox)$conf.int[1,4],2),")",sep = ""))+
                                        ggplot2::annotate("text",x = max(ms$OS.Time, na.rm =TRUE)* 0.45, y = 0.05,
                                                        label = pv1)
                        return(p)
                        })    
}) 


res1 <- arrange_ggsurvplots(
c(sig.not.signi.can_sur[[1]],
sig.not.signi.can_sur[[2]], 
sig.not.signi.can_sur[[3]],
sig.not.signi.can_sur[[4]], 
sig.not.signi.can_sur[[5]], 
sig.not.signi.can_sur[[6]],
sig.not.signi.can_sur[[7]],
sig.not.signi.can_sur[[8]],
sig.not.signi.can_sur[[9]],
sig.not.signi.can_sur[[10]]), ncol = 5)

p1 <- ggarrange(plotlist =res1, ncol =1, nrow = 5)
pdf("~/signature.pos.cancer_survival_analysis_cancertype.median_single.cox_not.significant_1017.pdf", height = 18, width = 18, onefile = FALSE)
ggarrange(p1,ncol =1, nrow =1)
dev.off()




