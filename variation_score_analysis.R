##Calculation of the inter-tumour variation score of methylation signatures
load("~/hyper_hypo_coef_H.RData")   
var_can_var <- NULL
for (i in colnames(hyper_hypo_coef)[1:10]){
          aov_sig <- summary(aov(lm(hyper_hypo_coef[,i]~hyper_hypo_coef[,12])))
          var_re <- data.frame(sig = i, var =  aov_sig[[1]][1,2]/sum( aov_sig[[1]][1,2]+ aov_sig[[1]][2,2]))
          var_can_var <- rbind(var_can_var, var_re)
}
write.csv(var_can_var, file = "~/var_sig/variation_between_cancer_with_sig.csv", row.names =F)
