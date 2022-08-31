## performing NMF to derive methylation signature
library(NMF)

### Hyper-DMPs
load("hyper_tumor_2969*2477.RData")  #hyper_tumor  2969 2477
hyper_dmp <- as.data.frame(t(hyper_tumor))
for (i in 1:dim(hyper_dmp)[2]){
	hyper_dmp[is.na(hyper_dmp[,i]),i] <- mean(hyper_dmp[,i], na.rm = TRUE)
}
hyper_dmp <- as.matrix(t(hyper_dmp))
hyper_dmp_nmf <- nmf(hyper_dmp,2:8,'ns',seed = 123456,nrun=100)

#### saving hypermethylation signature 
hyper_coef_H <- t(coef(hyper_dmp_nmf$fit[[6]]))
hyper_bas <- basis(hyper_dmp_nmf$fit[[6]])
save(hyper_coef_H, file = "hyper_coef_H.RData")
save(hyper_bas, file = "hypper_basis_W.RData")

#### saving the information of hypermethylation signature classfication
hyper_dmp_subty <- predict(hyper_dmp_nmf,what= "features",prob = TRUE)
hyper_dmp_subty <- do.call(cbind,hyper_dmp_subty)
hyper_dmp_subty <- as.data.frame(hyper_dmp_subty)
hyper_sample_subty <- predict(hyper_dmp_nmf,what= "samples",prob = TRUE)
hyper_sample_subty <- do.call(cbind,hyper_sample_subty)
hyper_sample_subty <- as.data.frame(hyper_sample_subty)
write.table(hyper_sig_subty, file = "hyper_signature_dmp_subtype.txt", sep = "\t", quote =F, row.names = TRUE)
write.table(hyper_sample_subty, file = "hyper_signature_samples_subtype.txt", sep = "\t", quote =F, row.names = TRUE)



### Hypo-DMPs
load("hypo_tumor_4349*2477.RData")  #hypo_tumor 
library(NMF)
hypo_dmp <- as.data.frame(t(hypo_tumor))
for (i in 1:dim(hypo_dmp)[2]){
	hypo_dmp[is.na(hypo_dmp[,i]),i] <- mean(hypo_dmp[,i], na.rm = TRUE)
}
hypo_dmp <- as.matrix(t(hypo_dmp))
hypo_dmp_nmf <- nmf(hypo_dmp,2:18,'ns',seed = 123456,nrun=100)

## saving hypomethylation signature 
hypo_coef_H <- t(coef(hypo_dmp_nmf$fit[[6]]))
hypo_bas <- basis(hypo_dmp_nmf$fit[[6]])
save(hypo_coef_H, file = "hypo_coef_H.RData")
save(hypo_bas, file = "hypo_basis_W.RData")

#### saving the information of hypomethylation signature classfication
hypo_dmp_subty <- predict(hypo_dmp_nmf,what= "features",prob = TRUE)
hypo_dmp_subty <- do.call(cbind,hypo_dmp_subty)
hypo_dmp_subty <- as.data.frame(hypo_dmp_subty)
hypo_sample_subty <- predict(hypo_dmp_nmf,what= "samples",prob = TRUE)
hypo_sample_subty <- do.call(cbind,hypo_sample_subty)
hypo_sample_subty <- as.data.frame(hypo_sample_subty)
write.table(hypo_sig_subty, file = "hypo_signature_dmp_subtype.txt", sep = "\t", quote =F, row.names = TRUE)
write.table(hypo_sample_subty, file = "hypo_signature_samples_subtype.txt", sep = "\t", quote =F, row.names = TRUE)















