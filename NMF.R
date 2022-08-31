## performing NMF to derive methylation signature
library(NMF)

load("/home/qingqing/TCGA_meth_pancancer/pancancer_hyper_hypo_nmf/hyper/hyper_tumor_2969*2477.RData")  #hyper_tumor  2969 2477
hyper_dmp <- as.data.frame(t(hyper_dmp))
for (i in 1:dim(hyper_dmp)[2]){
	hyper_dmp[is.na(hyper_dmp[,i]),i] <- mean(hyper_dmp[,i], na.rm = TRUE)
}
hyper_dmp <- as.matrix(t(hyper_dmp))
hyper_dmp_nmf <- nmf(hyper_dmp,2:8,'ns',seed = 123456,nrun=100)























