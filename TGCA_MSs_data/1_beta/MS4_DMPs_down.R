########
######### hypo_tumor: hypo-DMPs tumor matrix
##### hypo_bas:hypomethylation E matrix
#####hyper_hypo_coef: Hyper- and Hypo-MSs activity
library(MASS)
library(NMF)
load("/data/scRNA/beta_minus/hypo_down_nmf/hyper_hypo_basis.RData")  #####TCGA-DMPs E matrix
load("/data/scRNA/beta_minus/hypo_down_nmf/hyper_hypo_coef_H.RData")   #####TCGA-DMPs MSs activities
load("/data/scRNA/beta_minus/hypo_down_nmf/hypo_tumor_4349*2477.RData")   ######TCGA input matrix

V <- as.matrix(hypo_tumor)
#fill NA's with mean beta values
V <- apply(V, 2, function(v) {v[is.na(v)] <- mean(v, na.rm=T); v})

E <- as.matrix(hypo_bas[, -8])
H <- t(as.matrix(hyper_hypo_coef[, 4:10]))
#verify nmf
summary(as.vector(V - E %*% H)) 

H_inv <- ginv(H)
E_inv <- ginv(E)
#verify VH+ == E
summary(as.vector(V%*%H_inv-E))
U <- (1.2-V) %*% H_inv
nfac_H <- nmf(U, 7, 'ns', seed = 123456, nrun=100)
Y <- coef(nfac_H)
W <- basis(nfac_H)
summary(as.vector(U - W %*% Y)) 

tcga_hypo_E_down <- W   
tcga_hypo_H_down <- Y %*% H
