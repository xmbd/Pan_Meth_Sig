## Instrument variable analysis 
=====================================================================================================================================================
# Mutation type : missense mutation
library(stringr)
library(ivpack)
library(parallel)
# Load DNA methylation signatures
load("~/hyper_hypo_coef_H.RData") # hyper_hypo_coef(methylation signature scores)
# Input gene expression profile 
gene_matrix <- read.table("~/reference/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv.gz", header = T, sep = "\t") #20531 11070
gene <- gene_matrix
gene$gene_id <- str_split(gene$gene_id, '[|]', simplify = T)[,1]
gene <- gene[!duplicated(gene$gene_id),]
rownames(gene) <- gene$gene_id
gene <- gene[,-1]
gene <- as.data.frame(t(gene))
gene <- log2(gene + 1)

## Gene containg both gene expressions and mutations 
genelist <- read.csv("~/intersect_genelist.csv")

## Load missnese mutation
load("~/reference/mutation/TCGA_snv_missense.RData")    ##snv_mis   
snv_domin <- snv_mis[colnames(snv_mis) %in% genelist$x]
snv_domin$sample_ID <- str_sub(rownames(snv_domin), 1, 12)
cof_snv <- merge(hyper_hypo_coef, snv_domin)            

### Association analysis of gene mutaions with DNA methylation signatures
snv_cor_ge <- mclapply(colnames(cof_snv)[13:length(cof_snv)], function(p){
	sig <- lapply(colnames(cof_snv)[2:11], function(q){
		cat(paste(p,q), "\n")
		df <- cof_snv[c(q, "sample_ID", "cancertype", p)]
        df[,4] <- as.factor(df[,4])
		if(dim(table(df[,4])) == 1|length(which(df[,4] == 1)) < 10){
        re_p <- data.frame(sig = q,
        	    gene = p, 
        	    mut_rate = NA,
                p.values= NA) 
		}else{
            tes <- wilcox.test(df[,1] ~ df[,4])
            re_p <- data.frame(sig = q,
            	gene = p, 
        	    mut_rate = length(which(df[,4]== 1))/dim(df)[1],
                p.values= tes$p.value)}
        return(re_p)
    })		
}, mc.cores = 5)

snv_re <- lapply(1:length(snv_cor_ge), function(i){
	df <- do.call(rbind, snv_cor_ge[[i]])
	return(df)
	})

snv_mis_cor <- do.call(rbind, snv_re)
snv_mis_cor <- na.omit(snv_mis_cor)

# adjust p values
mis_cor_re <- lapply(unique(snv_mis_cor$sig), function(i){
    df <- subset(snv_mis_cor, sig == i)
    df$fdr <- p.adjust(df$p.value, method = "BH")
    return(df)
})
snv_mis_cor <- do.call(rbind, mis_cor_re)

snv_mis_cor_signi <- subset(snv_mis_cor, fdr < 0.1)    
write.csv(snv_mis_cor_signi, file = "~/snv_mis_correlation_signature_gene_result_significant0.1.csv", row.names = F)
write.csv(snv_mis_cor, file = "~/snv_mis_correlation_signature_gene_all_result.csv", row.names = F)



##################################
### Instrument variable analysis of frame shift mutations
## Instrumental variable analysis of missense mutations of those genes with significantly correlated with each methylation signature
# Input above associated genes 
snv_mis_cor_signi <- read.csv("~/snv_mis_correlation_signature_gene_result_significant0.1.csv")
snv_cor_sig <- as.data.frame(table(snv_mis_cor_signi$sig))

name <- c("Hyper.MS1", "Hyper.MS2","Hyper.MS3","Hypo.MS1", "Hypo.MS2", "Hypo.MS3", "Hypo.MS4", "Hypo.MS5", "Hypo.MS6", "Hypo.MS7")

library(parallel)
snv_iv_test <- mclapply(name, function(p) {
	sg <- subset(snv_mis_cor_signi, sig == p)
	gene_iv <- gene[colnames(gene) %in% sg$gene]
	gene_iv$sample_ID <- gsub("\\.", "-", str_sub(rownames(gene_iv), 1, 12))
	snv_iv <- snv_mis[colnames(snv_mis) %in% sg$gene]
	snv_iv$sample_ID <- str_sub(rownames(snv_iv), 1, 12)
    cb <- lapply(colnames(gene_iv)[1:length(gene_iv)-1], function(q) {
        print(paste(p,q))
        df <- hyper_hypo_coef[c("sample_ID", "cancertype", p)]
        gene_df <- gene_iv[c("sample_ID", q)]
        snv_df <- snv_iv[,c("sample_ID",q)]
        colnames(snv_df)[2] <- "snv"
        da <- merge(df, merge(gene_df, snv_df))
        if(length(which(da[,4] == 0)) > 0.2 * length(da[,4])){
            re_p = data.frame(sig = colnames(da[3]), 
                gene        = colnames(da[4]),
                coef        = NA,
                wu_Hausman  = NA,
                iv_p.value  = NA, 
                p.value     = NA,
                adjusted_r_squard = NA)}
        else{
        re2 <-  ivreg(formula = log2(da[,3]+1) ~ da[,4] + da[,2] -1 | da[,2] + da[, 5] - 1, x = TRUE, data = da)
        try(su <- summary(re2, vcov = sandwich, df = Inf, diagnostics = TRUE), TRUE)
        try(re_p <- data.frame(sig = colnames(da[3]), 
                  gene       = colnames(da[4]),
                  coef       = su$coefficients[1,1], 
                  wu_Hausman = su$diagnostics[2,4],
                  iv_p.value = su$diagnostics[1,4], 
                  p.value    = su$coefficients[1,4],
                  adjusted_r_squard = su$adj.r.squared), TRUE)}
        try(return(re_p), TRUE)
      })
  }, mc.cores = 5)

snv_iv_re <- lapply(1:10, function(i){
	df <- do.call(rbind, snv_iv_test[[i]])
	return(df)
	})
snv_mis_iv <- do.call(rbind, snv_iv_re)
snv_mis_iv  <- na.omit(snv_mis_iv)

snv_re <- lapply(unique(snv_mis_iv$sig), function(i){
    df <- subset(snv_mis_iv, sig == i)
    df$fdr <- p.adjust(df$p.value, method = "BH")
    df$wu.fdr <- p.adjust(df$wu_Hausman, method = "BH")
    df$weak.fdr <- p.adjust(df$iv_p.value, method = "BH")   
    return(df)
})
snv_mis_iv_new <- do.call(rbind, snv_re)
snv_mis_iv_signifi <- subset(snv_mis_iv, wu_Hausman < 0.05 &iv_p.value < 0.05&fdr < 0.05)
write.csv(snv_mis_iv_signifi_fdr, file = "~/snv_mis_iv_fdr0.05.csv", row.names = F)
write.csv(snv_mis_iv, file = "~/snv_iv_mis_all_reslut.csv", row.names = F)








# Mutation type : Frame Shift Mutation 
# Load DNA methylation signatures
load("~/hyper_hypo_coef_H.RData")

# Input gene expression profile 
gene_matrix <- read.table("~/reference/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv.gz", header = T, sep = "\t") #20531 11070
gene <- gene_matrix
gene$gene_id <- str_split(gene$gene_id, '[|]', simplify = T)[,1]
gene <- gene[!duplicated(gene$gene_id),]
rownames(gene) <- gene$gene_id
gene <- gene[,-1]
gene <- as.data.frame(t(gene))
gene <- log2(gene + 1)

## Gene containg both gene expressions and mutations 
genelist <- read.csv("~/intersect_genelist.csv")

## Load frame shift mutation
load("/home/qingqing/reference/mutation/TCGA_frame_shift.RData") # snv_fram
snv_domin <- snv_fram[colnames(snv_fram) %in% genelist$x]
snv_domin$sample_ID <- str_sub(rownames(snv_domin), 1, 12)
cof_snv <- merge(hyper_hypo_coef, snv_domin)

## Association analysis of gene mutaions with DNA methylation signatures
snv_cor_ge <- mclapply(colnames(cof_snv)[13:length(cof_snv)], function(p){
    sig <- lapply(colnames(cof_snv)[2:11], function(q){
        cat(paste(p,q), "\n")
        df <- cof_snv[c(q, "sample_ID", "cancertype", p)]
        df[,4] <- as.factor(df[,4])
        if(dim(table(df[,4])) == 1|length(which(df[,4] == 1)) < 10){
        re_p <- data.frame(sig = q,
                gene = p, 
                mut_rate = NA,
                p.values= NA) 
        }else{
            tes <- wilcox.test(df[,1] ~ df[,4])
            re_p <- data.frame(sig = q,
                gene = p, 
                mut_rate = length(which(df[,4]== 1))/dim(df)[1],
                p.values= tes$p.value)}
        return(re_p)
    })      
}, mc.cores = 5)

snv_re <- lapply(1:length(snv_cor_ge), function(i){
    df <- do.call(rbind, snv_cor_ge[[i]])
    return(df)
    })

snv_fram_cor <- do.call(rbind, snv_re)
snv_fram_cor <- na.omit(snv_fram_cor)
# Adjust p value
frame_cor_re <- lapply(unique(snv_fram_cor$sig), function(i){
    df <- subset(snv_fram_cor, sig == i)
    df$fdr <- p.adjust(df$p.value, method = "BH")
    return(df)
})
snv_fram_cor <- do.call(rbind, frame_cor_re)
snv_fram_cor_signi <- subset(snv_fram_cor, fdr < 0.1)
write.csv(snv_fram_cor_signi, file = "~/snv_frame_correlation_signature_gene_result_significant_fdr0.1.csv", row.names = F)
write.csv(snv_fram_cor, file = "~/snv_frame_correlation_signature_gene_all_result.csv", row.names = F)


##################################
## Instrument variable analysis of frame shift mutations
### Instrumental variable analysis of missense mutations of those genes with significantly correlated with each methylation signature
# Input above associated genes 
snv_fram_cor_signi <- read.csv("~/snv_frame_correlation_signature_gene_result_significant_fdr0.1.csv")

name <- c("Hyper.MS1", "Hyper.MS2","Hyper.MS3","Hypo.MS1", "Hypo.MS2", "Hypo.MS3", "Hypo.MS4", "Hypo.MS5", "Hypo.MS6", "Hypo.MS7")

snv_iv <- mclapply(name, function(p) {
    sg <- subset(snv_fram_cor_signi, sig == p)
    gene_iv <- gene[colnames(gene) %in% sg$gene]
    gene_iv$sample_ID <- gsub("\\.", "-", str_sub(rownames(gene_iv), 1, 12))
    snv_iv <- snv_fram[colnames(snv_fram) %in% sg$gene]
    snv_iv$sample_ID <- str_sub(rownames(snv_iv), 1, 12)
    cb <- lapply(colnames(gene_iv)[1:length(gene_iv)-1], function(q) {
        print(paste(p,q))
        df <- hyper_hypo_coef[c("sample_ID", "cancertype", p)]
        gene_df <- gene_iv[c("sample_ID", q)]
        snv_df <- snv_iv[,c("sample_ID",q)]
        colnames(snv_df)[2] <- "snv"
        da <- merge(df, merge(gene_df, snv_df))
        if(dim(table(da[,4])) == 1 | dim(table(da[,5])) == 1|length(which(da[,4] == 0)) > 0.2 * length(da[,4])){
            re_p = data.frame(sig = colnames(da[3]), 
                gene        = colnames(da[4]),
                coef        = NA,
                wu_Hausman  = NA,
                iv_p.value  = NA, 
                p.value     = NA,
                adjusted_r_squard = NA)}
        else{
        re2 <-  ivreg(formula = log2(da[,3]+1) ~ da[,4] + da[,2] -1 | da[,2] + da[, 5] - 1, x = TRUE, data = da)
        try(su <- summary(re2, vcov = sandwich, df = Inf, diagnostics = TRUE), TRUE)
        try(re_p <- data.frame(sig = colnames(da[3]), 
                  gene       = colnames(da[4]),
                  coef       = su$coefficients[1,1], 
                  wu_Hausman = su$diagnostics[2,4],
                  iv_p.value = su$diagnostics[1,4], 
                  p.value    = su$coefficients[1,4],
                  adjusted_r_squard = su$adj.r.squared), TRUE)}
        try(return(re_p), TRUE)
      })
  }, mc.cores = 5)

snv_iv_re <- lapply(1:10, function(i){
    df <- do.call(rbind, snv_iv[[i]])
    return(df)
    })
snv_fram_iv <- do.call(rbind, snv_iv_re)
snv_fram_iv  <- na.omit(snv_fram_iv)
snv_fram_iv$p.values <- as.numeric(snv_fram_iv$p.value)

frame_re <- lapply(unique(snv_fram_iv$sig), function(i){
    df <- subset(snv_fram_iv, sig == i)
    df$fdr <- p.adjust(df$p.value, method = "BH")
    df$wu.fdr <- p.adjust(df$wu_Hausman, method = "BH")
    df$weak.fdr <- p.adjust(df$iv_p.value, method = "BH")   
    return(df)
})
snv_fram_iv_new <- do.call(rbind, frame_re)

snv_fram_iv_signifi_fdr  <- subset(snv_fram_iv_new, wu.fdr < 0.05 &weak.fdr < 0.05&fdr < 0.05)
write.csv(snv_fram_iv_signifi_fdr, file = "~/snv_fram_iv_fdr0.05.csv", row.names = F)
write.csv(snv_fram_iv, file = "~/snv_iv_fram_all_reslut.csv", row.names = F)
