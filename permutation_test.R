#### Permutation test for paired tumor samples 

library(predictmeans)
library(foreach)

## Input arguments
load("~/pancancer_paired_450k_filter.RData")  ###pancancer_filter

status <- str_sub(rownames(pancancer_filter),14,15)
status <- ifelse(status == "01", 1, 0)
barcode <- strsplit(rownames(pancancer_filter),".",fixed = TRUE)
TSS <- sapply(barcode,function(i) b <- i[2])
part <- sapply(barcode,function(i) b <- i[3])
subject <- paste(TSS,".",part,sep="")
cancertype <- sapply(barcode,function(i) b <- i[length(i)])

##Setting multi-threaded parameters
numCores <- detectCores() - 4
registerDoParallel(numCores)

####permutation test  
perm_result <- foreach(i = colnames(pancancer_filter)) %dopar% {
					cat(i,"\n")
					betaval <- pancancer_filter[,i]
					m1 <- lmer(betaval ~ status * cancertype + (1|TSS), control=lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-8)))
					perm <- permanova.lmer(m1, nperm = 1E3)
					re <- data.frame(cg = i,status_p = perm$Perm.p[1], cancertype_p = perm$Perm.p[2],
					status_cancertype_p = perm$Perm.p[3])
					write.table(re,file="perm_probe_result.txt",row.names=F,sep = "\t",append = T, col.names =F,quote = F)
					stopImplicitCluster(numCores) 
}
























