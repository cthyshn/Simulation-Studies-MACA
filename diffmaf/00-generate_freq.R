# generates the freq file for different ancestries needed by simGWAS 
library(data.table)
chr <- 2

ref_hap <- fread("../data/chr2_mac_gt0.005.hap.gz") 
ref_leg <- fread("../data/chr2_maf_gt0.005_chr2.leg.gz") 
ref_sam <- fread("../data/1000GP_Phase3.sample2")
EUR_ind <- which(ref_sam$GROUP == "EUR") 
AFR_ind <- which(ref_sam$GROUP == "AFR")
EAS_ind <- which(ref_sam$GROUP == "EAS")
SAS_ind <- which(ref_sam$GROUP == "SAS")
AMR_ind <- which(ref_sam$GROUP == "AMR")


# Process data to simGWAS format
# FROM ref_hap:
#                 samples
#       s1hap1 s1hap2  s2hap1 s2hap2  etc.
# snp1     0     1       1       0
# snp2     1     0       1       1    ...
# etc.             ...
#
# TO freq:
#           SNPS
#        snp1 snp2 etc.
# s1hap1    0    1
# s1hap2    1    0
# s2hap1    1    1
# s2hap2    0    1
h1 <- ref_hap[, .SD, .SDcols=seq(1,ncol(ref_hap)-1,by=2)]  # extract haplotype data into 2 matrices
h2 <- ref_hap[, .SD, .SDcols=seq(2,ncol(ref_hap)  ,by=2)]  # split up by individuals (every other column)


h1_EUR <- h1[,..EUR_ind]
h2_EUR <- h2[,..EUR_ind]

freq_EUR <- as.data.frame(t(cbind(h1_EUR,h2_EUR))+1)
colnames(freq_EUR) <- ref_leg$id

freq_EUR$Probability <- 1/nrow(freq_EUR)

maf_eur <- (apply(freq_EUR, 2, mean)) - 1


save(freq_EUR, file = "../data/simGWAS_freq_EUR.RData")


h1_AFR <- h1[,..AFR_ind]
h2_AFR <- h2[,..AFR_ind]
freq_AFR <- as.data.frame(t(cbind(h1_AFR, h2_AFR))+1)
colnames(freq_AFR) <- ref_leg$id
freq_AFR$Probability <- 1/nrow(freq_AFR)

maf_afr <- (apply(freq_AFR, 2, mean)) - 1
length(which(ref_leg$AFR-maf_afr < 1e-3))

save(freq_AFR, file = "../data/simGWAS_freq_AFR.RData")


h1_EAS <- h1[,..EAS_ind]
h2_EAS <- h2[,..EAS_ind]
freq_EAS <- as.data.frame(t(cbind(h1_EAS, h2_EAS))+1)
colnames(freq_EAS) <- ref_leg$id
freq_EAS$Probability <- 1/nrow(freq_EAS)

maf_eas <- (apply(freq_EAS, 2, mean)) - 1
length(which(ref_leg$EAS-maf_eas < 1e-3))

save(freq_EAS, file = "../data/simGWAS_freq_EAS.RData")


h1_SAS <- h1[,..SAS_ind]
h2_SAS <- h2[,..SAS_ind]
freq_SAS <- as.data.frame(t(cbind(h1_SAS, h2_SAS))+1)
colnames(freq_SAS) <- ref_leg$id
freq_SAS$Probability <- 1/nrow(freq_SAS)

save(freq_SAS, file = "../data/simGWAS_freq_SAS.RData")


h1_AMR <- h1[,..AMR_ind]
h2_AMR <- h2[,..AMR_ind]
freq_AMR <- as.data.frame(t(cbind(h1_AMR, h2_AMR))+1)
colnames(freq_AMR) <- ref_leg$id
freq_AMR$Probability <- 1/nrow(freq_AMR)

save(freq_AMR, file = "../data/simGWAS_freq_AMR.RData")




