setwd(paste0("../results/", ancestries, "/", split, "/", cvsall$id[i]))
ref_leg_EUR$do_sim <- NA
ref_leg_AFR$do_sim <- NA
ref_leg_EAS$do_sim <- NA
start <- cvsall$pos[i] - 125000
stop <- cvsall$pos[i] + 125000


ref_leg_EUR$do_sim <- ref_leg_EUR$position < stop   &  
  ref_leg_EUR$position > start & 
  (ref_leg_EUR$EUR > 0.01 & ref_leg_EUR$EUR < 0.99) & 
  (ref_leg_EUR$AFR > 0.01 & ref_leg_EUR$AFR < 0.99) & 
  (ref_leg_EUR$EAS > 0.01 & ref_leg_EUR$EAS < 0.99) & 
  ref_leg_EUR$TYPE=="Biallelic_SNP"

ref_leg_AFR$do_sim <- ref_leg_AFR$position < stop   &
  ref_leg_AFR$position > start &  
  (ref_leg_AFR$AFR > 0.01 & ref_leg_AFR$AFR < 0.99) &
  (ref_leg_AFR$EUR > 0.01 & ref_leg_AFR$EUR < 0.99) &
  (ref_leg_AFR$EAS > 0.01 & ref_leg_AFR$EAS < 0.99) & 
  ref_leg_AFR$TYPE=="Biallelic_SNP"

ref_leg_EAS$do_sim <- ref_leg_EAS$position < stop   &
  ref_leg_EAS$position > start &  
  (ref_leg_EAS$EUR > 0.01 & ref_leg_EAS$EUR < 0.99) &
  (ref_leg_EAS$AFR > 0.01 & ref_leg_EAS$AFR < 0.99) & 
  (ref_leg_EAS$EAS > 0.01 & ref_leg_EAS$EAS < 0.99) & 
  ref_leg_EAS$TYPE=="Biallelic_SNP"


LD_EUR <- simGWAS:::wcor2(as.matrix( freq_EUR[,colnames(freq_EUR) %in% ref_leg_EUR[do_sim==T,id]] ), freq_EUR$Probability)
diag(LD_EUR) <- 1.000


LD_AFR <- simGWAS:::wcor2(as.matrix( freq_AFR[,colnames(freq_AFR) %in% ref_leg_AFR[do_sim==T,id]] ), freq_AFR$Probability)
diag(LD_AFR) <- 1.000

LD_EAS <- simGWAS:::wcor2(as.matrix( freq_EAS[,colnames(freq_EAS) %in% ref_leg_EAS[do_sim==T,id]] ), freq_EAS$Probability)
diag(LD_EAS) <- 1.000

write.table(LD_AFR, file = "MseCAVIAR/AFR.ld", sep = " ", row.names = F, col.names = F)
write.table(LD_EUR, file = "MseCAVIAR/EUR.ld", sep = " ", row.names = F, col.names = F)
write.table(LD_EAS, file = "MseCAVIAR/EAS.ld", sep = " ", row.names = F, col.names = F)



afrss <- fread(paste0("coloc+SuSiEx/AFR_gwas", i, ".txt"))
eurss <- fread(paste0("coloc+SuSiEx/EUR_gwas", i, ".txt"))
easss <- fread(paste0("coloc+SuSiEx/EAS_gwas", i, ".txt"))

keep <- intersect(afrss$snp, intersect(eurss$snp, easss$snp))

for (j in 1:100) {
  
  afrss <- fread(paste0("coloc+SuSiEx/AFR_gwas", j, ".txt"))
  eurss <- fread(paste0("coloc+SuSiEx/EUR_gwas", j, ".txt"))
  easss <- fread(paste0("coloc+SuSiEx/EAS_gwas", j, ".txt"))
  
  
  
  afrss2 <- subset(afrss, snp %in% keep)
  eurss2 <- subset(eurss, snp %in% keep)
  easss2 <- subset(easss, snp %in% keep)
  
  afrss2 <- afrss2[,c(2,8)]
  eurss2 <- eurss2[,c(2,8)]
  easss2 <- easss2[,c(2,8)]
  
  
  write.table(afrss2, file=paste0("MseCAVIAR/AFR_gwas", j, ".zscores"), sep = "\t", row.names = F, col.names = F)
  write.table(eurss2, file=paste0("MseCAVIAR/EUR_gwas", j, ".zscores"), sep = "\t", row.names = F, col.names = F)
  write.table(easss2, file=paste0("MseCAVIAR/EAS_gwas", j, ".zscores"), sep = "\t", row.names = F, col.names = F)
  
  
}

afrss <- fread(paste0("coloc+SuSiEx/AFR_eqtl", i, ".txt"))
eurss <- fread(paste0("coloc+SuSiEx/EUR_eqtl", i, ".txt"))
easss <- fread(paste0("coloc+SuSiEx/EAS_eqtl", i, ".txt"))


keep <- intersect(afrss$snp, intersect(eurss$snp, easss$snp))

for (j in 1:100) {
  
  afrss <- fread(paste0("coloc+SuSiEx/AFR_eqtl", j, ".txt"))
  eurss <- fread(paste0("coloc+SuSiEx/EUR_eqtl", j, ".txt"))
  easss <- fread(paste0("coloc+SuSiEx/EAS_eqtl", j, ".txt"))
  
  
  afrss2 <- subset(afrss, snp %in% keep)
  eurss2 <- subset(eurss, snp %in% keep)
  easss2 <- subset(easss, snp %in% keep)
  
  
  
  afrss2 <- afrss2[,c(2,8)]
  eurss2 <- eurss2[,c(2,8)]
  easss2 <- easss2[,c(2,8)]
  
  
  write.table(afrss2, file=paste0("MseCAVIAR/AFR_eqtl", j, ".zscores"), sep = "\t", row.names = F, col.names = F)
  write.table(eurss2, file=paste0("MseCAVIAR/EUR_eqtl", j, ".zscores"), sep = "\t", row.names = F, col.names = F)
  write.table(easss2, file=paste0("MseCAVIAR/EAS_eqtl", j, ".zscores"), sep = "\t", row.names = F, col.names = F)
  
  
}


for (j in 1:100) {
  
  file_name <- paste0("MseCAVIAR/gwas", j, ".txt")
  
  content <- c(paste0("AFR_gwas", j, ".zscores"), paste0("EUR_gwas", j, ".zscores"), paste0("EAS_gwas", j, ".zscores"))
  
  writeLines(content, file_name)
  
}

for (j in 1:100) {
  
  file_name <- paste0("MseCAVIAR/eqtl", j, ".txt")
  
  content <- c(paste0("AFR_eqtl", j, ".zscores"), paste0("EUR_eqtl", j, ".zscores"), paste0("EAS_eqtl", j, ".zscores"))
  
  writeLines(content, file_name)
  
}


writeLines(c("AFR.ld", "EUR.ld", "EAS.ld"), "MseCAVIAR/ldfiles.txt")
#}
