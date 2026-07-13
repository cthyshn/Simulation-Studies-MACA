setwd(paste0("../results/", ancestries, "/", split, "/", cvsall$id[i], "/coloc+SuSiEx"))
ref_leg_EUR <- ref_leg
ref_leg_AFR <- ref_leg
ref_leg_EAS <- ref_leg
ref_leg_EUR$do_sim <- NA
ref_leg_AFR$do_sim <- NA
ref_leg_EAS$do_sim <- NA
start <- cvsall$pos[i] - 125000
stop <- cvsall$pos[i] + 125000

ref_leg_EUR$do_sim <- ref_leg_EUR$position < stop   &
  ref_leg_EUR$position > start & (ref_leg_EUR$EUR > 0.01 & ref_leg_EUR$EUR < 0.99)  & ref_leg_EUR$TYPE=="Biallelic_SNP"

ref_leg_AFR$do_sim <- ref_leg_AFR$position < stop   &
  ref_leg_AFR$position > start &  (ref_leg_AFR$AFR > 0.01 & ref_leg_AFR$AFR < 0.99) & ref_leg_AFR$TYPE=="Biallelic_SNP"

ref_leg_EAS$do_sim <- ref_leg_EAS$position < stop   &
  ref_leg_EAS$position > start &  (ref_leg_EAS$EAS > 0.01 & ref_leg_EAS$EAS < 0.99) & ref_leg_EAS$TYPE=="Biallelic_SNP"


FP_EUR <- make_GenoProbList(snps=ref_leg_EUR[do_sim==T,id], W=cvsall$id[i], freq=freq_EUR)
LD_EUR <- simGWAS:::wcor2(as.matrix( freq_EUR[,colnames(freq_EUR) %in% ref_leg_EUR[do_sim==T,id]] ), freq_EUR$Probability)
diag(LD_EUR) <- 1.0001

FP_AFR <- make_GenoProbList(snps=ref_leg_AFR[do_sim==T,id], W=cvsall$id[i], freq=freq_AFR)
LD_AFR <- simGWAS:::wcor2(as.matrix( freq_AFR[,colnames(freq_AFR) %in% ref_leg_AFR[do_sim==T,id]] ), freq_AFR$Probability)
diag(LD_AFR) <- 1.0001# LD cannot contain NAs 

FP_EAS <- make_GenoProbList(snps=ref_leg_EAS[do_sim==T,id], W=cvsall$id[i], freq=freq_EAS)
LD_EAS <- simGWAS:::wcor2(as.matrix( freq_EAS[,colnames(freq_EAS) %in% ref_leg_EAS[do_sim==T,id]] ), freq_EAS$Probability)
diag(LD_EAS) <- 1.0001# LD cannot contain NAs 

n_rep <- 100

set.seed(1234)
z_sims_EUR <- simulated_z_score(
  N0=eurcontrols, N1=eurcases, # N controls, cases
  snps=ref_leg_EUR[do_sim==T,id], # column names in freq of SNPs for which Z scores should be generated
  W=cvsall$id[i], # CVs, subset of snps
  gamma.W=log(or), # CVs' log(OR)s
  freq=freq_EUR[,names(freq_EUR) %in% c(ref_leg_EUR[do_sim==T,id],"Probability")], # Ref haplotypes
  GenoProbList = FP_EUR,
  #  LD=LD_EUR,
  #  rmvnorm_method="chol",
  nrep=n_rep
)

z_sims_AFR <- simulated_z_score(
  N0=afrcontrols, N1=afrcases, # N controls, cases
  snps=ref_leg_AFR[do_sim==T,id], # column names in freq of SNPs for which Z scores should be generated
  W=cvsall$id[i], # CVs, subset of snps
  gamma.W=log(or), # CVs' log(OR)s
  freq=freq_AFR[,names(freq_AFR) %in% c(ref_leg_AFR[do_sim==T,id],"Probability")], # Ref haplotypes
  GenoProbList = FP_AFR,
  #  LD=LD_AFR,
  #  rmvnorm_method="chol",
  nrep=n_rep
)

z_sims_EAS <- simulated_z_score(
  N0=eascontrols, N1=eascases, # N controls, cases
  snps=ref_leg_EAS[do_sim==T,id], # column names in freq of SNPs for which Z scores should be generated
  W=cvsall$id[i], # CVs, subset of snps
  gamma.W=log(or), # CVs' log(OR)s
  freq=freq_EAS[,names(freq_EAS) %in% c(ref_leg_EAS[do_sim==T,id],"Probability")], # Ref haplotypes
  GenoProbList = FP_EAS,
  #  LD=LD_EAS,
  #  rmvnorm_method="chol",
  nrep=n_rep
)

set.seed(1234)
se_sims_EUR <- sqrt(simulated_vbeta(
  N0=eurcontrols, N1=eurcases, # N controls, cases
  snps=ref_leg_EUR[do_sim==T,id], # column names in freq of SNPs for which Z scores should be generated
  W=cvsall$id[i], # CVs, subset of snps
  gamma.W=log(or), # CVs' log(OR)s
  freq=freq_EUR, # Ref haplotypes
  GenoProbList = FP_EUR,
  nrep=n_rep
))

se_sims_AFR <- sqrt(simulated_vbeta(
  N0=afrcontrols, N1=afrcases, # N controls, cases
  snps=ref_leg_AFR[do_sim==T,id], # column names in freq of SNPs for which Z scores should be generated
  W=cvsall$id[i], # CVs, subset of snps
  gamma.W=log(or), # CVs' log(OR)s
  freq=freq_AFR, # Ref haplotypes
  GenoProbList = FP_AFR,
  nrep=n_rep
))

se_sims_EAS <- sqrt(simulated_vbeta(
  N0=eascontrols, N1=eascases, # N controls, cases
  snps=ref_leg_EAS[do_sim==T,id], # column names in freq of SNPs for which Z scores should be generated
  W=cvsall$id[i], # CVs, subset of snps
  gamma.W=log(or), # CVs' log(OR)s
  freq=freq_EAS, # Ref haplotypes
  GenoProbList = FP_EAS,
  nrep=n_rep
))

b_sims_EUR <- z_sims_EUR*se_sims_EUR
p_sims_EUR <- pnorm(-abs(z_sims_EUR))*2

b_sims_AFR <- z_sims_AFR*se_sims_AFR
p_sims_AFR <- pnorm(-abs(z_sims_AFR))*2

b_sims_EAS <- z_sims_EAS*se_sims_EAS
p_sims_EAS <- pnorm(-abs(z_sims_EAS))*2

for (j in 1:n_rep) {
  
  EUR_sumstats <- as.data.frame(matrix(ncol = 9, nrow = ncol(b_sims_EUR)))
  AFR_sumstats <- as.data.frame(matrix(ncol = 9, nrow = ncol(b_sims_AFR)))
  EAS_sumstats <- as.data.frame(matrix(ncol = 9, nrow = ncol(b_sims_EAS)))
  colnames(EUR_sumstats) <- colnames(AFR_sumstats) <- colnames(EAS_sumstats) <- c("chr", "snp", "bp", "a1", "a2", "BETA", "se", "stat", "p")
  EUR_sumstats$chr <- AFR_sumstats$chr <- EAS_sumstats$chr <- 2
  EUR_sumstats$snp <- ref_leg_EUR$id[which(ref_leg_EUR$do_sim==T)]
  EUR_sumstats$bp <- ref_leg_EUR$position[which(ref_leg_EUR$do_sim==T)]
  EUR_sumstats$a1 <- ref_leg_EUR$a0[which(ref_leg_EUR$do_sim==T)]
  EUR_sumstats$a2 <- ref_leg_EUR$a1[which(ref_leg_EUR$do_sim==T)]
  EUR_sumstats$BETA <- as.vector(b_sims_EUR[j,])
  EUR_sumstats$se <- as.vector(se_sims_EUR[j,])
  EUR_sumstats$stat <- as.vector(z_sims_EUR[j,])
  EUR_sumstats$p <- as.vector(p_sims_EUR[j,])
  
  AFR_sumstats$snp <- ref_leg_AFR$id[which(ref_leg_AFR$do_sim==T)]
  AFR_sumstats$bp <- ref_leg_AFR$position[which(ref_leg_AFR$do_sim==T)]
  AFR_sumstats$a1 <- ref_leg_AFR$a0[which(ref_leg_AFR$do_sim==T)]
  AFR_sumstats$a2 <- ref_leg_AFR$a1[which(ref_leg_AFR$do_sim==T)]
  AFR_sumstats$BETA <- as.vector(b_sims_AFR[j,])
  AFR_sumstats$se <- as.vector(se_sims_AFR[j,])
  AFR_sumstats$stat <- as.vector(z_sims_AFR[j,])
  AFR_sumstats$p <- as.vector(p_sims_AFR[j,])
  
  EAS_sumstats$snp <- ref_leg_EAS$id[which(ref_leg_EAS$do_sim==T)]
  EAS_sumstats$bp <- ref_leg_EAS$position[which(ref_leg_EAS$do_sim==T)]
  EAS_sumstats$a1 <- ref_leg_EAS$a0[which(ref_leg_EAS$do_sim==T)]
  EAS_sumstats$a2 <- ref_leg_EAS$a1[which(ref_leg_EAS$do_sim==T)]
  EAS_sumstats$BETA <- as.vector(b_sims_EAS[j,])
  EAS_sumstats$se <- as.vector(se_sims_EAS[j,])
  EAS_sumstats$stat <- as.vector(z_sims_EAS[j,])
  EAS_sumstats$p <- as.vector(p_sims_EAS[j,])
  
  fwrite(AFR_sumstats, paste("AFR_gwas", j, ".txt", sep=""), sep = "\t")
  fwrite(EUR_sumstats, paste("EUR_gwas", j, ".txt", sep=""), sep = "\t")
  fwrite(EAS_sumstats, paste("EAS_gwas", j, ".txt", sep=""), sep = "\t")
}

