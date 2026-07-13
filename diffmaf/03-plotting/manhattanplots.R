### generates figures 3 and 5 in the main text
library(simGWAS)
library(data.table)
library(readr)
library(qqman)
library(GWASBrewer)


or <- 1.12


load("../../data/simGWAS_freq_EUR.RData")
load("../../data/simGWAS_freq_AFR.RData")

cvsall <- read.table("../cvsall2.txt")
cvsall$or <- 1.12 
ref_leg <- fread("../../data/chr2_maf_gt0.005.leg.gz")
ref_leg_EUR <- ref_leg
ref_leg_AFR <- ref_leg

freq_EUR$Probability <- 1/nrow(freq_EUR)
freq_AFR$Probability <- 1/nrow(freq_AFR)


i <- 1

ref_leg_EUR$do_sim <- NA
ref_leg_AFR$do_sim <- NA
start <- cvsall$pos[i] - 125000
stop <- cvsall$pos[i] + 125000

ref_leg_EUR$do_sim <- ref_leg_EUR$position < stop   &
  ref_leg_EUR$position > start & (ref_leg_EUR$EUR > 0.01 & ref_leg_EUR$EUR < 0.99)  

ref_leg_AFR$do_sim <- ref_leg_AFR$position < stop   &
  ref_leg_AFR$position > start &  (ref_leg_AFR$AFR > 0.01 & ref_leg_AFR$AFR < 0.99)

FP_EUR <- make_GenoProbList(snps=ref_leg_EUR[do_sim==T,id], W=cvsall$id[i], freq=freq_EUR)
LD_EUR <- simGWAS:::wcor2(as.matrix( freq_EUR[,colnames(freq_EUR) %in% ref_leg_EUR[do_sim==T,id]] ), freq_EUR$Probability)
diag(LD_EUR) <- 1.0001

FP_AFR <- make_GenoProbList(snps=ref_leg_AFR[do_sim==T,id], W=cvsall$id[i], freq=freq_AFR)
LD_AFR <- simGWAS:::wcor2(as.matrix( freq_AFR[,colnames(freq_AFR) %in% ref_leg_AFR[do_sim==T,id]] ), freq_AFR$Probability)
diag(LD_AFR) <- 1.0001# LD cannot contain NAs 


generatedata <- function(n, ref_leg, freq, FP) {  
  
  ncases <- n*0.20
  ncontrols <- n * 0.80
  
  z_exp<- expected_z_score(
    N0=ncontrols, N1=ncases, # N controls, cases
    snps=ref_leg[do_sim==T,id], # column names in freq of SNPs for which Z scores should be generated
    W=cvsall$id[i], # CVs, subset of snps
    gamma.W=log(cvsall$or[i]), # CVs' log(OR)s
    freq=freq, # Ref haplotypes
    GenoProbList=FP
  )
  
  se_exp <- sqrt(expected_vbeta(
    N0=ncontrols, N1=ncases,
    snps=ref_leg[do_sim==T,id],
    W=cvsall$id[i],
    gamma.W=log(cvsall$or[i]),
    freq=freq,
    GenoProbList=FP
  ))
  
  
  b_exp <- z_exp*se_exp
  p_exp <- pnorm(-abs(z_exp))*2
  
  sumstats <- as.data.frame(matrix(ncol = 6, nrow = length(b_exp)))
  colnames(sumstats) <- c("CHR", "SNP", "BP", "P", "BETA", "SE")
  sumstats$SNP <- ref_leg$id[which(ref_leg$do_sim==T)]
  sumstats$CHR <- as.numeric(2)
  sumstats$BP <- as.numeric(ref_leg$position[which(ref_leg$do_sim==T)])
  sumstats$P <- as.numeric(as.vector(p_exp))
  sumstats$BETA <- as.numeric(as.vector(b_exp))
  sumstats$SE <- as.numeric(as.vector(se_exp))
  
  return(sumstats)
  
}


sumstats_eur_5050 <- generatedata(n = 50000, ref_leg = ref_leg_EUR, freq = freq_EUR, FP = FP_EUR)
sumstats_eur_8020 <- generatedata(n = 80000, ref_leg = ref_leg_EUR, freq = freq_EUR, FP = FP_EUR)
sumstats_afr_5050 <- generatedata(n = 50000, ref_leg = ref_leg_AFR, freq = freq_AFR, FP = FP_AFR)
sumstats_afr_8020 <- generatedata(n = 20000, ref_leg = ref_leg_AFR, freq = freq_AFR, FP = FP_AFR)

snps <- intersect(sumstats_eur_5050$SNP, sumstats_afr_5050$SNP)

metaanalysis <- function(snps, sumstats_eur, sumstats_afr) {
  eurkeep <- which(sumstats_eur$SNP %in% snps)
  afrkeep <- which(sumstats_afr$SNP %in% snps)
  w_eur <- 1/(sumstats_eur$SE[eurkeep]^2)
  w_afr <- 1/(sumstats_afr$SE[afrkeep]^2)
  beta <- (w_eur*sumstats_eur$BETA[eurkeep] + w_afr*sumstats_afr$BETA[afrkeep])/(w_eur+w_afr)
  se <- sqrt(1/(w_eur+w_afr))
  z <- beta/se
  p <- 2 * pnorm(-abs(z))
  
  metares <- as.data.frame(matrix(ncol = 6, nrow = length(eurkeep)))
  colnames(metares) <- c("CHR", "SNP", "BP", "P", "BETA", "SE")
  metares$SNP <- sumstats_eur$SNP[eurkeep]
  metares$CHR <- as.numeric(2)
  metares$BP <- as.numeric(sumstats_eur$BP[eurkeep])
  metares$P <- as.numeric(as.vector(p))
  metares$BETA <- as.numeric(as.vector(beta))
  metares$SE <- as.numeric(as.vector(se))
  
  return(metares)
}

metagwas5050 <- metaanalysis(snps, sumstats_eur_5050, sumstats_afr_5050)
metagwas8020 <- metaanalysis(snps, sumstats_eur_8020, sumstats_afr_8020)


#par(mfrow = c(3, 2), mar = c(2.5,2.5,1,1), oma = rep(0.5,4), mgp = c(1.5, 0.5, 0))
par(mfrow = c(1,2), mar = c(2.5,2.5,1,1), oma = rep(0.5,4), mgp = c(1.5, 0.5, 0))

manhattan(metagwas5050, xlim = c(start, stop), ylim = c(0,15), highlight = cvsall$id[i])
# mtext("A. 50/50 EUR:AFR", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 
mtext("A.", side = 3, line = 0.2, adj = 0, font = 2, cex = 1)
mtext("50/50 EUR:AFR", side = 3, line = 0.2, adj = 0.5, font = 1, cex = 1)
manhattan(metagwas8020, xlim = c(start, stop), ylim = c(0,15), highlight = cvsall$id[i])
# mtext("B. 80/20 EUR:AFR", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 
mtext("B.", side = 3, line = 0.2, adj = 0, font = 2, cex = 1)
mtext("80/20 EUR:AFR", side = 3, line = 0.2, adj = 0.5, font = 1, cex = 1)

tiff("Figure3.tiff", units="in", width=6, height=3, res=300)

#par(mfrow = c(3, 2), mar = c(2.5,2.5,1,1), oma = rep(0.5,4), mgp = c(1.5, 0.5, 0))
par(mfrow = c(1,2), mar = c(2.5,2.5,1,1), oma = rep(0.5,4), mgp = c(1.5, 0.5, 0))

manhattan(metagwas5050, xlim = c(start, stop), ylim = c(0,15), highlight = cvsall$id[i])
# mtext("A. 50/50 EUR:AFR", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 
mtext("A.", side = 3, line = 0.2, adj = 0, font = 2, cex = 1)
mtext("50/50 EUR:AFR", side = 3, line = 0.2, adj = 0.5, font = 1, cex = 1)
manhattan(metagwas8020, xlim = c(start, stop), ylim = c(0,15), highlight = cvsall$id[i])
# mtext("B. 80/20 EUR:AFR", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 
mtext("B.", side = 3, line = 0.2, adj = 0, font = 2, cex = 1)
mtext("80/20 EUR:AFR", side = 3, line = 0.2, adj = 0.5, font = 1, cex = 1)

dev.off()


tiff("Figure5.tiff", units="in", width=6, height=6, res=300)

par(mfrow = c(2,2), mar = c(2.5,2.5,2,2), oma = c(0.5, 0.5, 0.5, 0.5), mgp = c(1.5, 0.5, 0))
manhattan(sumstats_afr_5050, xlim = c(start, stop), ylim = c(0, 15), highlight = cvsall$id[i], xlab = "") 
#mtext("A. AFR N=50000", side = 3, line = 0.2, adj = 0, font = 2, cex = 0.85) 
mtext("A.", side = 3, line = 0.2, adj = 0, font = 2, cex = 1)
mtext("AFR N=50000", side = 3, line = 0.2, adj = 0.5, font = 1, cex = 1)
#mtext("50/50 EUR:AFR", side = 3, line = 1.2, adj = 0.5, font = 2, cex = 1)
manhattan(sumstats_afr_8020, xlim = c(start, stop), ylim = c(0, 15), highlight = cvsall$id[i], xlab = "")
mtext("B.", side = 3, line = 0.2, adj = 0, font = 2, cex = 1)
mtext("AFR N=20000", side = 3, line = 0.2, adj = 0.5, font = 1, cex = 1)

#mtext("80/20 EUR:AFR", side = 3, line = 1.2, adj = 0.5, font = 2, cex = 1)
manhattan(sumstats_eur_5050, xlim = c(start, stop), ylim = c(0, 15), highlight = cvsall$id[i])
mtext("C.", side = 3, line = 0.2, adj = 0, font = 2, cex = 1)
mtext("EUR N=50000", side = 3, line = 0.2, adj = 0.5, font = 1, cex = 1)
manhattan(sumstats_eur_8020, xlim = c(start, stop), ylim = c(0, 15), highlight = cvsall$id[i])
mtext("D.", side = 3, line = 0.2, adj = 0, font = 2, cex = 1)
mtext("EUR N=80000", side = 3, line = 0.2, adj = 0.5, font = 1, cex = 1)

dev.off()





start <- cvsall$pos[i] - 125000
stop <- cvsall$pos[i] + 125000

LD_AFR <- simGWAS:::wcor2(as.matrix( freq_AFR[,colnames(freq_AFR) %in% ref_leg_AFR[do_sim==T,id]] ), freq_AFR$Probability)
diag(LD_AFR) <- 1.0001# LD cannot contain NAs 
LD_AFR <- list(as.matrix(LD_AFR))

LD_EUR <- simGWAS:::wcor2(as.matrix( freq_EUR[,colnames(freq_EUR) %in% ref_leg_EUR[do_sim==T,id]] ), freq_EUR$Probability)
diag(LD_EUR) <- 1.0001# LD cannot contain NAs 
LD_EUR <- list(as.matrix(LD_EUR))

beta <- 0.046

generateeqtldata <- function(n, ref_leg, LD, af) {
  nrep = 1
  G <- matrix(0, nrow = nrep, ncol = nrep)
  effects <- matrix(0, nrow = length(ref_leg[do_sim==T,id]) , ncol = nrep)
  cvind <- which(ref_leg[do_sim==T,id] == cvsall$id[i])
  effects[cvind, ] <- beta
  effects <- as.matrix(effects)
  sim <- sim_mv_determined(N = n, 
                           direct_SNP_effects_joint = effects, 
                           geno_scale = "allele", 
                           pheno_sd = 1, 
                           G=G, 
                           est_s = T, 
                           R_LD = LD, 
                           af = af)
  
  zexp <- sim$beta_marg / sim$se_beta_hat
  pexp <- 2*pnorm(-abs(zexp))
  
  sumstats <- as.data.frame(matrix(ncol = 6, nrow = length(zexp)))
  colnames(sumstats) <- c("CHR", "SNP", "BP", "P", "BETA", "SE")
  sumstats$SNP <- ref_leg$id[which(ref_leg$do_sim==T)]
  sumstats$CHR <- as.numeric(2)
  sumstats$BP <- as.numeric(ref_leg$position[which(ref_leg$do_sim==T)])
  sumstats$P <- as.numeric(as.vector(pexp))
  sumstats$BETA <- as.numeric(as.vector(sim$beta_marg))
  sumstats$SE <- as.numeric(as.vector(sim$se_beta_hat))
  
  
  return(sumstats)
}


eqtl_eur_5050 <- generateeqtldata(n= 50000, ref_leg = ref_leg_EUR, LD = LD_EUR, af = ref_leg_EUR$EUR[ref_leg_EUR$do_sim==T])
eqtl_eur_8020 <- generateeqtldata(n= 80000, ref_leg = ref_leg_EUR, LD = LD_EUR, af = ref_leg_EUR$EUR[ref_leg_EUR$do_sim==T])
eqtl_afr_5050 <- generateeqtldata(n= 50000, ref_leg = ref_leg_AFR, LD = LD_AFR, af = ref_leg_AFR$AFR[ref_leg_AFR$do_sim==T])
eqtl_afr_8020 <- generateeqtldata(n= 20000, ref_leg = ref_leg_AFR, LD = LD_AFR, af = ref_leg_AFR$AFR[ref_leg_AFR$do_sim==T])

snps <- intersect(eqtl_eur_5050$SNP, eqtl_afr_5050$SNP)
metaeqtl5050 <- metaanalysis(snps, eqtl_eur_5050, eqtl_afr_5050)
metaeqtl8020 <- metaanalysis(snps, eqtl_eur_8020, eqtl_afr_8020)


par(mfrow = c(3, 2), mar = c(2.5,2.5,1,1), oma = rep(0.5,4), mgp = c(1.5, 0.5, 0))
manhattan(metaeqtl5050, xlim = c(start, stop), ylim = c(0,15), highlight = cvsall$id[i])
mtext("A", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 
manhattan(metaeqtl8020, xlim = c(start, stop), ylim = c(0,15), highlight = cvsall$id[i])
mtext("B", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 

manhattan(eqtl_eur_5050, xlim = c(start, stop), ylim = c(0, 15), highlight = cvsall$id[i])
mtext("C", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 
manhattan(eqtl_eur_8020, xlim = c(start, stop), ylim = c(0, 15), highlight = cvsall$id[i])
mtext("D", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 

manhattan(eqtl_afr_5050, xlim = c(start, stop), ylim = c(0, 15), highlight = cvsall$id[i]) 
mtext("E", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 
manhattan(eqtl_afr_8020, xlim = c(start, stop), ylim = c(0, 15), highlight = cvsall$id[i])
mtext("F", side = 3, line = 0.2, adj = 0, font = 2, cex = 1) 
