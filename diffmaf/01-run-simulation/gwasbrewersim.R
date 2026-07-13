setwd(paste0("../results/", split, "/", cvsall$id[i], "/coloc+SuSiEx"))


LD_AFR <- simGWAS:::wcor2(as.matrix( freq_AFR[,colnames(freq_AFR) %in% ref_leg_AFR[do_sim==T,id]] ), freq_AFR$Probability)
diag(LD_AFR) <- 1.0001# LD cannot contain NAs 
LD_AFR <- list(as.matrix(LD_AFR))

LD_EUR <- simGWAS:::wcor2(as.matrix( freq_EUR[,colnames(freq_EUR) %in% ref_leg_EUR[do_sim==T,id]] ), freq_EUR$Probability)
diag(LD_EUR) <- 1.0001# LD cannot contain NAs 
LD_EUR <- list(as.matrix(LD_EUR))


### simulate data for AFR 
set.seed(1234)
nrep = 100
G <- matrix(0, nrow = nrep, ncol = nrep) # simulating 100 traits such that there is no causal relationship between them (distinct eQTL studies)
# and they have identical effect sizes, same CV, and all in AFR ancestry
effects <- matrix(0, nrow = length(ref_leg_AFR[do_sim==T,id]) , ncol = nrep)
cvind <- which(ref_leg_AFR[do_sim==T,id] == cvsall$id[i])
effects[cvind, ] <- beta
effects <- as.matrix(effects)
afrsim <- sim_mv_determined(N = rep(afrn,nrep), 
                            direct_SNP_effects_joint = effects, 
                            geno_scale = "allele", 
                            pheno_sd = 1, 
                            G=G, 
                            est_s = T, 
                            R_LD = LD_AFR, 
                            af = ref_leg_AFR$AFR[ref_leg_AFR$do_sim==T])

afrzsim <- afrsim$beta_hat / afrsim$s_estimate
afrpsim <- 2*pnorm(-abs(afrzsim))


### simulate data for EUR
set.seed(1234)
nrep = 100
G <- matrix(0, nrow = nrep, ncol = nrep) # simulating 100 traits such that there is no causal relationship between them (distinct eQTL studies)
# and they have identical effect sizes, same CV, and all in EUR ancestry
effects <- matrix(0, nrow = length(ref_leg_EUR[do_sim==T,id]) , ncol = nrep)
cvind <- which(ref_leg_EUR[do_sim==T,id] == cvsall$id[i])
effects[cvind, ] <- beta
effects <- as.matrix(effects)
eursim <- sim_mv_determined(N = rep(eurn,nrep), 
                            direct_SNP_effects_joint = effects, 
                            geno_scale = "allele", 
                            pheno_sd = 1, 
                            G=G, 
                            est_s = T, 
                            R_LD = LD_EUR, 
                            af = ref_leg_EUR$EUR[ref_leg_EUR$do_sim==T])

eurzsim <- eursim$beta_hat / eursim$s_estimate
eurpsim <- 2*pnorm(-abs(eurzsim))


for (j in 1:nrep) {
  EUR_sumstats <- as.data.frame(matrix(ncol = 9, nrow = length(ref_leg_EUR[do_sim==T,id])))
  AFR_sumstats <- as.data.frame(matrix(ncol = 9, nrow = length(ref_leg_AFR[do_sim==T,id])))
  colnames(EUR_sumstats) <- colnames(AFR_sumstats) <- c("chr", "snp", "bp", "a1", "a2", "BETA", "se", "stat", "p")
  EUR_sumstats$chr <- AFR_sumstats$chr <- 2
  EUR_sumstats$snp <- ref_leg_EUR$id[which(ref_leg_EUR$do_sim==T)]
  EUR_sumstats$bp <- ref_leg_EUR$position[which(ref_leg_EUR$do_sim==T)]
  EUR_sumstats$a1 <- ref_leg_EUR$a0[which(ref_leg_EUR$do_sim==T)]
  EUR_sumstats$a2 <- ref_leg_EUR$a1[which(ref_leg_EUR$do_sim==T)]
  EUR_sumstats$BETA <- as.vector(eursim$beta_hat[,j])
  EUR_sumstats$se <- as.vector(eursim$s_estimate[,j])
  EUR_sumstats$stat <- as.vector(eurzsim[,j])
  EUR_sumstats$p <- as.vector(eurpsim[,j])
  
  AFR_sumstats$snp <- ref_leg_AFR$id[which(ref_leg_AFR$do_sim==T)]
  AFR_sumstats$bp <- ref_leg_AFR$position[which(ref_leg_AFR$do_sim==T)]
  AFR_sumstats$a1 <- ref_leg_AFR$a0[which(ref_leg_AFR$do_sim==T)]
  AFR_sumstats$a2 <- ref_leg_AFR$a1[which(ref_leg_AFR$do_sim==T)]
  AFR_sumstats$BETA <- as.vector(afrsim$beta_hat[,j])
  AFR_sumstats$se <- as.vector(afrsim$s_estimate[,j])
  AFR_sumstats$stat <- as.vector(afrzsim[,j])
  AFR_sumstats$p <- as.vector(afrpsim[,j])
  
  fwrite(AFR_sumstats, paste("AFR_eqtl", j, ".txt", sep=""), sep = "\t")
  fwrite(EUR_sumstats, paste("EUR_eqtl", j, ".txt", sep=""), sep = "\t")
}

#}





