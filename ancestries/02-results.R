library(MACA)
library(data.table)

cvsall <- read.table("cvsall.txt")
cvsall <- na.omit(cvsall)
ancestries <- "EUR:EAS:AFR"
 split <- "1:1:1"
# split <- "75:15:10"

res <- data.frame(matrix(ncol=17, nrow = 100))
colnames(res) <- c("iter", "cs_incs", "cs_size", "cs_vl", "cs_ll",
                   "es_incs", "es_size", "es_vl", "es_ll",
                   "cm_incs", "cm_size", "cm_vl", "cm_ll",
                   "em_incs", "em_size", "em_vl", "em_ll")

allres <- data.frame(matrix(ncol=17, nrow = nrow(cvsall)))
colnames(allres) <- c("cv", "cs_incs", "cs_size", "cs_vl", "cs_ll",
                      "es_incs", "es_size", "es_vl", "es_ll",
                      "cm_incs", "cm_size", "cm_vl", "cm_ll",
                      "em_incs", "em_size", "em_vl", "em_ll")

for (j in 1:nrow(cvsall)) {
  cv <- cvsall$id[j]
  setwd(paste0("results/", ancestries, "/", split, "/", cv))
  for (i in 1:100) {
    res$iter[i] <- i
    sgwas <- fread(paste0("coloc+SuSiEx/gwas", i, ".SuSiEx.EUR.AFR.EAS.output.cs95.snp"))
    seqtl <- fread(paste0("coloc+SuSiEx/eqtl", i, ".SuSiEx.EUR.AFR.EAS.output.cs95.snp"))
    cs_res <- coloc_susiex(t1=sgwas, t2=seqtl)  
    if (length(cs_res)==3) {
      rowind <- which.max(cs_res$summary$PP.H4.abf) 
      cs_cs <- make_coloc_susiex_credible_set(cs_res, rowind, thresh = 0.95)
      res$cs_incs[i] <- as.numeric(cv %in% cs_cs$snp)
      res$cs_ll[i] <- cs_res$summary$PP.H4.abf[(rowind)]
      k <- rowind + 1
      res$cs_vl[i] <- as.numeric(cs_res$results[(cs_res$results$snp == cv), ..k])
      res$cs_size[i] <- nrow(cs_cs)
    } else {
      res$cs_incs[i] <- res$cs_vl[i] <- res$cs_size[i] <- res$cs_ll[i] <- NA
      
    }
    
    
    es_res <- ecaviar_susiex(t1=sgwas, t2=seqtl)
    if (length(es_res) ==3) {
      rowind <- which.max(es_res$locus_level_clpp)
      es_cs <- make_ecaviar_susiex_credible_set(es_res, rowind, thresh = 0.95)
      res$es_incs[i] <- as.numeric(cv %in% es_cs$snps)
      res$es_ll[i] <- as.numeric(es_res$locus_level_clpp[rowind])
      k <- rowind + 1
      res$es_vl[i] <- as.numeric(es_res$scaled_clpp[which(es_res$scaled_clpp$snps==cv), k])
      res$es_size[i] <- nrow(es_cs)
    } else {
      res$es_incs[i] <- res$es_vl[i] <- res$es_size[i] <- res$es_ll[i] <- NA
    }
    
    mgwas <- fread(paste0("MseCAVIAR/gwas", i, "_", "eqtl", i, ".msecaviar.EUR.AFR.EAS_1_post.txt"))
    meqtl <- fread(paste0("MseCAVIAR/gwas", i, "_", "eqtl", i, ".msecaviar.EUR.AFR.EAS_2_post.txt"))
    
    em_res <- emscaviar(t1=mgwas, t2=meqtl)
    if (length(em_res) == 2) {
      em_cs <- make_emscaviar_credible_set(em_res, thresh = 0.95)
      res$em_incs[i] <- as.numeric(cv %in% em_cs$snps)
      res$em_ll[i] <- as.numeric(em_res$locus_level_clpp)
      res$em_vl[i] <- as.numeric(em_res$variant_level$scaled_clpp[which(em_res$variant_level$snps == cv)])
      res$em_size[i] <- nrow(em_cs)
    }else {
      res$em_incs[i] <- res$em_vl[i] <- res$em_size[i] <- res$em_ll[i] <- NA
    }
    
    cm_res <- coloc_mscaviar(t1=mgwas, t2=meqtl)
    if (length(cm_res)==3) {
      cm_cs <- make_coloc_mscaviar_credible_set(cm_res, thresh = 0.95)
      res$cm_incs[i] <- as.numeric(cv %in% cm_cs$snp)
      res$cm_ll[i] <- as.numeric(cm_res$summary$PP.H4.abf)
      res$cm_vl[i] <- as.numeric(cm_res$results$SNP.PP.H4.abf[which(cm_res$results$snp == cv)])
      res$cm_size[i] <- nrow(cm_cs)
    } else {
      res$cm_incs[i] <- res$cm_vl[i] <- res$cm_size[i] <- res$cm_ll[i] <- NA
    }
    
    
  }
  
  fwrite(res, "results.tsv")
  allres[j,] <- c(cv, sum(res$cs_incs, na.rm = T), median(res$cs_size, na.rm = T), mean(res$cs_vl, na.rm = T), mean(res$cs_ll, na.rm = T),
                  sum(res$es_incs, na.rm = T), median(res$es_size, na.rm = T), mean(res$es_vl, na.rm = T), mean(res$es_ll, na.rm = T),
                  sum(res$cm_incs, na.rm = T), median(res$cm_size, na.rm = T), mean(res$cm_vl, na.rm = T), mean(res$cm_ll, na.rm = T),
                  sum(res$em_incs, na.rm = T), median(res$em_size, na.rm = T), mean(res$em_vl, na.rm = T), mean(res$em_ll, na.rm = T))

}

setwd("../../../results")
write.csv(allres, "1:1:1allres.csv")
