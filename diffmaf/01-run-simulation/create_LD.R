
setwd(paste0("../results/", split, "/", cvsall$id[i], "/coloc+SuSiEx"))
  
snps_AFR <- ref_leg_AFR[do_sim==T, id]
snps_EUR <- ref_leg_EUR[do_sim==T, id]


writeLines(snps_AFR, "snps_AFR_to_extract.txt")
writeLines(snps_EUR, "snps_EUR_to_extract.txt")


system(paste("/plink2", 
             "--haps   ../../data/chr2_maf_gt0.005.hap.gz",
             "--legend ../../data/chr2_maf_gt0.005.leg.gz 2", # *<!>* Don't forget the chromosome code! 
             "--sample ../../data/1000GP_Phase3.sample2", 
             "--keep", "../../data/sample_ids-AFR.txt",
             "--extract snps_AFR_to_extract.txt",
             "--make-bed",
             "--out AFR"
))

system(paste("/plink2", 
             "--haps   ../../data/chr2_maf_gt0.005.hap.gz",
             "--legend ../../data/chr2_maf_gt0.005.leg.gz 2", # *<!>* Don't forget the chromosome code! 
             "--sample ../../data/1000GP_Phase3.sample2", 
             "--keep", "../../data/sample_ids-EUR.txt",
             "--extract snps_EUR_to_extract.txt",
             "--make-bed",
             "--out EUR"
))



#}









