cvsall <- read.table("../cvsall.txt")
cvsall <- na.omit(cvsall)

for (i in 1:nrow(cvsall)) {
dir.create(paste0("../results/EUR:EAS:AFR/75:15:10/",cvsall$id[i], "/", "coloc+SuSiEx"), 
           recursive = T)
dir.create(paste0("../results/EUR:EAS:AFR/75:15:10/",cvsall$id[i], "/", "MseCAVIAR"))
}

