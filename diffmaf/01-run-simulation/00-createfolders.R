cvsall <- read.table("../cvsall2.txt")
cvsall <- na.omit(cvsall)

for (i in 1:nrow(cvsall)) {
dir.create(paste0("../results/50:50/",cvsall$id[i], "/", "coloc+SuSiEx"), 
           recursive = T)
dir.create(paste0("../data/50:50/",cvsall$id[i], "/", "MseCAVIAR"))
}


for (i in 1:nrow(cvsall)) {
  dir.create(paste0("../results/80:20/",cvsall$id[i], "/", "coloc+SuSiEx"), 
             recursive = T)
  dir.create(paste0("../results/80:20/",cvsall$id[i], "/", "MseCAVIAR"))
}