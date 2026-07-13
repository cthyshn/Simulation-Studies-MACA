
# 0. global parameters and packages --------------------------------------
#this code simulates data for selected regions, according to the SuSiEx format
#that is, I have not filtered out SNPs that are present in only one ancestry but not in others
library(simGWAS)
library(GWASBrewer)
library(data.table)
library(readr)
library(qqman)
library(doParallel)
library(foreach)


load("../../data/simGWAS_freq_EUR.RData")
load("../../data/simGWAS_freq_AFR.RData")

ref_leg <- fread("../../data/chr2_maf_gt0.005.leg.gz")
ref_sam <- read.table("../../data/1000GP_Phase3.sample2")

freq_EUR$Probability <- 1/nrow(freq_EUR)
freq_AFR$Probability <- 1/nrow(freq_AFR)


# 1. Analysis  --------------------------------------


cvsall <- read.table("../cvsall2.txt")
cvsall <- na.omit(cvsall)
split <- "80:20" # or "50:50"

n <- 100000


for (i in 1:nrow(cvsall)){
  print(i)
  

    if (split == "50:50") {
      afrn <- 0.5*n
      eurn <- 0.5*n
    }
    if (split == "80:20") {
      afrn <- 0.2*n
      eurn <- 0.8*n
    }
  
  
  afrcases <- afrn*0.2
  eurcases <- eurn*0.2
  afrcontrols <- afrn*0.8
  eurcontrols <- eurn*0.8
  or <- 1.12
  beta <- 0.046



    source("simgwassim.R")
    source("gwasbrewersim.R")
    source("create_LD.R")
    source("msecaviar_sumstat.R")
    source("createscripts.R")
  
}
  
  


