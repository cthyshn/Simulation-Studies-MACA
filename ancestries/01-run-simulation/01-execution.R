
# 0. global parameters and packages --------------------------------------
#this code simulates data for selected regions, according to the SuSiEx format
#that is, I have not filtered out SNPs that are present in only one ancestry but not in others
library(simGWAS)
library(GWASBrewer)
library(data.table)
library(readr)
library(qqman)

load("../../data/simGWAS_freq_EUR.RData")
load("../../data/simGWAS_freq_AFR.RData")
load("../../data/simGWAS_freq_EAS.RData")

ref_leg <- fread("../../data/chr2_maf_gt0.005.leg.gz")
ref_sam <- read.table("../../data/1000GP_Phase3.sample2")

freq_EUR$Probability <- 1/nrow(freq_EUR)
freq_AFR$Probability <- 1/nrow(freq_AFR)
freq_EAS$Probability <- 1/nrow(freq_EAS)

# 1. Analysis  --------------------------------------

cvsall <- read.table("cvsall.txt")
cvsall <- na.omit(cvsall)
ancestries <- "EUR:EAS:AFR"
split <- "75:15:10" # "1:1:1"
n <- 100000


for (i in 1:nrow(cvsall)){
  print(i)
  
  if (ancestries == "EUR:EAS:AFR") {
    if (split == "1:1:1" & n == 100000) {
      eurcases <- afrcases <- eascases <- 6667
      afrcontrols <- eascontrols <- 26666
      eurcontrols <- 26667
      eurn <- 33334
      afrn <- easn <- 33333
      }else {
      afrn <- easn <- eurn <- as.integer(n/3)}
    if (split == "75:15:10") {
      afrn <- 0.1*n
      easn <- 0.15*n
      eurn <- 0.75*n
      eascases <- easn*0.2
      eascontrols <- easn*0.8
      afrcases <- afrn*0.2
      eurcases <- eurn*0.2
      afrcontrols <- afrn*0.8
      eurcontrols <- eurn*0.8
    }
    }
  
  or <- 1.12
  beta <- 0.046

  if (ancestries == "EUR:EAS:AFR") {
    source("simgwassim.R")
    source("gwasbrewersim.R")
    source("create_LD.R")
    source("msecaviar_sumstat.R")
    source("createscripts.R")
  }
  }
  
  

