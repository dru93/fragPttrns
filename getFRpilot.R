library(parallel)
library(stringr)

setwd("~/Documents/Andreas/Bordet/Projects/Pearl/")
source("~/Documents/Andreas/Code/scriptsR/fragmentationAnalysis.r")

# get bam file names of all samples
bamPlasmaDirs = paste0(list.files("/Volumes/Transparent/Bioinfo/Data/NIPTpearlPlasma", full.names = T)[1:15], "/aligned")
IDs = sapply(bamPlasmaDirs, function(i) strsplit(i, "/")[[1]][7])
bamFilenames = sapply(IDs, function(x) paste0("/Volumes/Transparent/Bioinfo/Data/NIPTpearlPlasma/",
                                              x, "/aligned/", x, ".sorted.markDup.bam"))
names(bamFilenames) = IDs

# bamFilenames = bamFilenames[1]
AB = getHg19bins100kb()

Sys.time()
mclapply(bamFilenames, function(i){
  ID = strsplit(strsplit(i, "/")[[1]][7], "-PP")[[1]][1]
  outFile = paste0("data/FRpilot100kbNoTrim/", ID, ".RDS")
  fragmentRatio = getFragmentRatioFromBAM(i, AB, bins5Mb = F, print = T, mapqFilt = 30, trimFRs = F,
                                          FRthres = 150) # ~12mins
  saveRDS(fragmentRatio, file = outFile)
}, mc.cores = 4, mc.preschedule = F)
Sys.time()