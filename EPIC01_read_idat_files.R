gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(minfi)
library(comprehenr)
library(limma)

BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest")
BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")

maindir <- "/home/hieunguyen/CRC1382/data"
outdir <- file.path(maindir, "outdir")
PROJECT <- "EPIC"
output.version <- "20240805"
data.version <- "20240805"
path.to.main.input <- file.path(maindir, PROJECT, data.version)
path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

if (file.exists(file.path(path.to.01.output, "idat.obj.rds")) == FALSE){
  idat.obj <- read.metharray.exp(verbose = TRUE, 
                                 force = TRUE, 
                                 base = path.to.main.input)
  
  saveRDS(idat.obj, file.path(path.to.01.output, "idat.obj.rds"))
} else {
  idat.obj <- readRDS(file.path(path.to.01.output, "idat.obj.rds"))
}

idat.obj@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

if (file.exists(file.path(path.to.01.output, "detP.sampling_test.rds")) == FALSE){
  detP <- detectionP(idat.obj)
  saveRDS(detP, file.path(path.to.01.output, "detP.sampling_test.rds"))
} else {
  print("reading in detP ...")
  detP <- readRDS(file.path(path.to.01.output, "detP.sampling_test.rds"))
}

if (file.exists(file.path(path.to.01.output, "idat.obj.sampling_test.preprocessQuantile.rds")) == FALSE){
  # filter high detection probes
  keep <- colMeans(detP) < 0.05
  idat.obj <- idat.obj[,keep]
  
  detP <- detP[,keep]
  
  # preprocessing quantile for idat object
  mSetSq <- preprocessQuantile(idat.obj)  
  
  detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
  
  # remove any probes that have failed in one or more samples
  keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
  mSetSqFlt <- mSetSq[keep,]
}
