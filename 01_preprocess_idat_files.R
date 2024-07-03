gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(minfi)
library(comprehenr)
library(limma)

infodir <- "/media/hieunguyen/HNSD01/src/TCGA_data_analysis/TCGA_database_DMR"
outdir <- "/media/hieunguyen/HNSD01/outdir"
raw.data.dir <- "/media/hieunguyen/HNHD01/TCGA_all_idat"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/TCGA_data_analysis"

source(file.path(path.to.project.src, "check_methylation_array_type.R"))

PROJECT <- "TCGA_methyl_panel"
output.version <- "20240620"
data.version <- "raw"

sampling.version <- "sampling_20240620"

path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.01.output <- file.path(path.to.main.output, "01_output", sampling.version)
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.metadata.dir <- file.path(path.to.main.output, "metadata")
meta.data <- read.csv(file.path(path.to.metadata.dir, "collected_TCGA_metadata.csv"))

path.to.pool.data <- file.path(path.to.main.output, sprintf("pool_%s", output.version))

meta.data.normal <- subset(meta.data, meta.data$tumor_or_normal == "normal")
meta.data.tumor <- subset(meta.data, meta.data$tumor_or_normal == "tumor")

sample.list <- list()
for (input.label in unique(meta.data.tumor$label)){
  sample.list[[input.label]] <- unique(subset(meta.data.tumor, meta.data.tumor$label == input.label)$Sample.ID)
}

if (file.exists(file.path(path.to.01.output, "bVals.csv")) == FALSE){
  if (file.exists(file.path(path.to.01.output, sprintf("metadata.%s.csv", sampling.version))) == FALSE){
    sampling.data <- list()
    for (i in names(sample.list)){
      if (length(sample.list[[i]]) >= 100){
        sampling.data <- c(sampling.data, unlist(sample(sample.list[[i]], 100)))    
      } else {
        sampling.data <- c(sampling.data, sample.list[[i]])
      }
    }
    sampling.data <- unlist(sampling.data)
    sampling.metadata <- rbind(meta.data.normal, 
                               subset(meta.data.tumor, meta.data.tumor$Sample.ID %in% sampling.data))
    write.csv(sampling.metadata, file.path(path.to.01.output, sprintf("metadata.%s.csv", sampling.version)))
  } else {
    sampling.metadata <- read.csv(file.path(path.to.01.output, sprintf("metadata.%s.csv", sampling.version)))
  }
  
  if (file.exists(file.path(path.to.01.output, "idat.obj.sampling_test.rds")) == FALSE){
    idat.obj <- read.metharray.exp(targets = data.frame(Basename = unique(basename(sampling.metadata$base.path))), 
                                   verbose = TRUE, 
                                   force = TRUE, 
                                   base = path.to.pool.data)
    
    saveRDS(idat.obj, file.path(path.to.01.output, "idat.obj.sampling_test.rds"))
    write.csv(sampling.metadata, file.path(path.to.01.output, "sampling_metadata.csv"))
  } else {
    print("reading in idat object...")
    idat.obj <- readRDS(file.path(path.to.01.output, "idat.obj.sampling_test.rds"))
    sampling.metadata <- read.csv(file.path(path.to.01.output, "sampling_metadata.csv"))
  }
  
  ##### DETECTION P VALUE
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
    
    # get 450k CpG annotation    
    ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    
    # remove CpGs/probe in sex chromosomes
    keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
    mSetSqFlt <- mSetSqFlt[keep,]
    
    # remove probes with SNPs at CpG site
    mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
    
    # exclude cross reactive probes 
    dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
    xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                           "48639-non-specific-probes-Illumina450k.csv",
                                           sep="/"), stringsAsFactors=FALSE)
    keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
    mSetSqFlt <- mSetSqFlt[keep,] 
    
    # save object mSetSqFlt
    saveRDS(mSetSqFlt, file.path(path.to.01.output, "mSetSqFlt"))
    
  } else {
    print("reading preprocessed quantile data, mSetSq")
    mSetSqFlt <- readRDS(file.path(path.to.01.output, "mSetSqFlt"))
  }
  
  mVals <- getM(mSetSqFlt)
  bVals <- getBeta(mSetSqFlt)
  
  write.csv(bVals, file.path(path.to.01.output, "bVals.csv"))
  saveRDS(bVals, file.path(path.to.01.output, "bVals.rds"))
  
  write.csv(mVals, file.path(path.to.01.output, "mVals.csv"))
  saveRDS(mVals, file.path(path.to.01.output, "mVals.rds"))
} else {
  print(sprintf("reading in bVals object ... at %s", file.path(path.to.01.output, "bVals.rds")))
  bVals <- readRDS(file.path(path.to.01.output, "bVals.rds"))
  sampling.metadata <- read.csv(file.path(path.to.01.output, sprintf("metadata.%s.csv", sampling.version)))
}

