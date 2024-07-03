gc()
rm(list = ls())

new.pkgs <- c("minfi", 
              "IlluminaHumanMethylation450kanno.ilmn12.hg19",
              "IlluminaHumanMethylation450kmanifest",
              "DMRcate",
              "missMethyl", 
              "HDF5Array")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    BiocManager::install(pkg, update = FALSE)
  }
}

library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)
library(dplyr)
library(tidyverse)

if ("methylationArrayAnalysis" %in% installed.packages() == FALSE){
  BiocManager::install("methylationArrayAnalysis", update = FALSE)
}

inputdir <- "/media/hieunguyen/GSHD_HN01/idat_TCGA"
outdir <- "/media/hieunguyen/HNSD01/outdir"

output.version <- "20240612"
path.to.main.input <- file.path(outdir, sprintf("TCGA_%s", output.version), sprintf("idat_%s", output.version))
# path.to.main.input <- "/media/hieunguyen/HNSD01/outdir/TCGA_20240612/test_dir"
path.to.main.output <- file.path(outdir, sprintf("TCGA_%s", output.version))
path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)
meta.data <- read.csv(file.path(inputdir, "metadata", output.version, "gdc_sample_sheet.2024-06-11.idat_only.primary_tumor_normal_only.tsv"))

raw.arraytypedf <- read.csv(file.path(path.to.01.output, "array_typedf.csv")) 
arraytypedf <- raw.arraytypedf %>%
  subset(arrayType == "IlluminaHumanMethylation450k") %>%
  rowwise() %>%
  mutate(base.name = str_replace(basename(path), "_Grn.idat", ""))

meta.data <- meta.data %>%
  mutate(base.name = str_replace(str_replace(File.Name, "_Grn.idat", ""), "_Red.idat", ""))

meta.data <- subset(meta.data, meta.data$base.name %in% arraytypedf$base.name) %>%
  rowwise() %>%
  mutate(c = ifelse(grepl("Grn", File.Name), "Green", "Red")) %>%
  mutate(true_sampleid = sprintf("%s_%s_%s", Case.ID, Sample.ID, Sample.Type))

count.channels.samples <- table(meta.data$true_sampleid)
meta.data <- subset(meta.data, meta.data$true_sampleid %in% names(count.channels.samples[count.channels.samples == 2]))

if (file.exists(file.path(path.to.01.output, "bVals.rds")) == FALSE){
  meta.data.normal <- subset(meta.data, meta.data$Sample.Type == "Solid Tissue Normal")
  meta.data.tumor <- subset(meta.data, meta.data$Sample.Type == "Primary Tumor")
  sampling.metadata.tumor <- data.frame()
  for (label in unique(meta.data.tumor$Project.ID)){
    all.samples <- subset(meta.data.tumor, meta.data.tumor$Project.ID == label)$true_sampleid
    if (length(all.samples) >= 200){
      selected.samples <- sample(all.samples, 200)
    } else {
      selected.samples <- all.samples 
    }
    sampling.metadata.tumor <- rbind(sampling.metadata.tumor, subset(meta.data.tumor, meta.data.tumor$Project.ID == label & true_sampleid %in% selected.samples))
  }
  sampling.metadata <- rbind(sampling.metadata.tumor, meta.data.normal)
  if (file.exists(file.path(path.to.01.output, "idat.obj.sampling_test.rds")) == FALSE){
    idat.obj <- read.metharray.exp(targets = data.frame(Basename = unique(sampling.metadata$base.name)), verbose = TRUE, force = TRUE, base = path.to.main.input)
    
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
  print("reading in data ...")
  bVals <- readRDS(file.path(path.to.01.output, "bVals.rds"))
  sampling.metadata <- read.csv(file.path(path.to.01.output, "sampling_metadata.csv"))
}

