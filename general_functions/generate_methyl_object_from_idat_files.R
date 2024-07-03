library("dplyr")
library("tidyverse")
library("minfi")
library("comprehenr")
library("limma")
library("illuminaio")
library("BiocParallel")
library("DelayedArray")
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library("IlluminaHumanMethylation450kmanifest")

path.to.input.data <- "/media/hieunguyen/HNSD01/outdir/TCGA_methyl_panel/data_raw/20240620/pool_20240620"
path.to.output <- "/media/hieunguyen/HNSD01/outdir/TCGA_methyl_panel/data_raw/20240620/test_output"

all.files <- Sys.glob(file.path(path.to.input.data, "*.idat"))
all.file.names <- unlist(
  lapply(
    all.files, function(x){
      str_replace(str_replace(basename(x), "_Red.idat", ""), "_Grn.idat", "")
    }
  )
)

all.file.names <- unique(all.file.names)

##### generate pseudo metadata
input.metadata <- data.frame(Basename = unique(all.file.names))

##### generate main data objects
if (file.exists(file.path(path.to.output, "idat.obj.sampling_test.rds")) == FALSE){
  idat.obj <- read.metharray.exp(targets = input.metadata, 
                                 verbose = TRUE, 
                                 force = TRUE, 
                                 base = path.to.input.data)
  
  saveRDS(idat.obj, file.path(path.to.output, "idat.obj.sampling_test.rds"))
  write.csv(sampling.metadata, file.path(path.to.output, "sampling_metadata.csv"))
} else {
  print("reading in idat object...")
  idat.obj <- readRDS(file.path(path.to.output, "idat.obj.sampling_test.rds"))
  sampling.metadata <- read.csv(file.path(path.to.output, "sampling_metadata.csv"))
}

##### DETECTION P VALUE
if (file.exists(file.path(path.to.output, "detP.sampling_test.rds")) == FALSE){
  detP <- detectionP(idat.obj)
  saveRDS(detP, file.path(path.to.output, "detP.sampling_test.rds"))
} else {
  print("reading in detP ...")
  detP <- readRDS(file.path(path.to.output, "detP.sampling_test.rds"))
}

if (file.exists(file.path(path.to.output, "idat.obj.sampling_test.preprocessQuantile.rds")) == FALSE){
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
  saveRDS(mSetSqFlt, file.path(path.to.output, "mSetSqFlt"))
  
} else {
  print("reading preprocessed quantile data, mSetSq")
  mSetSqFlt <- readRDS(file.path(path.to.output, "mSetSqFlt"))
}

mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

write.csv(bVals, file.path(path.to.output, "bVals.csv"))
saveRDS(bVals, file.path(path.to.output, "bVals.rds"))

write.csv(mVals, file.path(path.to.output, "mVals.csv"))
saveRDS(mVals, file.path(path.to.output, "mVals.rds"))