gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(minfi)
library(comprehenr)
library(limma)
if ("liftOver" %in% installed.packages() == FALSE){
  BiocManager::install("liftOver", update = FALSE)  
}
library(liftOver)

# BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest", update = FALSE)
# install.packages("/media/hieunguyen/HNSD01/Downloads/IlluminaHumanMethylationEPICv2anno.20a1.hg38-master", repos = NULL, type = "source")

maindir <- "/media/hieunguyen/HNSD01"
outdir <- file.path(maindir, "outdir")
PROJECT <- "EPIC"
output.version <- "20240805"
data.version <- "20241007"
path.to.main.input <- file.path(maindir, "storage", PROJECT, data.version)
path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.02.output <- file.path(path.to.main.output, "EPIC_02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

path.to.project.src <- "/media/hieunguyen/HNSD01/src/TCGA_data_mining"
if (file.exists(file.path(path.to.02.output, "idat.obj.rds")) == FALSE){
  idat.obj <- read.metharray.exp(verbose = TRUE, 
                                 force = TRUE, 
                                 base = path.to.main.input)
  
  saveRDS(idat.obj, file.path(path.to.02.output, "idat.obj.rds"))
} else {
  idat.obj <- readRDS(file.path(path.to.02.output, "idat.obj.rds"))
}

idat.obj@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

if (file.exists(file.path(path.to.02.output, "detP.sampling_test.rds")) == FALSE){
  detP <- detectionP(idat.obj)
  saveRDS(detP, file.path(path.to.02.output, "detP.sampling_test.rds"))
} else {
  print("reading in detP ...")
  detP <- readRDS(file.path(path.to.02.output, "detP.sampling_test.rds"))
}

if (file.exists(file.path(path.to.02.output, "idat.obj.sampling_test.preprocessQuantile.rds")) == FALSE){
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

if (file.exists(file.path(path.to.02.output, "check.csv")) == FALSE){
  mVals <- getM(mSetSqFlt)
  bVals <- getBeta(mSetSqFlt)
  
  write.csv(bVals, file.path(path.to.02.output, "bVals.csv"))
  saveRDS(bVals, file.path(path.to.02.output, "bVals.rds"))
  
  write.csv(mVals, file.path(path.to.02.output, "mVals.csv"))
  saveRDS(mVals, file.path(path.to.02.output, "mVals.rds"))
  
  write.csv(data.frame(status = c("finished saving bVals and mVals data to files")),
            file.path(path.to.02.output, "check.csv"))  
} else {
  bVals <- readRDS(file.path(path.to.02.output, "bVals.rds"))
  mVals <- readRDS(file.path(path.to.02.output, "mVals.rds"))
}

atlas <- read.csv(file.path(path.to.project.src, "Figure1_atlas_heatmap.final.csv")) 
colnames(atlas) <- to_vec(for (item in colnames(atlas)) str_replace(item, "X", "")) 
atlas <- atlas %>% column_to_rownames("label")
atlas.cpg <- read.csv(file.path(path.to.project.src, "atlas_cpg.csv")) %>% subset(select = -c(X)) %>%
  rowwise() %>%
  mutate(pos = sprintf("%s_%s", chrom, pos))

if (file.exists(file.path(path.to.02.output, "ann.epic.hg19.csv")) == FALSE){
  ann.epic <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)[, c("chr", "pos")] %>% as.data.frame() %>% rownames_to_column("probe")
  ann.epic <- ann.epic %>% rowwise() %>%
    mutate(end = pos + 1)
  ann.epic <- ann.epic[c("chr", "pos", "end", "probe")]
  ann.epic <- makeGRangesFromDataFrame(df = ann.epic, seqnames.field = "chr", start.field = "pos", end.field = "end", keep.extra.columns = TRUE)
  chain <- import.chain(file.path(path.to.project.src, "hg38ToHg19.over.chain"))
  ann.epic.hg19 <- liftOver(ann.epic, chain) %>% as.data.frame()  %>% subset(select = c(seqnames, start, probe))
  colnames(ann.epic.hg19) <- c("chr", "pos", "probe")
  ann.epic.hg19 <- (ann.epic.hg19) %>% rowwise() %>%
    mutate(pos = sprintf("%s_%s", str_replace(chr, "chr", ""), pos))
  write.csv(ann.epic.hg19, file.path(path.to.02.output, "ann.epic.hg19.csv"))
} else {
  ann.epic.hg19 <- read.csv(file.path(path.to.02.output, "ann.epic.hg19.csv"))
}

ann450k <- read.csv(file.path(path.to.project.src, "ann450kdf.txt"), sep = "\t", header = FALSE) %>%
  rowwise() %>%
  mutate(pos = sprintf("%s_%s", V1, V2))

##### modify the input bVals matrix and the atlas, get shared regions/probes 
# to do deconvolution

intersect.probes <- intersect(ann.epic.hg19$pos, atlas.cpg$pos)
intersect.regions <- subset(atlas.cpg, atlas.cpg$pos %in% intersect.probes)
intersect.regions <- merge(intersect.regions, ann.epic.hg19, by.x = "pos", by.y = "pos")

intersect.bVals <- bVals[intersect(intersect.regions$probe, rownames(bVals)), ] %>% as.data.frame() %>%
  rownames_to_column("probe") 

intersect.bVals <- merge(intersect.bVals, subset(intersect.regions, select = c(probe, region)), by.x = "probe", by.y = "probe") %>%
  column_to_rownames("probe")

intersect.bVals <- intersect.bVals %>% group_by(region) %>%
  summarize(across(everything(), ~ mean(.x, na.rm = TRUE)))

intersect.bVals <- as.data.frame(intersect.bVals) %>% column_to_rownames("region")
tmp.atlas <- atlas[, row.names(intersect.bVals)] %>% t() %>% as.data.frame()
mat <- merge(tmp.atlas %>% rownames_to_column("region"), intersect.bVals %>% rownames_to_column("region"), by.x = "region", by.y = "region") %>%
  column_to_rownames("region")

library(nnls)

tissue.types <- c("Liver", "Breast", "Lung", "Gastric", "CRC", "WBC")
resdf <- data.frame()
for (sample.id in colnames(bVals)){
  nnls.res <- nnls(mat[, tissue.types] %>% as.matrix(), 
                   mat[[sample.id]])
  res <- nnls.res$x/sum(nnls.res$x)
  tmpdf <- data.frame(type = tissue.types, res = res)
  colnames(tmpdf) <- c("type", sample.id)
  tmpdf <- tmpdf %>% column_to_rownames("type") %>% t()
  resdf <- rbind(resdf, tmpdf)
}


writexl::write_xlsx(resdf, file.path(path.to.02.output, "deconvolution_results_8_Breast_samples.xlsx"))


