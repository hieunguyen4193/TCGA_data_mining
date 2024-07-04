gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(minfi)
library(comprehenr)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(mltools)
library(DMRcate)
library(ggpubr)
library(RColorBrewer)
library(hash)

# path.to.project.src <- "/home/hieunguyen/CRC1382/src/TCGA_data_mining"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/TCGA_data_mining"

infodir <- file.path(path.to.project.src, "TCGA_database_DMR")

# outdir <- "/media/outdir"
outdir <- "/media/hieunguyen/HNSD01/outdir"

PROJECT <- "TCGA_methyl_panel"
output.version <- "20240704"
data.version <- "full"
thres.depth <- 5

# path.to.main.input <- "/media/data/TCGA"

path.to.main.input <- "/media/hieunguyen/HNSD01/outdir/TCGA_methyl_panel/input"
path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.09.output <- file.path(path.to.main.output, "09_output")
dir.create(path.to.09.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.09.output, sprintf("filtered_cov_%s", thres.depth)), showWarnings = FALSE, recursive = TRUE)
path.to.cov.files <- file.path(path.to.main.input, data.version, "cov_filtered450k")

all.df <- hash()
all.cov.files <- Sys.glob(file.path(path.to.cov.files, "*.cov"))
for (file in all.cov.files){
  print(sprintf("Working on sample %s", basename(file)))
  if (file.exists(file.path(path.to.09.output, sprintf("filtered_cov_%s", thres.depth), basename(file))) == FALSE){
    tmpdf <- read.table(file)[, c("V4", "V5", "V6", "V7", "V8", "V9")]
    colnames(tmpdf) <- c("chrom", "start", "end", "meth", "countC", "countT")
    tmpdf$total <- tmpdf$countC + tmpdf$countT
    tmpdf <- subset(tmpdf, tmpdf$total >= thres.depth)
    write.csv(tmpdf, file.path(path.to.09.output, sprintf("filtered_cov_%s", thres.depth), basename(file)))    
  } else {
    tmpdf <- read.csv(file.path(path.to.09.output, sprintf("filtered_cov_%s", thres.depth), basename(file))) %>%
      rowwise() %>%
      mutate(region_name = sprintf("%s_%s_%s", chrom, start, end)) %>%
      mutate(meth = meth/100)
    all.df[[str_replace(basename(file), ".deduplicated.bedGraph.gz.bismark.zero.filtered450k.cov", "")]] <- tmpdf
  }
}

if (file.exists(file.path(path.to.09.output, "high_depth_bismark.cov.csv")) == FALSE){
  maindf <- all.df[[names(all.df)[[1]]]][, c("region_name", "meth")]
  colnames(maindf) <- c("region_name", names(all.df)[[1]])
  
  for (i in names(all.df)[2:length(names(all.df))]){
    tmpdf <- all.df[[i]][, c("region_name", "meth")]
    colnames(tmpdf) <- c("region_name", i)
    maindf <- merge(maindf, tmpdf, by.x = "region_name", by.y = "region_name")
  }
  write.csv(maindf, file.path(path.to.09.output, "high_depth_bismark.cov.csv"))
} else {
  maindf <- read.csv(file.path(path.to.09.output, "high_depth_bismark.cov.csv"))
}

meta.data <- read.csv("/media/hieunguyen/HNSD01/ECD_metadata/metadata_GW_highdepth.csv")
umap.res <- umap(maindf %>% column_to_rownames("region_name") %>% t() %>% as.matrix())
umapdf <- data.frame(umap.res$layout) %>% rownames_to_column("Sample") %>%
  rowwise() %>%
  mutate(Sample = str_split(Sample, "_")[[1]][[1]]) %>%
  mutate(Sample = ifelse(grepl("-", Sample), str_split(Sample, "-")[[1]][[2]], Sample)) %>%
  mutate(Sample = ifelse(grepl("ME", Sample), str_split(Sample, "ME")[[1]][[1]], Sample)) %>%
  mutate(label = ifelse(length(subset(meta.data, meta.data$SampleID == Sample)$Label) != 0, subset(meta.data, meta.data$SampleID == Sample)$Label, "none"))

umapdf %>% ggplot(aes(x = X1, y = X2, color = label)) + geom_point(size = 3)

