gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(minfi)
library(comprehenr)
library(limma)
library(mltools)
library(ggpubr)

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

count.sample <- table(meta.data$tumor_or_normal, meta.data$label) %>% as.data.frame()
colnames(count.sample) <- c("Sample type", "Label", "Count")
count.sample %>% ggplot(aes(x = Label, y = Count, fill = `Sample type`)) + geom_bar(stat = "identity", position = "dodge") + theme_pubr()
