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

if ("caret" %in% installed.packages() == FALSE){
  install.packages("caret")
}
library(caret)

# path.to.project.src <- "/home/hieunguyen/CRC1382/src/TCGA_data_mining"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/TCGA_data_mining"
source(file.path(path.to.project.src, "helper_functions.R"))
infodir <- file.path(path.to.project.src, "TCGA_database_DMR")

# outdir <- "/media/outdir"
outdir <- "/media/hieunguyen/HNSD01/outdir"

PROJECT <- "TCGA_methyl_panel"
output.version <- "20240704"
data.version <- "full"

# path.to.main.input <- "/media/data/TCGA"
path.to.main.input <- "/media/hieunguyen/HNSD01/outdir/TCGA_methyl_panel/input"
path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.07.output <- file.path(path.to.main.output, "07_output_CV")
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### Read input data
#####----------------------------------------------------------------------#####
bVals <- readRDS(file.path(path.to.main.input, data.version, "bVals.rds"))

meta.data <- read.csv(file.path(path.to.project.src, "metadata", "collected_TCGA_metadata.csv")) %>%
  rowwise() %>%
  mutate(sampleid = str_replace(str_replace(File.Name, "_Grn.idat", ""), "_Red.idat", ""))
meta.data <- meta.data[!duplicated(meta.data$sampleid),]
meta.data <- subset(meta.data, meta.data$label != "Ovary")

tumor.samples <- list()
for (input.label in unique(subset(meta.data, meta.data$tumor_or_normal == "tumor")$label)){
  tumor.samples[[input.label]] <- subset(meta.data, meta.data$label == input.label & meta.data$tumor_or_normal == "tumor")$sampleid
}

train.tumor.samples <- list()
test.tumor.samples <- list()

for (fold in seq(1,10)){
  train.tumor.samples[[fold]] <- list()
  test.tumor.samples[[fold]] <- list()
  for (input.label in unique(subset(meta.data, meta.data$tumor_or_normal == "tumor")$label)){
    all.samples <- subset(meta.data, meta.data$label == input.label & meta.data$tumor_or_normal == "tumor")$sampleid
    N <- length(all.samples)
    train.tumor.samples[[fold]][[input.label]] <- sample(all.samples, round(0.8*N))
    test.tumor.samples[[fold]][[input.label]] <- setdiff(all.samples, train.tumor.samples[[fold]][[input.label]])
  }
}


# for (fold in seq(1, 10)){
#   dir.create(file.path(path.to.07.output, sprintf("Fold%s", fold)), showWarnings = FALSE, recursive = TRUE)
#   all.training.samples <- train.tumor.samples[[fold]]
#   print(sprintf("Working on fold %s", fold))
#   for (input.label in names(all.training.samples)){
#     if (file.exists(file.path(path.to.07.output, sprintf("Fold%s", fold), sprintf("test_%s_vs_others.fold%s.rds", input.label, fold))) == FALSE){
#       print(sprintf("working on testing %s vs others", input.label))
#       output <- run_test_group1_vs_group2(input.mat = bVals[, unlist(all.training.samples)], 
#                                           group1 = all.training.samples[[input.label]], 
#                                           group2 = setdiff(unlist(all.training.samples), all.training.samples[[input.label]]))
#       saveRDS(object = output, file.path(path.to.07.output, sprintf("Fold%s", fold), sprintf("test_%s_vs_others.fold%s.rds", input.label, fold)))      
#     }
#   }
# }

run_test <- function(input.label, fold, train.tumor.samples){
  dir.create(file.path(path.to.07.output, sprintf("Fold%s", fold)), showWarnings = FALSE, recursive = TRUE)
  all.training.samples <- train.tumor.samples[[fold]]
  print(sprintf("Working on fold %s", fold))
  if (file.exists(file.path(path.to.07.output, sprintf("Fold%s", fold), sprintf("test_%s_vs_others.fold%s.rds", input.label, fold))) == FALSE){
    print(sprintf("working on testing %s vs others", input.label))
    output <- run_test_group1_vs_group2(input.mat = bVals[, unlist(all.training.samples)],
                                        group1 = all.training.samples[[input.label]],
                                        group2 = setdiff(unlist(all.training.samples), all.training.samples[[input.label]]))
    saveRDS(object = output, file.path(path.to.07.output, sprintf("Fold%s", fold), sprintf("test_%s_vs_others.fold%s.rds", input.label, fold)))
  }  
}

for (fold in seq(1, 10)){
  all.training.samples <- train.tumor.samples[[fold]]
  for (input.label in names(all.training.samples)){
    run_test(input.label, fold, train.tumor.samples)
  }
}
