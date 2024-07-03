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

path.to.project.src <- "/home/hieunguyen/CRC1382/src/TCGA_data_mining"
infodir <- file.path(path.to.project.src, "TCGA_database_DMR")
outdir <- "/media/outdir"

PROJECT <- "TCGA_methyl_panel"
output.version <- "20240704"
data.version <- "full"

path.to.main.input <- "/media/data/TCGA"
path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.07.output <- file.path(path.to.main.output, "07_output")
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)


#####----------------------------------------------------------------------#####
##### Read input data
#####----------------------------------------------------------------------#####
bVals <- readRDS(file.path(path.to.main.input, data.version, "bVals.rds"))
meta.data <- read.csv(file.path(path.to.project.src, "metadata", "collected_TCGA_metadata.csv")) %>%
  rowwise() %>%
  mutate(sampleid = str_replace(str_replace(File.Name, "_Grn.idat", ""), "_Red.idat", ""))
meta.data <- meta.data[!duplicated(meta.data$sampleid),]

tumor.samples <- list()
for (input.label in unique(subset(meta.data, meta.data$tumor_or_normal == "tumor")$label)){
  tumor.samples[[input.label]] <- subset(meta.data, meta.data$label == input.label & meta.data$tumor_or_normal == "tumor")$sampleid
}

normal.samples <- list()
for (input.label in unique(subset(meta.data, meta.data$tumor_or_normal == "normal")$label)){
  normal.samples[[input.label]] <- subset(meta.data, meta.data$label == input.label & meta.data$tumor_or_normal == "normal")$sampleid
}

count.sampledf <- table(meta.data$label, meta.data$tumor_or_normal) %>% as.data.frame() %>%
  pivot_wider(names_from = "Var1", values_from = "Freq") %>% column_to_rownames("Var2") %>% as.data.frame() %>% t() %>% 
  data.frame() %>%
  rownames_to_column("label") %>%
  rowwise() %>%
  mutate(test = ifelse(tumor >= 30, "yes", "no"))

#####----------------------------------------------------------------------#####
##### RUN MAIN TEST BY LIMMA, DML
#####----------------------------------------------------------------------#####
for (group in subset(count.sampledf, count.sampledf$test == "yes")$label){
  print(sprintf("working on group %s, test group1 = %s vs others", group, group))  
  group1 <- tumor.samples[[group]]
  group2 <- setdiff(colnames(bVals), group1)
  
  input.metadata <- data.frame(sample = c(group1, group2), 
                               label = c(
                                 to_vec(for(item in seq(1, length(group1))) "group1"),
                                 to_vec(for(item in seq(1, length(group2))) "group2")
                               ))
  
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  # this is the factor of interest
  g <- factor(input.metadata$label, levels = c("group1", "group2"))
  # use the above to create a design matrix
  design <- model.matrix(~0+label, data=input.metadata)
  colnames(design) <- levels(g)
  fit <- lmFit(bVals[, input.metadata$sample], design)
  # create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(group1-group2,
                              levels=design)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  
  ann450kSub <- ann450k[match(rownames(bVals[, input.metadata$sample]), ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
  DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
  
  
  plot_cpg_two_groups <- function(cpg){
    # cpg <- rownames(DMPs)[1]
    plotdf <- bVals[cpg, input.metadata$sample] %>% as.data.frame()
    plotdf$label <- input.metadata$label
    colnames(plotdf) <- c("cpg", "label")
    p <- plotdf %>% ggplot(aes(x = label, y = cpg, color = label)) + geom_boxplot() + geom_jitter(width = 0.1) + theme_pubr()
    return(p)  
  }
  
  #####----------------------------------------------------------------------#####
  ##### RUN MAIN TEST BY LIMMA, DMR
  #####----------------------------------------------------------------------#####
  myAnnotation <- cpg.annotate(object = bVals[, input.metadata$sample], datatype = "array", what = "Beta", 
                               analysis.type = "differential", design = design, 
                               contrasts = TRUE, cont.matrix = contMatrix, 
                               coef = "group1 - group2", arraytype = "450K")
  
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, min.cpgs = 5)
  results.ranges <- extractRanges(DMRs)
  DMRdf <- data.frame(results.ranges)
  
  writexl::write_xlsx(DMPs, file.path(path.to.07.output, sprintf("DMP_%s.xlsx", group)))
  writexl::write_xlsx(DMRdf, file.path(path.to.07.output, sprintf("DMR_%s.xlsx", group)))
}
