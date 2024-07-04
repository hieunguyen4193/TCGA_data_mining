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
library(umap)
# library(sva)


# path.to.project.src <- "/home/hieunguyen/CRC1382/src/TCGA_data_mining"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/TCGA_data_mining"

infodir <- file.path(path.to.project.src, "TCGA_database_DMR")

# outdir <- "/media/outdir"
outdir <- "/media/hieunguyen/HNSD01/outdir"

PROJECT <- "TCGA_methyl_panel"
output.version <- "20240704"
data.version <- "full"

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# path.to.main.input <- "/media/data/TCGA"
path.to.main.input <- "/media/hieunguyen/HNSD01/outdir/TCGA_methyl_panel/input"
path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.umap.output <- file.path(path.to.main.output, "UMAP_all_samples")
dir.create(path.to.umap.output, showWarnings = FALSE, recursive = TRUE)

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

if (file.exists(file.path(path.to.project.src, "ann450k_in_atlas.csv")) == FALSE){
  atlascpg <- read.csv(file.path(path.to.project.src, "atlas_cpg.csv")) %>%
    rowwise() %>%
    mutate(probe_coord = sprintf("%s_%s", chrom, pos))
  
  ann450k <- ann450k  %>% as.data.frame() %>% rowwise() %>%
    mutate(probe_coord = sprintf("%s_%s", str_replace(chr, "chr", ""), pos))
  
  ann450k.in.atlas <- subset(ann450k, ann450k$probe_coord %in% atlascpg$probe_coord)
  write.csv(ann450k.in.atlas, file.path(path.to.project.src, "ann450k_in_atlas.csv"))
} else {
  ann450k.in.atlas <- read.csv(file.path(path.to.project.src, "ann450k_in_atlas.csv"))
}

bvals.atlas <- bVals[intersect(ann450k.in.atlas$Name, row.names(bVals)), ]
umap.res <- umap(t(bvals.atlas) %>% as.matrix())
umapdf <- data.frame(umap.res$layout) %>% rownames_to_column("SampleID")
umapdf <- merge(umapdf, subset(meta.data, select = c(sampleid, label, tumor_or_normal)), by.x = "SampleID", by.y = "sampleid")
p.atlas <- umapdf %>% ggplot(aes(x = X1, y = X2, color = label)) + geom_point() + facet_wrap(~tumor_or_normal)

bvals.random <- bVals[sample(row.names(bVals), nrow(bvals.atlas)), ]
umap.res <- umap(t(bvals.random) %>% as.matrix())
umapdf <- data.frame(umap.res$layout) %>% rownames_to_column("SampleID")
umapdf <- merge(umapdf, subset(meta.data, select = c(sampleid, label)), by.x = "SampleID", by.y = "sampleid")
p.random <- umapdf %>% ggplot(aes(x = X1, y = X2, color = label)) + geom_point()

ggarrange(p.atlas, p.random)
# 
# run_pca <- function(inputdf){
#   pc <- prcomp(inputdf,
#                center = TRUE,
#                scale. = TRUE)
#   pcadf <- data.frame(pc$rotation[, 1:2]) %>% as.data.frame() %>% rownames_to_column("SampleID")
#   pcadf <- merge(pcadf, subset(meta.data, select = c(sampleid, label)), by.x = "SampleID", by.y = "sampleid")
#   
#   p <- pcadf %>% ggplot(aes(x = PC1, y = PC2, color = label)) + geom_point()
#   return(p)
# }
# 
# p.atlas <- run_pca(bvals.atlas)
# p.random <- list()
# for (count in seq(1,5)){
#   bvals.random <- bVals[sample(row.names(bVals), nrow(bvals.atlas)), ]
#   p.random[[count]] <- run_pca(bvals.random)
# }
# 
# ggarrange(
#   p.atlas, p.random
# )
# 
# #####
# path.to.main.input <- "/media/hieunguyen/HNSD01/outdir/TCGA_methyl_panel/input"
# path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
# path.to.07.output <- file.path(path.to.main.output, "07_output")
# path.to.08.output <- file.path(path.to.main.output, "08_output")
# 
# group <- "Breast"
# 
# DMPs <- readxl::read_excel(file.path(path.to.07.output, sprintf("DMP_%s.xlsx", group)))
# DMRs <- readxl::read_excel(file.path(path.to.07.output, sprintf("DMR_%s.xlsx", group)))
# 
# DMPs$abs.logFC <- abs(DMPs$logFC)
# DMPs <- subset(DMPs, DMPs$adj.P.Val < 0.001) %>% arrange(desc(abs.logFC))

# 
# library(clusterProfiler)
# get_pathway_descb <- function(i){
#   genes <- str_split(DMRs[i, ]$overlapping.genes, ",")[[1]]
#   genes <- to_vec(for (item in genes) str_replace(item, " ", ""))
#   convertdf <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")    
#   test.res <- data.frame(enrichKEGG(gene         = convertdf$ENTREZID,
#                                     organism     = 'hsa',
#                                     pvalueCutoff = 0.05))
# 
#   if (nrow(test.res) != 0){
#     test.res <- subset(test.res, test.res$p.adjust <= 0.05) %>% arrange(desc(Count))
#     output.str <- paste(head(test.res, 5)$Description, collapse = ", ")  
#   } else {
#     output.str <- ""
#   }
#   return(output.str)
# }
# 
# DMRs$pathway <- unlist(lapply(
#   seq(1, nrow(DMRs)), function(x){
#     return(get_pathway_descb(x))
#   }
# ))



