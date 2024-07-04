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
library(heatmaply)
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
path.to.heatmap.output <- file.path(path.to.main.output, "heatmap")
dir.create(path.to.heatmap.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.input <- "/media/hieunguyen/HNSD01/outdir/TCGA_methyl_panel/input"
path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.07.output <- file.path(path.to.main.output, "07_output")
#####----------------------------------------------------------------------#####
##### Read input data
#####----------------------------------------------------------------------#####
bVals <- readRDS(file.path(path.to.main.input, data.version, "bVals.rds"))
meta.data <- read.csv(file.path(path.to.project.src, "metadata", "collected_TCGA_metadata.csv")) %>%
  rowwise() %>%
  mutate(sampleid = str_replace(str_replace(File.Name, "_Grn.idat", ""), "_Red.idat", ""))
meta.data <- meta.data[!duplicated(meta.data$sampleid),]

topN <- 500
if (file.exists(file.path(path.to.heatmap.output, sprintf("hypo_cpgs_tumor_samples.top%s.rds", topN))) == FALSE){
  DMPs <- list()
  hypo.cpgs <- list()
  hyper.cpgs <- list()
  for (group in setdiff(unique(meta.data$label), c("Ovary"))){
    print(sprintf("working on %s", group))
    tmpdf <- readxl::read_excel(file.path(path.to.07.output, sprintf("DMP_%s.xlsx", group)))  
    tmpdf <- subset(tmpdf, tmpdf$adj.P.Val <= 0.01)
    DMPs[[group]] <- tmpdf
    hyper.cpgs[[group]] <- subset(tmpdf, tmpdf$logFC >= 0) %>% arrange(desc(logFC)) %>% head(topN) %>% pull(Name)
    hypo.cpgs[[group]] <- subset(tmpdf, tmpdf$logFC <= 0) %>% arrange(desc(logFC)) %>% tail(topN) %>% pull(Name)
  }
  saveRDS(hyper.cpgs, file.path(path.to.heatmap.output, sprintf("hyper_cpgs_tumor_samples.top%s.rds", topN)))
  saveRDS(hypo.cpgs, file.path(path.to.heatmap.output, sprintf("hypo_cpgs_tumor_samples.top%s.rds", topN)))  
} else {
  hyper.cpgs <- readRDS(file.path(path.to.heatmap.output, sprintf("hyper_cpgs_tumor_samples.top%s.rds", topN)))
  hypo.cpgs <- readRDS(file.path(path.to.heatmap.output, sprintf("hypo_cpgs_tumor_samples.top%s.rds", topN)))
}

diff.hyper.cpgs <- data.frame()
for (group in names(hyper.cpgs)){
  tmpdf <- data.frame(probe = hyper.cpgs[[group]])
  tmpdf$label <- group
  diff.hyper.cpgs <- rbind(diff.hyper.cpgs, tmpdf)
}

diff.hypo.cpgs <- data.frame()
for (group in names(hypo.cpgs)){
  tmpdf <- data.frame(probe = hypo.cpgs[[group]])
  tmpdf$label <- group
  diff.hypo.cpgs <- rbind(diff.hypo.cpgs, tmpdf)
}

filterTSMA.diff.hyper.cpgs <- subset(diff.hyper.cpgs, diff.hyper.cpgs$label %in% c("Breast", "Colon", "Gastric", "Liver", "Lung", "Rectum"))
filterTSMA.diff.hypo.cpgs <- subset(diff.hypo.cpgs, diff.hypo.cpgs$label %in% c("Breast", "Colon", "Gastric", "Liver", "Lung", "Rectum"))

# order.samples <- subset(meta.data, meta.data$tumor_or_normal == "tumor") %>% arrange(desc(label)) %>% pull(sampleid)
# input.heatmap.hyper <- bVals[diff.hyper.cpgs$probe, order.samples]
# input.heatmap.hypo <- bVals[diff.hypo.cpgs$probe, order.samples]

# pheatmap::pheatmap(input.heatmap.hyper, show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE)
# pheatmap::pheatmap(input.heatmap.hypo, show_rownames = FALSE, show_colnames = FALSE, cluster_rows = FALSE, cluster_cols = FALSE)

moss.atlas <- readxl::read_excel(file.path(path.to.project.src, "41467_2018_7466_MOESM3_ESM.xlsx"), sheet = "Table 1")
filterTSMA.moss.atlas <- subset(moss.atlas, moss.atlas$name %in% c("Hepatocytes", "Lung cells", "Colon epithelial cells", "Breast", "Upper GI" ))
ann450k.in.tsma <- read.csv(file.path(path.to.project.src, "ann450k_in_atlas.csv"))
library(VennDiagram)

venn.diagram(
  x = list(unique(moss.atlas$acc), 
           unique(c(filterTSMA.diff.hyper.cpgs$probe, filterTSMA.diff.hypo.cpgs$probe)), 
           unique(ann450k.in.tsma$Name)),
  category.names = c("Moss atlas" , "TCGA " , "TSMA"),
  filename = file.path(path.to.project.src, "venn_diagram_all_probes.png"),
  output=TRUE, disable.logging = TRUE 
)
