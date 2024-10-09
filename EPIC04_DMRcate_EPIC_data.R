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
# install.packages("/media/hieunguyen/GSHD_HN01/storage/offline_pkgs/IlluminaHumanMethylationEPICv2anno.20a1.hg38_1.0.0.tar.gz", type = "sources", repos = NULL)

## Create standard colour pallet first
if ("ggthemes" %in% installed.packages() == FALSE){
  install.packages("ggthemes")
}
library("ggthemes")
library("ggrepel")
# install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/DMRcate_3.0.10.tar.gz", type = "source", repos = NULL)
# install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/AnnotationHub_3.12.0.tar.gz", type = "source", repos = NULL)
library(DMRcate)

maindir <- "/media/hieunguyen/HNSD01"
outdir <- file.path(maindir, "outdir")
PROJECT <- "EPIC"
output.version <- "20240805"
data.version <- "20241007"
path.to.storage <- "/media/hieunguyen/GSHD_HN01/storage"
path.to.main.input <- file.path(path.to.storage, PROJECT, data.version)
path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.02.output <- file.path(path.to.main.output, "EPIC_02_output")
path.to.03.output <- file.path(path.to.main.output, "EPIC_03_output")
path.to.04.output <- file.path(path.to.main.output, "EPIC_04_output")
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

resdf <- readxl::read_excel(file.path(path.to.02.output, "deconvolution_results.xlsx"))
bVals <- readRDS(file.path(path.to.02.output, "bVals.rds"))
tissue.types <- c("Liver", "Breast", "Lung", "Gastric", "CRC", "WBC")

resdf$prediction <- unlist(lapply(seq(1, nrow(resdf)), function(i){
  x <- resdf[i, ][tissue.types]
  return(names(x)[x == max(x)])
}))

pbmc.control <- subset(resdf, resdf$SampleType == "PBMC" & resdf$Label == "Control")$SampleCode
pbmc.cancer <- subset(resdf, resdf$SampleType == "PBMC" & resdf$Label != "Control")$SampleCode
pbmc.crc <- subset(resdf, resdf$SampleType == "PBMC" & resdf$Label == "CRC-cancer")$SampleCode
pbmc.lung <- subset(resdf, resdf$SampleType == "PBMC" & resdf$Label == "Lung-cancer")$SampleCode

breast.cancer <- subset(resdf, resdf$SampleType == "Tissue" & resdf$Label == "Breast-cancer")$SampleCode
breast.control <- subset(resdf, resdf$SampleType == "Tissue" & resdf$Label != "Breast-cancer")$SampleCode

#####---------------------------------------------------------######
##### RUN DMR
#####---------------------------------------------------------######
## Define your variable of interest and convert to a factor
group1 <- breast.cancer
group2 <- breast.control

input.metadata <- data.frame(sample = c(group1, group2), 
                             label = c(
                               to_vec(for(item in seq(1, length(group1))) "cancer"),
                               to_vec(for(item in seq(1, length(group2))) "control")))
                               
variable.of.choice <- factor(input.metadata$label, level = c("control", "cancer"))

myColors <- gdocs_pal()(length(unique(variable.of.choice)))

names(myColors) <- levels(factor(levels = unique(variable.of.choice)))
colScale <- scale_colour_manual(
  name = "Cell Line",
  values = myColors,
  drop = T,
  limits = force
)
##### FILTER
betaMat.clean <- rmSNPandCH(bVals[, input.metadata$sample])

## Note: If you would like to keep all of the SNPs then you can set parameters in the above code as follows. rmSNPandCH(beta_v2.detP2, dist=0, mafcut = 1)

## We can also remove the replicate probes and collapse them based off preference. Mean, sensitivity, specificity and random are available options. We will use mean in this scenario

#####---------------------------------------------------------######
##### density plot
#####---------------------------------------------------------#####
x <- data.table::melt(betaMat.clean)

x$Variable <- variable.of.choice[match(x$Var2, input.metadata$sample)]

p1 <- ggplot(x, aes(x = value, color = Variable)) +
  geom_density() +
  theme_minimal() +
  ylab("Density") +
  theme(legend.position = "bottom") +
  colScale +
  ggtitle("Density Plot by Group Average")

p2 <- ggplot(x, aes(x = value, color = Var2)) +
  geom_density() +
  theme_minimal() +
  ylab("Density") +
  theme(legend.position = "bottom") +
  ggtitle("Density Plot by Sample")

plot <- limma::plotMDS(
  betaMat.clean,
  top = nrow(betaMat.clean),
  cex = 0.5,
  main = "All Samples",
  labels = colnames(betaMat.clean),
  plot = F
)

toplot <- data.frame(
  distance.matrix = plot$distance.matrix,
  x = plot$x,
  y = plot$y
)

ggplot(
  toplot$distance.matrix,
  aes(
    x = toplot$x,
    y = toplot$y,
    colour = variable.of.choice,
    label = input.metadata$sample
  )
) +
  ggtitle("MDS - All Probes") +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  ylab("Dim 2") + xlab("Dim 1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom") +
  colScale +
  geom_text_repel(min.segment.length = 0, box.padding = 1)

#####---------------------------------------------------------#####
##### MDS
#####---------------------------------------------------------#####
plot <- limma::plotMDS(
  betaMat.clean,
  top = nrow(betaMat.clean),
  cex = 0.5,
  main = "All Samples",
  labels = colnames(betaMat.clean),
  plot = F
)

toplot <- data.frame(
  distance.matrix = plot$distance.matrix,
  x = plot$x,
  y = plot$y
)

ggplot(
  toplot$distance.matrix,
  aes(
    x = toplot$x,
    y = toplot$y,
    colour = variable.of.choice,
    label = input.metadata$sample
  )
) +
  ggtitle("MDS - All Probes") +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  ylab("Dim 2") + xlab("Dim 1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom") +
  colScale +
  geom_text_repel(min.segment.length = 0, box.padding = 1)

#####---------------------------------------------------------#####
##### DMR
#####---------------------------------------------------------#####
design <- model.matrix( ~ variable.of.choice)
myannotation <- cpg.annotate(
  datatype = "array",
  object = betaMat.clean,
  arraytype = "EPICv2",
  epicv2Remap = FALSE,
  analysis.type = "differential",
  design = design,
  coef = 2, 
  what = "Beta", 
  epicv2Filter = "mean"
)

myannotation <- changeFDR(myannotation, 1e-10)
dmrcoutput <- dmrcate(myannotation)

results.ranges <- extractRanges(dmrcoutput = dmrcoutput, genome = "hg38")
results.rangesdf <- data.frame(results.ranges) %>% rowwise() %>%
  mutate(abs.meandiff = abs(meandiff)) %>% 
  arrange(desc(abs.meandiff)) %>%
  separate_rows(overlapping.genes)

#####---------------------------------------------------------#####
##### ORA KEGG AND GO
#####---------------------------------------------------------#####
if ("clusterProfiler" %in% installed.packages() == FALSE){
  BiocManager::install("clusterProfiler", update = FALSE)
}
library("clusterProfiler")

sig.genes <- head(results.rangesdf, 1000)$overlapping.genes 
sig.genes <- sig.genes[is.na(sig.genes) == FALSE]

convertdf <- bitr(results.ranges$overlapping.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

go.res <- enrichGO(universe = convertdf$ENTREZID,
                   gene = subset(convertdf, convertdf$SYMBOL %in% sig.genes)$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable = TRUE)
kegg.res <- enrichKEGG(gene = subset(convertdf, convertdf$SYMBOL %in% sig.genes)$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)

dotplot(kegg.res, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(go.res, showCategory=30) + ggtitle("dotplot for ORA")

#####---------------------------------------------------------#####
##### plot
#####---------------------------------------------------------#####
# cols <- myColors[as.character(variable.of.choice)]
# DMR.plot(
#   ranges = results.ranges,
#   dmr = 1,
#   CpGs = betaMat.clean,
#   what = "Beta",
#   arraytype = "EPICv2",
#   genome = "hg38",
#   phen.col = cols
# )
