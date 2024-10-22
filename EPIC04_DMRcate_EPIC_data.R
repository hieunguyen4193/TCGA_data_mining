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

#####----------------------------------------------------------------------#####
# install some new packages
# BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest", update = FALSE)
# install.packages("/media/hieunguyen/GSHD_HN01/storage/offline_pkgs/IlluminaHumanMethylationEPICv2anno.20a1.hg38_1.0.0.tar.gz",
#                  type = "sources", repos = NULL)
# BiocManager::install("methylationArrayAnalysis", update = FALSE)
# install.packages("nnls")
# install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/DMRcate_3.0.10.tar.gz", type = "source", repos = NULL)
# install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/AnnotationHub_3.12.0.tar.gz", type = "source", repos = NULL)
#####----------------------------------------------------------------------#####

## Create standard colour pallet first
if ("ggthemes" %in% installed.packages() == FALSE){
  install.packages("ggthemes")
}
library("ggthemes")
library("ggrepel")
library(DMRcate)

maindir <- "/media/hieunguyen/HNSD01"
outdir <- file.path(maindir, "outdir")
PROJECT <- "EPIC"
output.version <- "20241020"
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
mVals <- readRDS(file.path(path.to.02.output, "mVals.rds"))
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

perform_DMR_with_DMRcate <- function(group1, group2, inputdf, outputdir, pval_thres, beta_or_M, run_GO_KEGG = FALSE){
  dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
  input.metadata <- data.frame(sample = c(group1, group2), 
                               label = c(
                                 to_vec(for(item in seq(1, length(group1))) "cancer"),
                                 to_vec(for(item in seq(1, length(group2))) "control")))
  
  variable.of.choice <- factor(input.metadata$label, levels = c("control", "cancer"))
  
  myColors <- gdocs_pal()(length(unique(variable.of.choice)))
  
  names(myColors) <- levels(factor(levels = c("control", "cancer")))
  colScale <- scale_colour_manual(
    name = "Cell Line",
    values = myColors,
    drop = T,
    limits = force
  )
  ##### FILTER
  betaMat.clean <- rmSNPandCH(inputdf[, input.metadata$sample])
  
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
    what = beta_or_M, 
    epicv2Filter = "mean"
  )
  
  myannotation <- changeFDR(myannotation, pval_thres)
  dmrcoutput <- dmrcate(myannotation)
  
  results.ranges <- extractRanges(dmrcoutput = dmrcoutput, genome = "hg38")
  results.rangesdf <- data.frame(results.ranges) %>% rowwise() %>%
    mutate(abs.meandiff = abs(meandiff)) %>% 
    arrange(desc(abs.meandiff)) %>%
    separate_rows(overlapping.genes)
  
  writexl::write_xlsx(results.rangesdf, file.path(outputdir, "DMR_results.xlsx"))
  print(sprintf("DMR result table is saved at %s", file.path(outputdir, "DMR_results.xlsx")))
  
  #####---------------------------------------------------------#####
  ##### ORA KEGG AND GO
  #####---------------------------------------------------------#####
  if (run_GO_KEGG == TRUE){
    if ("clusterProfiler" %in% installed.packages() == FALSE){
      BiocManager::install("clusterProfiler", update = FALSE)
    }
    library("clusterProfiler")
    
    sig.genes <- head(results.rangesdf, 2000)$overlapping.genes 
    sig.genes <- sig.genes[is.na(sig.genes) == FALSE]
    sig.genes <- sig.genes[is.na(sig.genes) == FALSE]
    
    all.genes <- results.ranges$overlapping.genes
    all.genes <- all.genes[is.na(all.genes) == FALSE]
    convertdf <- bitr(all.genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    
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
    write.csv(data.frame(kegg.res), file.path(outputdir, "KEGG_res.csv"))
    write.csv(data.frame(go.res), file.path(outputdir, "GO_res.csv"))
    dotplot.kegg <- dotplot(kegg.res, showCategory=30) + ggtitle("dotplot for ORA")
    ggsave(plot = dotplot.kegg, filename = "dotplot_KEGG_DMR_genes.svg", path = outputdir, device = "svg", width = 14, height = 10)
    print(sprintf("DMR KEGG dotplot is saved at %s", file.path(outputdir, "dotplot_KEGG_DMR_genes.svg")))
    
  }
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
}

perform_DMR_with_DMRcate(pbmc.crc, pbmc.control, bVals, file.path(path.to.04.output, "bVals", "CRC_vs_control"), 0.05, "Beta", FALSE)
perform_DMR_with_DMRcate(pbmc.lung, pbmc.control, bVals, file.path(path.to.04.output, "bVals", "Lung_vs_control"), 0.05, "Beta", FALSE)
perform_DMR_with_DMRcate(pbmc.cancer, pbmc.control, bVals, file.path(path.to.04.output, "bVals", "CRC_Lung_vs_control"), 0.05, "Beta", FALSE)
perform_DMR_with_DMRcate(breast.cancer, breast.control, bVals, file.path(path.to.04.output, "bVals", "Breast_vs_BreastControl"), 1e-10, "Beta", TRUE)

perform_DMR_with_DMRcate(pbmc.crc, pbmc.control, log2(bVals), file.path(path.to.04.output, "bVals_log", "CRC_vs_control"), 0.05, "Beta", FALSE)
perform_DMR_with_DMRcate(pbmc.lung, pbmc.control, log2(bVals), file.path(path.to.04.output, "bVals_log", "Lung_vs_control"), 0.05, "Beta", FALSE)
perform_DMR_with_DMRcate(pbmc.cancer, pbmc.control, log2(bVals), file.path(path.to.04.output, "bVals_log", "CRC_Lung_vs_control"), 0.05, "Beta", FALSE)
perform_DMR_with_DMRcate(breast.cancer, breast.control, log2(bVals), file.path(path.to.04.output, "bVals_log", "Breast_vs_BreastControl"), 1e-10, "Beta", FALSE)

perform_DMR_with_DMRcate(pbmc.crc, pbmc.control, mVals, file.path(path.to.04.output, "mVals", "CRC_vs_control"), 0.05, "M", FALSE)
perform_DMR_with_DMRcate(pbmc.lung, pbmc.control, mVals, file.path(path.to.04.output, "mVals", "Lung_vs_control"), 0.05, "M", FALSE)
perform_DMR_with_DMRcate(pbmc.cancer, pbmc.control, mVals, file.path(path.to.04.output, "mVals", "CRC_Lung_vs_control"), 0.05, "M", FALSE)
perform_DMR_with_DMRcate(breast.cancer, breast.control, mVals, file.path(path.to.04.output, "mVals", "Breast_vs_BreastControl"), 1e-10, "M", TRUE)


