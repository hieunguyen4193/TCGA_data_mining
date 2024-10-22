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

BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2manifest", update = FALSE)
install.packages("/media/hieunguyen/GSHD_HN01/storage/offline_pkgs/IlluminaHumanMethylationEPICv2anno.20a1.hg38_1.0.0.tar.gz",
                 type = "sources", repos = NULL)
BiocManager::install("methylationArrayAnalysis", update = FALSE)
install.packages("nnls")

library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")

maindir <- "/media/hieunguyen/HNSD01"
# maindir <- "/media/hieunguyen/GSHD_HN01"

outdir <- file.path(maindir, "outdir")
PROJECT <- "EPIC"
output.version <- "20241020"
data.version <- "20241007"
path.to.storage <- "/media/hieunguyen/GSHD_HN01/storage"
path.to.main.input <- file.path(path.to.storage, PROJECT, data.version)
path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
path.to.02.output <- file.path(path.to.main.output, "EPIC_02_output")

resdf <- readxl::read_excel(file.path(path.to.02.output, "deconvolution_results.xlsx"))

input.vals <- list(
  bVals = readRDS(file.path(path.to.02.output, "bVals.rds")),
  mVals = readRDS(file.path(path.to.02.output, "mVals.rds"))
  )

input.vals[["bVals_log"]] <- log2(readRDS(file.path(path.to.02.output, "bVals.rds")))

for (input.type in names(input.vals)){
  path.to.03.output <- file.path(path.to.main.output, "EPIC_03_output", input.type)
  dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
  
  tissue.types <- c("Liver", "Breast", "Lung", "Gastric", "CRC", "WBC")
  
  pbmc.control <- subset(resdf, resdf$SampleType == "PBMC" & resdf$Label == "Control")$SampleCode
  pbmc.cancer <- subset(resdf, resdf$SampleType == "PBMC" & resdf$Label != "Control")$SampleCode
  pbmc.crc <- subset(resdf, resdf$SampleType == "PBMC" & resdf$Label == "CRC-cancer")$SampleCode
  pbmc.lung <- subset(resdf, resdf$SampleType == "PBMC" & resdf$Label == "Lung-cancer")$SampleCode
  
  breast.cancer <- subset(resdf, resdf$SampleType == "Tissue" & resdf$Label == "Breast-cancer")$SampleCode
  breast.control <- subset(resdf, resdf$SampleType == "Tissue" & resdf$Label != "Breast-cancer")$SampleCode
  run_test <- function(group1, group2, inputdf){
    input.metadata <- data.frame(sample = c(group1, group2), 
                                 label = c(
                                   to_vec(for(item in seq(1, length(group1))) "cancer"),
                                   to_vec(for(item in seq(1, length(group2))) "control")
                                 ))
    # this is the factor of interest
    # positive logFC  --> control > cancer
    # negative logFC --> cancer > control
    g <- factor(input.metadata$label, levels = c("control", "cancer"))
    # use the above to create a design matrix
    design <- model.matrix(~0+label, data=input.metadata)
    colnames(design) <- levels(g)
    fit <- lmFit(inputdf[, input.metadata$sample], design)
    # create a contrast matrix for specific comparisons
    contMatrix <- makeContrasts(cancer-control,
                                levels=design)
    fit2 <- contrasts.fit(fit, contMatrix)
    fit2 <- eBayes(fit2)
    
    DMPs <- topTable(fit2, num=Inf, coef=1) %>% as.data.frame() %>% rownames_to_column("probe")
    return(DMPs)  
  }
  
  annEPIC <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38) %>% as.data.frame() %>% rownames_to_column("probe")
  
  #####----------------------------------------------------------------------#####
  ##### RUN DML
  #####----------------------------------------------------------------------#####
  test_cancer_vs_control <- run_test(pbmc.cancer, pbmc.control, input.vals[[input.type]])
  test_cancer_vs_control <- merge(test_cancer_vs_control, annEPIC, by.x = "probe", by.y = "probe")
  test_cancer_vs_control.raw <- test_cancer_vs_control
  test_cancer_vs_control <- subset(test_cancer_vs_control, test_cancer_vs_control$adj.P.Val <= 0.05)
  writexl::write_xlsx(test_cancer_vs_control, file.path(path.to.03.output, "test_cancer_vs_control.xlsx"))
  writexl::write_xlsx(test_cancer_vs_control.raw, file.path(path.to.03.output, "test_cancer_vs_control.raw.xlsx"))
  
  test_crc_vs_control <- run_test(pbmc.crc, pbmc.control, input.vals[[input.type]])
  test_crc_vs_control <- merge(test_crc_vs_control, annEPIC, by.x = "probe", by.y = "probe")
  test_crc_vs_control.raw <- test_crc_vs_control 
  test_crc_vs_control <- subset(test_crc_vs_control, test_crc_vs_control$adj.P.Val <= 0.05)
  writexl::write_xlsx(test_crc_vs_control.raw, file.path(path.to.03.output, "test_crc_vs_control.raw.xlsx"))
  writexl::write_xlsx(test_crc_vs_control, file.path(path.to.03.output, "test_crc_vs_control.xlsx"))
  
  test_lung_vs_control <- run_test(pbmc.lung, pbmc.control, input.vals[[input.type]])
  test_lung_vs_control <- merge(test_lung_vs_control, annEPIC, by.x = "probe", by.y = "probe")
  test_lung_vs_control.raw <- test_lung_vs_control
  test_lung_vs_control <- subset(test_lung_vs_control, test_lung_vs_control$adj.P.Val <= 0.05)
  writexl::write_xlsx(test_lung_vs_control, file.path(path.to.03.output, "test_lung_vs_control.xlsx"))
  writexl::write_xlsx(test_lung_vs_control.raw, file.path(path.to.03.output, "test_lung_vs_control.raw.xlsx"))
  
  test_breast_vs_control <- run_test(breast.cancer, breast.control, input.vals[[input.type]])
  test_breast_vs_control <- merge(test_breast_vs_control, annEPIC, by.x = "probe", by.y = "probe")
  test_breast_vs_control.raw <- test_breast_vs_control
  test_breast_vs_control <- subset(test_breast_vs_control, test_breast_vs_control$adj.P.Val <= 0.05)
  
  test_breast_vs_control <- test_breast_vs_control %>% rowwise() %>%
    mutate(UCSC_RefGene_Accession = paste(unique(str_split(UCSC_RefGene_Accession, ";")[[1]]), collapse = ",") ) %>%
    mutate(UCSC_RefGene_Name = paste(unique(str_split(UCSC_RefGene_Name, ";")[[1]]), collapse = ",") ) %>%
    mutate(UCSC_RefGene_Group = paste(unique(str_split(UCSC_RefGene_Group, ";")[[1]]), collapse = ",") ) %>%
    mutate(GencodeV41_Group = paste(unique(str_split(GencodeV41_Group, ";")[[1]]), collapse = ",") ) %>%
    mutate(GencodeV41_Name = paste(unique(str_split(GencodeV41_Name, ";")[[1]]), collapse = ",") ) %>%
    mutate(GencodeV41_Accession = paste(unique(str_split(GencodeV41_Accession, ";")[[1]]), collapse = ",") )
  
  writexl::write_xlsx(test_breast_vs_control, file.path(path.to.03.output, "test_breast_vs_control.xlsx"))
  writexl::write_xlsx(test_breast_vs_control.raw, file.path(path.to.03.output, "test_breast_vs_control.raw.xlsx"))
  
  #####----------------------------------------------------------------------#####
  ##### ANNOTATION AND PATHWAY ANALYSIS
  #####----------------------------------------------------------------------#####
  if ("clusterProfiler" %in% installed.packages() == FALSE){
    BiocManager::install("clusterProfiler")
  }
  library("clusterProfiler")
  
  breats.de.genes <- table(test_breast_vs_control$UCSC_RefGene_Name) %>% sort()
  breats.de.genes <- breats.de.genes[names(breats.de.genes) != ""] %>% as.data.frame() %>%
    separate_rows(Var1)
  
  sig.genes <- subset(breats.de.genes, breats.de.genes$Freq >= 10)$Var1
  convertdf <- bitr(breats.de.genes$Var1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  go.res <- enrichGO(universe = convertdf$ENTREZID,
                     gene = subset(convertdf, convertdf$SYMBOL %in% sig.genes)$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
  kegg.res <- enrichKEGG(gene = subset(convertdf, convertdf$SYMBOL %in% sig.genes)$ENTREZID,
                         organism = 'hsa',
                         pvalueCutoff = 0.05)
  
  write.csv(data.frame(kegg.res), file.path(path.to.03.output, "KEGG_res_Breast_vs_Control.csv"))
  write.csv(data.frame(go.res), file.path(path.to.03.output, "GO_res_Breast_vs_Control.csv"))
  
  dotplot.kegg <- dotplot(kegg.res, showCategory=30) + ggtitle("dotplot for ORA")
  dotplot.go <- dotplot(go.res, showCategory=30) + ggtitle("dotplot for ORA")
  ggsave(plot = dotplot.kegg, filename = "dotplot_KEGG_DML_genes.Breast.svg", path = path.to.03.output, device = "svg", width = 14, height = 10)
  ggsave(plot = dotplot.go, filename = "dotplot_GO_DML_genes.Breast.svg", path = path.to.03.output, device = "svg", width = 14, height = 10)
}
