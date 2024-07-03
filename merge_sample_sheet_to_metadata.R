gc()
rm(list = ls())

library(dplyr)
library(tidyverse)

infodir <- "/media/hieunguyen/HNSD01/src/TCGA_data_analysis/TCGA_database_DMR"
outdir <- "/media/hieunguyen/HNSD01/outdir"
raw.data.dir <- "/media/hieunguyen/HNHD01/TCGA_all_idat"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/TCGA_data_analysis"

source(file.path(path.to.project.src, "check_methylation_array_type.R"))

PROJECT <- "TCGA_methyl_panel"
output.version <- "20240620"
data.version <- "raw"

path.to.main.output <- file.path(outdir, PROJECT, sprintf("data_%s", data.version), output.version)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.pool.data <- file.path(path.to.main.output, sprintf("pool_%s", output.version))
dir.create(path.to.pool.data, showWarnings = FALSE, recursive = TRUE)

path.to.metadata.dir <- file.path(path.to.main.output, "metadata")
dir.create(path.to.metadata.dir, showWarnings = FALSE, recursive = TRUE)

all.sample.sheets <- Sys.glob(file.path(infodir, "*idat", "*sample_sheet*450K*"))
meta.data <- data.frame()
for (f in all.sample.sheets){
  filename <- basename(f)
  if (grepl("Tumor", filename)){
    sample.type <- "tumor"
  } else {
    sample.type <- "normal"
  }
  label <- dirname(f) %>% basename() %>% str_replace("_idat", "")
  tmpdf <- read.csv(f, sep = "\t")
  tmpdf$tumor_or_normal <- sample.type
  tmpdf$label <- label
  tmpdf$dirname <- dirname(f) %>% basename()
  tmpdf <- tmpdf %>% rowwise() %>%
    mutate(path = file.path(raw.data.dir, data.version, label, sample.type, File.ID, File.Name)) %>%
    mutate(check.file.exists = file.exists(path)) %>%
    mutate(color = str_replace(str_split(File.Name, "_noid_")[[1]][[2]], ".idat", ""))
  meta.data <- rbind(meta.data, tmpdf)
}

meta.data[[sprintf("pool_%s", output.version)]] <- unlist(
  lapply(meta.data$File.Name, function(x){
    return(file.path(path.to.main.output, sprintf("pool_%s", output.version), x))
  })
)

meta.data <- meta.data %>% rowwise() %>%
  mutate(base.path = file.path(path.to.pool.data, str_replace(str_replace(File.Name, "_Grn.idat", ""), "_Red.idat", "")))
write.csv(meta.data, file.path(path.to.metadata.dir, "collected_TCGA_metadata.csv"))

if (file.exists(file.path(path.to.metadata.dir, "check_array_type.csv")) == FALSE){
  check.array.type <- check_methylation_array_type(subset(meta.data, meta.data$color == "Grn")$path)
  write.csv(check.array.type, file.path(path.to.metadata.dir, "check_array_type.csv"))  
} else {
  check.array.type <- read.csv(file.path(path.to.metadata.dir, "check_array_type.csv"))
}

# for (i in seq(1, nrow(meta.data))){
#   orig.path <- meta.data[i, ]$path
#   dest.path <- path.to.pool.data
#   system(sprintf("cp %s %s", orig.path, dest.path))
# }

