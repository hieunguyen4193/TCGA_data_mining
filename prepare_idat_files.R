gc()
library(dplyr)
library(tidyverse)

inputdir <- "/media/hieunguyen/GSHD_HN01/idat_TCGA"
outdir <- "/media/hieunguyen/HNSD01/outdir"

output.version <- "20240612"
path.to.main.input <- file.path(outdir, sprintf("TCGA_%s", output.version), sprintf("idat_%s", output.version))
path.to.main.output <- file.path(outdir, sprintf("TCGA_%s", output.version))

dir.create(path.to.main.input, showWarnings = FALSE, recursive = TRUE)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(inputdir, "metadata", output.version), showWarnings = FALSE, recursive = TRUE)

raw.metadata <- read.table(file.path(inputdir, "metadata", "gdc_sample_sheet.2024-06-11.tsv"), sep = "\t", header = TRUE) 
meta.data <- raw.metadata %>%
  subset(Data.Category == "DNA Methylation" & grepl(".idat", File.Name) == TRUE) %>%
  subset(Sample.Type %in% c("Primary Tumor", "Solid Tissue Normal")) %>% arrange(desc(Case.ID))

write.csv(meta.data, file.path(inputdir, "metadata", output.version, "gdc_sample_sheet.2024-06-11.idat_only.primary_tumor_normal_only.tsv"))

##### Move data to main input dir
for (i in seq(1, nrow(meta.data))){
  file.id <- meta.data[i, ]$File.ID
  File.name <- meta.data[i, ]$File.Name
  system(sprintf("cp %s %s", 
                 file.path(inputdir, "idat_files", file.id, File.name),
                 file.path(path.to.main.input)))
}
