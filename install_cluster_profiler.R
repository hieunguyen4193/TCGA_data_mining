BiocManager::install("msigdbr", update = FALSE)
BiocManager::install("org.Hs.eg.db", update = FALSE)
remove.packages("clusterProfiler")
remove.packages("DOSE")
remove.packages("GOSemSim")
path.to.install.dir <- "/media/hieunguyen/HNHD01/offline_pkgs/clusterProfiler"
BiocManager::install("fgsea", update = FALSE)
install.packages("tidygraph")
BiocManager::install("enrichplot", update = FALSE)
for (pkg in c("HDO.db_0.99.1.tar.gz",
              "yulab.utils_0.1.4.tar.gz",
              "GO.db_3.18.0.tar.gz",
              "GOSemSim_2.28.1.tar.gz",
              "DOSE_3.28.2.tar.gz",
              "gson_0.1.0.tar.gz",
              "downloader_0.4.tar.gz",
              "clusterProfiler_4.10.1.tar.gz")){
  install.packages(file.path(path.to.install.dir, pkg), type = "source", repos = NULL)
}
