check_methylation_array_type <- function(input.files){
  library("illuminaio")
  library("BiocParallel")
  library("DelayedArray")
  library(dplyr)
  library(tidyverse)
  
  BACKEND <- getAutoRealizationBackend()
  BPREDO <- list()
  BPPARAM <- SerialParam()
  recursive <- FALSE
  extended <- FALSE
  force <- TRUE
  verbose <- TRUE
  
  G.Quants <- bplapply(input.files, function(xx) {
    if (verbose) 
      message("[read.metharray] Reading ", basename(xx))
    Quants <- readIDAT(xx)[["Quants"]]
    if (!extended) {
      Quants <- Quants[, "Mean", drop = FALSE]
    }
    if (!is.null(BACKEND)) {
      Quants <- realize(Quants, BACKEND = BACKEND)
    }
    Quants
  }, BPREDO = BPREDO, BPPARAM = BPPARAM)
  
  allNProbes <- vapply(G.Quants, nrow, integer(1L))
  
  guessArrayTypes <- function (nProbes) 
  {
    if (nProbes >= 622000 && nProbes <= 623000) {
      arrayAnnotation <- c(array = "IlluminaHumanMethylation450k")
    }
    else if (nProbes >= 1050000 && nProbes <= 1053000) {
      arrayAnnotation <- c(array = "IlluminaHumanMethylationEPIC")
    }
    else if (nProbes >= 1032000 && nProbes <= 1033000) {
      arrayAnnotation <- c(array = "IlluminaHumanMethylationEPIC")
    }
    else if (nProbes >= 55200 && nProbes <= 55400) {
      arrayAnnotation <- c(array = "IlluminaHumanMethylation27k")
    }
    else if (nProbes >= 54700 && nProbes <= 54800) {
      arrayAnnotation <- c(array = "IlluminaHumanMethylationAllergy")
    }
    else if (nProbes >= 41000 & nProbes <= 41100) {
      arrayAnnotation <- c(array = "HorvathMammalMethylChip40")
    }
    else if (nProbes >= 43650 & nProbes <= 43680) {
      arrayAnnotation <- c(array = "IlluminaHumanMethylationAllergy")
    }
    else {
      arrayAnnotation <- c(array = "Unknown", annotation = "Unknown")
    }
    arrayAnnotation
  }
  
  arrayTypes <- cbind(do.call(rbind, lapply(allNProbes, guessArrayTypes)), 
                      size = allNProbes) %>% as.data.frame()
  
  arrayTypedf <- data.frame(path = input.files, arrayType = arrayTypes$array)
  return(arrayTypedf)
}

