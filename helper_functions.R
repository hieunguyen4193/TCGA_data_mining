run_test_group1_vs_group2 <- function(input.mat, group1, group2){
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
  fit <- lmFit(input.mat[, input.metadata$sample], design)
  # create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(group1-group2,
                              levels=design)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  
  ann450kSub <- ann450k[match(rownames(input.mat[, input.metadata$sample]), ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
  DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
  
  #####----------------------------------------------------------------------#####
  ##### RUN MAIN TEST BY LIMMA, DMR
  #####----------------------------------------------------------------------#####
  myAnnotation <- cpg.annotate(object = input.mat[, input.metadata$sample], datatype = "array", what = "Beta", 
                               analysis.type = "differential", design = design, 
                               contrasts = TRUE, cont.matrix = contMatrix, 
                               coef = "group1 - group2", arraytype = "450K")
  
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, min.cpgs = 5)
  results.ranges <- extractRanges(DMRs)
  DMRdf <- data.frame(results.ranges)
  output <- list(
    dmr = DMRdf,
    dmp = DMPs
  )
}  
