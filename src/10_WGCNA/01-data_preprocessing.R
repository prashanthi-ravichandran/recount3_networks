rm(list = ls())
library_paths <- readRDS("/data/abattle4/prashanthi/recount3_paper/results/WGCNA/library_paths.rds")
.libPaths(library_paths)

library(WGCNA)

net <- "sra_normal"
nStudies <- 49
if(net == "GTEx"){
  stepSize <- 10
}
#nAgg <- unique(c(1, 3, 5, 7, seq(stepSize, nStudies, stepSize), nStudies))
nAgg <- c(20)
datDir <- "/data/abattle4/prashanthi/recount3_paper/data/"
resDir <- "/data/abattle4/prashanthi/recount3_paper/results/WGCNA/processed-data/sra_normal/"
studies <- list.files(paste0(datDir, net, "/corrected_expression"))
study_metadata <- readRDS(paste0(datDir, net, "/study_metaData.rds"))
study_metadata <- study_metadata[order(study_metadata$n), ]
# Read the expression values
expr <- lapply(study_metadata$study, function(istudy){
  cat(istudy, "\n")
  readRDS(paste0(datDir, net, "/corrected_expression/", istudy, ".rds"))
})
names(expr) <- study_metadata$study
common.genes <- rownames(expr[[1]])
for(i in c(1:length(expr))){
  common.genes <- intersect(common.genes, rownames(expr[[i]]))
}
expr <- lapply(expr, function(iexpr){
  iexpr <- iexpr[rownames(iexpr) %in% common.genes, ]
  iexpr[match(common.genes, rownames(iexpr)), ]
})
for(i in nAgg){
  select_studies <- study_metadata$study[1:i]
  select_expr <- expr[names(expr) %in% select_studies]
  nSets <- length(select_expr)
  setLabels = select_studies
  shortLabels = select_studies
  multiExpr = vector(mode = "list", length = nSets)
  for(i in c(1:length(select_expr))){
    multiExpr[[i]] <- list(data = as.data.frame(t(select_expr[[i]])))
  }
  exprSize = checkSets(multiExpr)
  gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
  gsg$allOK
  if (!gsg$allOK)
  {
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
      printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                                collapse = ", ")))
    for (set in 1:exprSize$nSets){
      if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",
                         paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
      # Remove the offending genes and samples
      multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
    # Update exprSize
    exprSize = checkSets(multiExpr)
  }
  # Define data set dimensions
  nGenes = exprSize$nGenes
  nSamples = exprSize$nSamples
  
  save(multiExpr, nGenes, nSamples, setLabels, shortLabels, exprSize,
       file = paste0(resDir, "dat_", i ,".RData"))
}



