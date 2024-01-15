# This script read in a set of genes and an expression matrix, which it subsets
# The expression matrix is then quantile normalized and covariances are calculated
# Input 1 (GTEx/ recount): Specifies the study 
# Input 2: The name of the expression file 
# Input 3: The directory the covariances should be saved to

rm(list = ls())
library(coop)
library(parallel)

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
datDir <- paste0(homeDir, "data/")
inputArgs <- commandArgs(TRUE)
agg_level <- inputArgs[1]
corrDir <- paste0(datDir, agg_level, "/corrected_expression/")
covDir <- paste0(datDir, agg_level, "/single_covariances/")

projects <- list.files(corrDir)
expr_list <- lapply(projects, function(iproj){
    cat(iproj, sep = "\n")
    readRDS(paste0(corrDir, iproj))
})
  
names(expr_list) <- gsub(".rds", "", projects)
  
nsamples <- lapply(expr_list, function(icount){
    dim(icount)[2]
})
nsamples <- unlist(nsamples)
min_sample_size <- 15
expr_list <- expr_list[nsamples >= min_sample_size]
  
common.genes <- rownames(expr_list[[1]])
for(i in c(2:length(expr_list))){
  print(length(common.genes))
  common.genes <- intersect(common.genes, rownames(expr_list[[i]]))
}
expr_list <- lapply(expr_list, function(iexpr){
  iexpr <- iexpr[rownames(iexpr) %in% common.genes, ]
  iexpr[match(common.genes, rownames(iexpr)), ]
})


nsamples <- lapply(expr_list, function(rse){
  dim(rse)[2]
})

nsamples <- unlist(nsamples)

for(i in c(1:length(expr_list))){
  cat(names(expr_list)[i], "\n")
  cat(dim(expr_list[[i]]), "\n")
  iexpr <- limma::normalizeQuantiles(expr_list[[i]])
  cov <- covar(scale(t(iexpr)))
  saveRDS(cov, paste0(covDir, names(expr_list)[i], ".rds"))
  rm(cov)
}
