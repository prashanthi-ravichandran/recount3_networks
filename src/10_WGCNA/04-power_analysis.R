rm(list = ls())
library_paths <- readRDS("/data/abattle4/prashanthi/recount3_paper/results/WGCNA/library_paths.rds")
.libPaths(library_paths)
library(WGCNA)
library(RColorBrewer)

type.tom = "unsigned"
minModuleSize = 10
reassignThreshold = 0
mergeCutHeight = 0.20
numericLabels = TRUE
pamRespectsDendro = FALSE
verbose = 3
cutheights = seq(0.95,1.0,length.out = 50)

datDir <- "/data/abattle4/prashanthi/recount3_paper/data/"
resDir <- "/data/abattle4/prashanthi/recount3_paper/WGCNA/"

inputArgs <- commandArgs(TRUE)
nStudies <- inputArgs[1]
agg_level <- "GTEx"
agg_level_w <- "GTEx_minimal_corrected"
threshold <- as.numeric(inputArgs[2])

saveDir <- paste0("/data/abattle4/prashanthi/recount3_paper/results/WGCNA/networks/", agg_level_w, "/net_", nStudies, "/", "power_analysis/")

study_metaData <- readRDS(paste0(datDir, agg_level, "/study_metaData.rds"))
study_metaData <- study_metaData[order(study_metaData$n), ]
common_genes <- readRDS(paste0(datDir, agg_level, "/common_genes.rds"))
study_metaData <- study_metaData[1:nStudies, ]

expr <- lapply(study_metaData$study, function(istudy){
  cat(istudy, "\n")
  readRDS(paste0(datDir, agg_level, "/minimal_corrected_expression/", istudy, ".rds"))
})
expr <- lapply(expr, function(iexpr){
  iexpr <- iexpr[rownames(iexpr) %in% common_genes, ]
  iexpr <- iexpr[match(common_genes, rownames(iexpr)), ]
  as.data.frame(t(iexpr)) 
})

expr <- do.call(rbind, expr)

adjacency <- adjacency(expr, power = threshold)
TOM <- TOMsimilarity(adjacency)

saveRDS(TOM, paste0(saveDir, "power_", threshold, "/TOM.rds"))

dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average") 

dynamicMods <- lapply(cutheights, function(c){
 cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize, cutHeight = c)
})

nUnassigned_genes <- lapply(dynamicMods, function(dynamicMod){
  sum(dynamicMod == 0)
})
fraction_unassigned <- unlist(nUnassigned_genes)/length(common_genes)

res_df <- data.frame("fraction_unassigned_genes" = fraction_unassigned, 
                     "cutHeights" = cutheights)

saveRDS(list("cutHeights" = cutheights, "dynamicMods" = dynamicMods, "summary" = res_df), 
        paste0(saveDir, "power_", threshold, "/results.rds"))

