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
cutheights = seq(0.9,1.0,length.out = 50)

datDir <- "/data/abattle4/prashanthi/recount3_paper/data/"
resDir <- "/data/abattle4/prashanthi/recount3_paper/WGCNA/"

inputArgs <- commandArgs(TRUE)
nStudies <- inputArgs[1]
start_index <- as.numeric(inputArgs[2])
#agg_level <- inputArgs[2]
# agg_level_w <- inputArgs[3]
#nStudies <- "1"
agg_level <- "GTEx"
agg_level_w <- "GTEx"
saveDir <- paste0("/data/abattle4/prashanthi/recount3_paper/results/WGCNA/networks/", agg_level_w, "/net_", nStudies, "/")

study_metaData <- readRDS(paste0(datDir, agg_level, "/study_metaData.rds"))
study_metaData <- study_metaData[order(study_metaData$n), ]
common_genes <- readRDS(paste0(datDir, agg_level, "/common_genes.rds"))
study_metaData <- study_metaData[1:nStudies, ]

expr <- lapply(study_metaData$study, function(istudy){
  cat(istudy, "\n")
  readRDS(paste0(datDir, agg_level, "/corrected_expression/", istudy, ".rds"))
})
expr <- lapply(expr, function(iexpr){
  iexpr <- iexpr[rownames(iexpr) %in% common_genes, ]
  iexpr <- iexpr[match(common_genes, rownames(iexpr)), ]
  as.data.frame(t(iexpr)) 
})

expr <- do.call(rbind, expr)

TOM_computed <- TRUE
if(TOM_computed){
  TOM <- readRDS(paste0(saveDir, "TOM_similarity.rds"))
}else{
  pick.power = function(dat, nType = "unsigned"){
    pickSoftThreshold.output = pickSoftThreshold(dat, networkType = nType, RsquaredCut = 0.7)
    power <- pickSoftThreshold.output$powerEstimate
    if(is.na(power)){
      print(paste("no power reached r-suared cut-off, now choosing max r-squared based power"))
      power <- pickSoftThreshold.output$fitIndices$Power[which(pickSoftThreshold.output$fitIndices$SFT.R.sq == max(pickSoftThreshold.output$fitIndices$SFT.R.sq))]
    }
    power
  }
  soft_threshold <- pick.power(expr, nType = type.tom)
  
  adjacency = adjacency(expr, power = soft_threshold)
  # Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency)
  saveRDS(TOM, paste0(saveDir, "TOM_similarity.rds"))
}


dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average") 

for(i in c(start_index:length(cutheights))){
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize, cutHeight = cutheights[i])
  dynamicColors <- labels2colors(dynamicMods)
  moduleColors <- list()
  mergedMEs <- list()
  if(length(unique(dynamicColors)) > 1){
    MEList = moduleEigengenes(expr, colors = dynamicColors)
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs)
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average")
    merge = mergeCloseModules(expr, dynamicColors, cutHeight = mergeCutHeight, verbose = 3)
    moduleColors <- merge$colors
    mergedMEs <- merge$newMEs
  }else{
    moduleColors <- dynamicColors
    MEList = moduleEigengenes(expr, colors = dynamicColors)
    MEs = MEList$eigengenes
    mergedMEs <- MEs
  }
  colorOrder = c("grey", standardColors(50))
  moduleLabels <- match(moduleColors, colorOrder)-1
  saveRDS(list("cutheight" = cutheights[i],
               "genes" = colnames(expr),
               "mergedMEs" = mergedMEs, 
               "moduleLabels" = moduleLabels,
               "moduleColors" = moduleColors), file = paste0(saveDir, "net_", i, ".rds"))
}