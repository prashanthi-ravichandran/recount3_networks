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
study_name <- inputArgs[1]
agg_level <- "GTEx"
agg_level_w <- "GTEx_minimal_corrected_single"
saveDir <- paste0("/data/abattle4/prashanthi/recount3_paper/results/WGCNA/networks/", agg_level_w, "/", study_name,"/")

expr <-  readRDS(paste0(datDir, agg_level, "/minimal_corrected_expression/", study_name, ".rds"))
common_genes <- readRDS(paste0(datDir, agg_level, "/common_genes.rds"))
expr <- expr[rownames(expr) %in% common_genes, ]
expr <- as.data.frame(t(expr)) 
pick.power = function(dat, nType = "unsigned"){
  pickSoftThreshold.output = pickSoftThreshold(dat, networkType = nType, RsquaredCut = 0.8)
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

dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average") 

for(i in c(1:length(cutheights))){
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

 