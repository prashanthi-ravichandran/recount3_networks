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

inputArgs <- commandArgs(TRUE)
net <- inputArgs[1]
nStudies <- inputArgs[2]
# net <- "sra_normal"
# nStudies <- "300"
load(paste0("/data/abattle4/prashanthi/recount3_paper/results/WGCNA/processed-data/", net, "/dat_", nStudies , ".RData"))
nSets <- checkSets(multiExpr)$nSets
saveDir <- paste0("/data/abattle4/prashanthi/recount3_paper/results/WGCNA/networks/", net, "/consensus/net_", nStudies, "/")
TOMDir <- paste0("/data/abattle4/prashanthi/recount3_paper/results/WGCNA/networks/", net, "/consensus/TOMs/")
cutheights = seq(0.98, 0.997, length.out = 50)

pick.power = function(dat, nType = "unsigned"){
  pickSoftThreshold.output = pickSoftThreshold(dat, networkType = nType, RsquaredCut = 0.8)
  power <- pickSoftThreshold.output$powerEstimate
  if(is.na(power)){
    print(paste("no power reached r-suared cut-off, now choosing max r-squared based power"))
    power <- pickSoftThreshold.output$fitIndices$Power[which(pickSoftThreshold.output$fitIndices$SFT.R.sq == max(pickSoftThreshold.output$fitIndices$SFT.R.sq))]
  }
  power
}
type.tom = "unsigned"

# for(set in 1:nSets){
#   cat(set, "\n")
#   soft_threshold <- pick.power(multiExpr[[set]]$data, nType = type.tom)
#   adjacency_i <- abs(cor(multiExpr[[set]]$data, use = "p"))^soft_threshold
#   TOM <- TOMsimilarity(adjacency_i)
#   saveRDS(TOM, paste0(TOMDir, "TOM_", set, ".rds"))
# }

TOM = array(0, dim = c(nSets, nGenes, nGenes))
for (set in 1:nSets){
  cat(set, "\n")
  TOM[set, , ] = readRDS(paste0(TOMDir, "TOM_", set, ".rds"))
}

# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets){
  cat(set, "\n")
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
   # Scale the male TOM
   if (set>1)
   {
     scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
     TOM[set, ,] = TOM[set, ,]^scalePowers[set]
   }
 }


consensusTOM = TOM[1, ,]
if(nSets > 1){
 for(set in c(2:nSets)){
     consensusTOM = consensusTOM + TOM[set, ,]
}
}
consensusTOM <- consensusTOM/nSets

saveRDS(consensusTOM, paste0(saveDir, "consensus_TOM.rds"))

 # Clustering
for(i in c(1:length(cutheights))){
   cat(i, "\n")
   consTree = hclust(as.dist(1-consensusTOM), method = "average")
   unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                                  deepSplit = 2, cutHeight = cutheights[i],
                                  minClusterSize = minModuleSize,
                                  pamRespectsDendro = FALSE )
   unmergedColors = labels2colors(unmergedLabels)
   # Calculate module eigengenes
   unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
   # Calculate consensus dissimilarity of consensus module eigengenes
   consMEDiss = consensusMEDissimilarity(unmergedMEs)
   # Cluster consensus modules
   consMETree = hclust(as.dist(consMEDiss), method = "average")
   merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)
   # Numeric module labels
   moduleLabels = merge$colors
   # Convert labels to colors
   moduleColors = labels2colors(moduleLabels)
   # Eigengenes of the new merged modules:
   consMEs = merge$newMEs
   saveRDS(list("cutheight" = cutheights[i],
               "genes" = colnames(multiExpr[[1]]$data),
               "mergedMEs" = consMEs,
               "moduleLabels" = moduleLabels,
               "moduleColors" = moduleColors), file = paste0(saveDir, "net_", i, ".rds"))
}

