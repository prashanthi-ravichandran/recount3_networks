rm(list = ls())
.libPaths(c("/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-markdown-1.1-g65guwxov2k2maedod2wedfr5jon6egu/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-knitr-1.28-cn7dhiz6mwl53op4gpecto35sljr4muz/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-yaml-2.2.0-ttbd4ipwa5jccbl733saf7c2zijmdixr/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-stringr-1.4.0-qbr2amu2xnxkicncfgvu7klzll4dg46v/rlib/R/library",
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-stringi-1.4.3-n22ruwgbot2i3becz3comesic75s47r6/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-magrittr-1.5-wy6q2ditqph62m4dib33mds3wsbluj7g/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-glue-1.4.1-5ejojwbzonzsvkmwa7o2h57gpbsutbl4/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-mime-0.7-n4hlpgd2g5fdbt3idl374eftb2utep3r/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-highr-0.8-rcw4hovy72yo4wj6mfrgtgbkovlqb3fg/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-evaluate-0.14-laasv6wsxntd3ht34ydl4ott7zzfthko/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-4.0.2-amdvcpog4ugspqwwx3ari7pzkmckelu6/rlib/R/library", 
            "/home/pravich2/rlibs/4.0.2/gcc/9.3.0"))
# Load required libraries
rm(list = ls())
library(Matrix)
library(matrixStats)
library(reshape2)
library(mvtnorm)
library(igraph)
library(org.Hs.eg.db)

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
net <- "sra_normal/consensus"
nStudies <- 566
res.dir <- paste0(homeDir, "results/WGCNA/networks/", net, "/net_", nStudies, "/")
dat.dir <- paste0(homeDir, "data/")

geneData <- readRDS(paste0(dat.dir, "geneData.rds"))

common_genes <- readRDS(paste0(dat.dir,"sra_normal/common_genes.rds"))
common_genes_symbol <- c()
for(i in c(1:length(common_genes))){
  common_genes_symbol[i] <- geneData$gene_name[geneData$gene_id == common_genes[i]]
}
# Read in the true edges
# canonical_edges <- readRDS(paste0(dat.dir, "canonical_pathways_edgelist.rds"))
# gene1 <- gsub("_.*", "", canonical_edges)
# gene2 <- gsub(".*_", "", canonical_edges)
# canonical_edges.df <- data.frame(gene1, gene2)
# canonical_edges.df <- canonical_edges.df[canonical_edges.df$gene1 %in% common_genes_symbol, ]
# canonical_edges.df <- canonical_edges.df[canonical_edges.df$gene2 %in% common_genes_symbol, ]
# true_edges <- paste(canonical_edges.df$gene1, canonical_edges.df$gene2, sep = "_")
# 
TF_targets_df <- read.delim("/data/abattle4/prashanthi/recount3_paper/data/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv")
TF_targets_df <- TF_targets_df[grepl("chromatin immunoprecipitation assay", TF_targets_df$Detection.method), ]
TF_targets_df <- TF_targets_df[TF_targets_df$Small.scale.evidence == "No", ]
TF_targets_df <- TF_targets_df[TF_targets_df$Name.TF %in% common_genes_symbol, ]
TF_targets_df <- TF_targets_df[TF_targets_df$Name.Target %in% common_genes_symbol, ]
TF_target_edges <- c()
for(i in c(1:dim(TF_targets_df)[1])){
  genes <- c(TF_targets_df$Name.TF[i], TF_targets_df$Name.Target[i])
  genes <- genes[order(genes)]
  TF_target_edges[i] <- paste0(genes[1], "_", genes[2])
}
# Read in network edges
WGCNA_files <- paste0("net_", 1:50, ".rds")
WGCNA_res <- lapply(WGCNA_files, function(i){
  cat(i, "\n")
  readRDS(paste0(res.dir, i))
})
cutHeights <- lapply(WGCNA_res, function(res){
  res$cutheight
})
cutHeights <- unlist(cutHeights)

modules <- lapply(WGCNA_res, function(res){
  res$moduleLabels
})

network_list <- list()
for(i in c(1:length(cutHeights))){
  network_list[[i]] <- data.frame("cutHeights" = cutHeights[[i]], 
                                  "modules" = modules[[i]], 
                                  "genes" = common_genes_symbol)
}

net_to_edgeList <- function(net){
  edges <- list()
  net <- net[!net$modules == 0, ]
  net <- net[!is.na(net$modules), ]
  index <- 1
  for(m in unique(net$modules)){
    cat(m, "\n")
    genes <- net$genes[net$modules == m]
    pathway.edge <- t(combn(sort(genes),2))
    pathway.edge <- paste(pathway.edge[,1], pathway.edge[,2], sep = "_")
    edges[[index]] <- pathway.edge
    index <- index + 1
  }
  edges <- unlist(edges)
  edges
}

edges_list <- lapply(network_list, function(inet){
  net_to_edgeList(inet)
})

nEdges <- lapply(edges_list, function(iedges){
  length(iedges)
})
nEdges <- unlist(nEdges)

true_edges <- TF_target_edges
# Compute summary metrics
true_positives <- lapply(edges_list, function(iedge){
  length(intersect(iedge, true_edges))
})
false_positives <- lapply(edges_list, function(iedge){
  length(iedge) - length(intersect(iedge, true_edges))
})
false_negatives <- lapply(edges_list, function(iedge){
  length(true_edges) - length(intersect(iedge, true_edges))
})

true_positives <- unlist(true_positives)
false_positives <- unlist(false_positives)
false_negatives <- unlist(false_negatives)

recall <- true_positives/(true_positives + false_negatives)
precision <- true_positives/(true_positives + false_positives)
FDR <- false_positives/(false_positives + true_positives)


res_df <- data.frame("cutHeights" = cutHeights, "nEdges" = nEdges, 
                     "true_positives" = true_positives, "false_positives" = false_positives, 
                     "false_negatives" = false_negatives, "precision" = precision, "recall" = recall, "FDR" = FDR)


saveRDS(res_df, paste0(res.dir, "fdr_results_TF.rds"))
#saveRDS(res_df, paste0(res.dir, study, "/fdr_results.rds"))


