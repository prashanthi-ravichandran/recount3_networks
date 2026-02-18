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

# This script analyses networks
lambda <- seq(0.10, 0.50, 0.02)
lambda <- sprintf("%.2f", lambda)

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
res.dir <- paste0(homeDir, "results/")
dat.dir <- paste0(homeDir, "data/")
geneData <- readRDS(paste0(dat.dir, "geneData.rds"))

edge.list <- function(net){
  edge.unformatted <- lapply(net, function(inet){
    if(!is.null(dim(inet))){
      inet <- as.matrix(inet)
      ngenes <- dim(inet)[1]
      genes <- rownames(inet)
      genes.ordered <- genes[order(genes)]
      inet <- inet[genes.ordered, genes.ordered]
      inet[lower.tri(inet, diag = T)] <- NA
      inet[inet == 0] <- NA
      iedge <- reshape2::melt(inet, na.rm = T)
      iedge <- paste(iedge[,1], iedge[,2], sep = "_")
    }
    else{
      iedge <- NULL
    }
    iedge
  })
  edge.unformatted
}

#inputArgs <- commandArgs(TRUE)
study <- "lung/all"
nStudies <- "10"
agg_type <- "weighted_cov_networks"

index <- 1
net <- list()
for(ilambda in lambda){
  print(ilambda)
  inet <- readRDS(paste0(res.dir, agg_type, "/", study, "/net_", nStudies, "/lambda_", ilambda, ".rds"))
  if(is.list(inet)){
    net[[index]] <- inet[[1]]
  }else{
    net[[index]] <- inet
  }
  index <- index + 1
}
gene_ids <- rownames(net[[1]])
gene_symbols <- c()
for(i in c(1:length(gene_ids))){
  gene_symbols[i] <- geneData$gene_name[geneData$gene_id == gene_ids[i]]
}
net <- lapply(net, function(inet){
  rownames(inet) <- gene_symbols
  colnames(inet) <- gene_symbols
  inet
})
inferred.networks <- edge.list(net)
n.edges <- lapply(inferred.networks, function(iedge){
  length(iedge)
})
n.edges <- unlist(n.edges)

# Read in TF-target information
TF_targets_df <- read.delim("/data/abattle4/prashanthi/recount3_paper/data/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv")
TF_targets_df <- TF_targets_df[grepl("chromatin immunoprecipitation assay", TF_targets_df$Detection.method), ]
TF_targets_df <- TF_targets_df[TF_targets_df$Small.scale.evidence == "No", ]
TF_targets_df <- TF_targets_df[TF_targets_df$Name.TF %in% rownames(net[[1]]), ]
TF_targets_df <- TF_targets_df[TF_targets_df$Name.Target %in% rownames(net[[1]]), ]
TF_target_edges <- c()
for(i in c(1:dim(TF_targets_df)[1])){
  genes <- c(TF_targets_df$Name.TF[i], TF_targets_df$Name.Target[i])
  genes <- genes[order(genes)]
  TF_target_edges[i] <- paste0(genes[1], "_", genes[2])
}

# Create a network with all possible edges
all_genes <- unique(c(rownames(net[[1]]), TF_targets_df$Name.TF, TF_targets_df$Name.Target))
ngenes <- length(all_genes)
all.net <- matrix(1, nrow = ngenes, ncol = ngenes)
colnames(all.net) <- all_genes[order(all_genes)]
rownames(all.net) <- all_genes[order(all_genes)]
all.net[lower.tri(all.net, diag = T)] <- NA
all.edges <- reshape2::melt(all.net, na.rm = T)
all.edgelist <- paste(all.edges[ ,1], all.edges[ ,2], sep = "_")
print("All edges computed")

# Compute precision and recall
frac.in.true_edges <- lapply(inferred.networks, function(iedge){
  length(intersect(iedge, TF_target_edges))/ length(iedge)
})
frac.in.true_edges <- unlist(frac.in.true_edges)

true_positives <- lapply(inferred.networks, function(iedge){
  length(intersect(iedge, TF_target_edges))
})

false_positives <- lapply(inferred.networks, function(iedge){
  length(iedge) - length(intersect(iedge, TF_target_edges))
})

false_negatives <- lapply(inferred.networks, function(iedge){
  length(TF_target_edges) - length(intersect(iedge, TF_target_edges))
})

true_negatives <- lapply(inferred.networks, function(iedge){
  cat(length(iedge), "\n")
  length(intersect(all.edgelist[!all.edgelist %in% TF_target_edges],all.edgelist[!all.edgelist %in% iedge]))
})

true_positives <- unlist(true_positives)
false_positives <- unlist(false_positives)
false_negatives <- unlist(false_negatives)
true_negatives <- unlist(true_negatives)

recall <- true_positives/(true_positives + false_negatives)
precision <- true_positives/(true_positives + false_positives)
fpr <- false_positives/(false_positives + true_negatives)
res_df <- data.frame(lambda, n.edges, frac.in.true_edges, recall, precision, fpr)

saveRDS(res_df, paste0(res.dir, agg_type, "/", study, "/net_", nStudies, "/TF_target_enrichment.rds"))



