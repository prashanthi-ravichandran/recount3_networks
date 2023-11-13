# Description: This script reads in a particular network for distinct densities 
# and computes 
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

rm(list = ls())
library(Matrix)
library(cowplot)
library(mvtnorm)
library(igraph)
library(plyr)
library(RColorBrewer)
library(pheatmap)
library(scales)

inputArgs <- commandArgs(TRUE)
study <- inputArgs[1]
nStudies <- inputArgs[2]
lambda <- inputArgs[3]

home.dir <- "/data/abattle4/prashanthi/recount3_paper/"
dat.dir <- paste0(home.dir, "data/")
res.dir <- paste0(home.dir, "results/weighted_cov_networks/")

# Read in the network
net <- readRDS(paste0(res.dir, 
                      study, "/net_", nStudies, "/lambda_", lambda, ".rds"))

geneData <- readRDS(paste0(dat.dir, "geneData.rds"))

get_graph_distance <- function(inet){
  # This function takes in a precision matrix and converts
  # it to a weighted graph 
  inet <- as.matrix(inet)
  genes <- rownames(inet)
  diag(inet) <- 0.0
  inet <- abs(inet)
  inet <- inet/max(inet)
  # Get a distance 
  inet[!inet == 0.0] <- 1.0000000000000000001 - inet[!inet == 0.0]
  inet <- as(inet, "sparseMatrix")
  g <- graph_from_adjacency_matrix(inet, mode = c("undirected"), diag = F, weighted = T)
  V(g)$label <- genes
  g
}

get_graph_cov <- function(inet){
  # This function takes in a precision matrix and converts
  # it to a weighted graph 
  inet <- as.matrix(inet)
  genes <- rownames(inet)
  diag(inet) <- 0
  inet <- abs(inet)
  g <- graph_from_adjacency_matrix(inet, mode = c("undirected"), diag = F, weighted = T)
  V(g)$label <- genes
  g
}

graph_distance <- get_graph_distance(net)
graph_cov <- get_graph_cov(net)
centrality.df <- data.frame(matrix(ncol = 8, nrow = length(V(graph_distance))))
colnames(centrality.df) <- c("node", "degree", "maximum_weight", "strength", 
                             "betweeness", "closeness", "eigen_centrality", "page_rank")
centrality.df$node <- V(graph_distance)$label
centrality.df$degree <- degree(graph_cov, v = V(graph_cov))
centrality.df$maximum_weight <- ldply(V(graph_cov), function(v) range(incident(graph_cov, v, mode='total')$weight))$V2
centrality.df$maximum_weight[centrality.df$maximum_weight < 0] <- 0
centrality.df$strength <- strength(graph_cov, v = V(graph_cov))
centrality.df$betweeness <- betweenness(graph_distance, v = V(graph_distance))
sample.dist.matrix <- distances(graph_distance, v = V(graph_distance), to = V(graph_distance), algorithm = "dijkstra")
sample.closeness.matrix <- 1/sample.dist.matrix
diag(sample.closeness.matrix) <- 0
centrality.df$closeness <- rowSums(sample.closeness.matrix)
centrality.df$eigen_centrality <- eigen_centrality(graph_cov)$vector
centrality.df$page_rank <- page_rank(graph_cov, v = V(graph_cov))$vector

fnorm <- function(a){
  (a - min(a))/ (max(a) - min(a))
}

centrality.df$degree <- fnorm(centrality.df$degree)
centrality.df$maximum_weight <- fnorm(centrality.df$maximum_weight)
centrality.df$strength <- fnorm(centrality.df$strength)
centrality.df$betweeness <- fnorm(centrality.df$betweeness)
centrality.df$closeness <- fnorm(centrality.df$closeness)
centrality.df$eigen_centrality <- fnorm(centrality.df$eigen_centrality)
centrality.df$page_rank <- fnorm(centrality.df$page_rank)

saveRDS(centrality.df, paste0(res.dir, study, "/net_", 
                              nStudies, "/centrality_", lambda, ".rds"))


# Read in the gene coordinate information
gene_loc <- readRDS(dat.dir, "gene_loc.rds")
gene_loc <- gene_loc[gene_loc$gene_id %in% centrality.df$node, ]
gene_loc <- gene_loc[match(centrality.df$node, gene_loc$gene_id), ]
window <- 100000
gene_loc$start <- gene_loc$start - window
gene_loc$end <- gene_loc$end + window
gene_loc$start[gene_loc$start  < 0] <- 0

centrality.df <- cbind(gene_loc[ ,1:3], centrality.df)

if(study == "blood/all"){
  for(i in c(5:11)){
    write.table(centrality.df[ ,c(1,2,3,i)], 
                paste0(home.dir, "results/s_LDSC/blood_consensus_", lambda,"/",
                       colnames(centrality.df)[i], "/bed/", 
                       "blood_consensus_", lambda, ".bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
}else{
  if(study == "central_nervous_system/all"){
    for(i in c(5:11)){
      write.table(centrality.df[ ,c(1,2,3,i)], 
                  paste0(home.dir, "results/s_LDSC/CNS_consensus_", lambda,"/",
                         colnames(centrality.df)[i], "/bed/", 
                         "CNS_consensus_", lambda, ".bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)}
  }else{
    if(study == "blood/GTEx"){
      for(i in c(5:11)){
        write.table(centrality.df[ ,c(1,2,3,i)], 
                    paste0(home.dir, "results/s_LDSC/blood_GTEx_", lambda,"/",
                           colnames(centrality.df)[i], "/bed/", 
                           "blood_GTEx_", lambda, ".bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)}
    }else{
      if(study == "central_nervous_system/GTEx"){
        for(i in c(5:11)){
          write.table(centrality.df[ ,c(1,2,3,i)], 
                      paste0(home.dir, "results/s_LDSC/CNS_GTEx_", lambda,"/",
                             colnames(centrality.df)[i], "/bed/", 
                             "CNS_GTEx_", lambda, ".bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)}
      }else{  
        for(i in c(5:11)){
          write.table(centrality.df[ ,c(1,2,3,i)], 
                      paste0(home.dir, "results/s_LDSC/", study, "_", lambda,"/",
                             colnames(centrality.df)[i], "/bed/", 
                             study, "_", lambda, ".bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)}
      }
    }}}


