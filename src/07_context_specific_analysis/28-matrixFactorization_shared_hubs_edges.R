# Description: Perform matrix factorization to discover patterns of 
# shared co-regulation among related tissues. 
rm(list = ls())
.libPaths(c("/home/pravich2/rlibs/4.0.2/gcc/9.3.0", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-4.0.2-amdvcpog4ugspqwwx3ari7pzkmckelu6/rlib/R/library", 
            "/data/apps/extern/tidyverse/1.3.1"))
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(rstatix)
library(RcppML)
library(igraph)
library(stringr)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

contexts <- c("adipose", "airway", "B_cells", "blood", "breast", "cardiac", 
              "central_nervous_system", "colon", "esophagus", "eye", "fibroblasts", 
              "hescs", "intestine", "ipscs", "kidney", "liver", "lung", 
              "multipotent_cells", "myeloid_cells", "nervous_system", "pancreas", 
              "pbmcs_t_cells", "prostate", "skeletal_muscle", "skin", "stomach", "vascular")
nAll_studies <- c()
nAll_samples <- c()
nAll_median_sample_size <- c()
homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
datDir <- paste0(homeDir, "data/")
for(i in c(1:length(contexts))){
  study_metaData <- readRDS(paste0(datDir, contexts[i], "/all/study_metaData.rds"))
  nAll_studies[i] <- dim(study_metaData)[1]
  nAll_samples[i] <- sum(study_metaData$n)
  nAll_median_sample_size[i] <- median(study_metaData$n)
}
nGTEx_studies <- c()
nGTEx_samples <- c()
nGTEx_median_sample_size <- c()
for(i in c(1:length(contexts))){
  if(contexts[i] %in% c("airway", "eye", "hescs", "ipscs", "multipotent_cells", "myeloid_cells", "pbmcs_t_cells")){
    nGTEx_studies[i] <- 0
    nGTEx_samples[i] <- 0
    nGTEx_median_sample_size[i] <- 0
  }
  else{
    study_metaData <- readRDS(paste0(datDir, contexts[i], "/GTEx/study_metaData.rds"))
    nGTEx_studies[i] <- dim(study_metaData)[1]
    nGTEx_samples[i] <- sum(study_metaData$n)
    nGTEx_median_sample_size[i] <- median(study_metaData$n)}
}

all_tissue_df <- data.frame(contexts = contexts,
                            nStudies = nAll_studies, 
                            nSamples = nAll_samples, 
                            nMedian_ss = nAll_median_sample_size)
all_tissue_df$Samples_agg <- "Tissue consensus"

GTEx_tissue_df <- data.frame(contexts = contexts, 
                             nStudies = nGTEx_studies, 
                             nSamples = nGTEx_samples, 
                             nMedian_ss = nGTEx_median_sample_size)
GTEx_tissue_df$Samples_agg <- "GTEx"

tissue_df <- rbind(all_tissue_df, GTEx_tissue_df)
tissue_df$contexts <- gsub("_", " ", tissue_df$contexts)
tissue_df$contexts <- str_to_sentence(tissue_df$contexts)
tissue_df$contexts[tissue_df$contexts == "Hescs"] <- "hESCs"
tissue_df$contexts[tissue_df$contexts == "Ipscs"] <- "iPSCs"
tissue_df$contexts[tissue_df$contexts == "Pbmcs t cells"] <- "PBMCs/T cells"

# Pick lambda for each context 
nEdges <- 7000
lambda <- c()
nEdges_actual <- c()
slope <- c()
rsquared <- c()
resDir <- paste0(homeDir, "results/weighted_cov_networks/")
for(i in c(1:length(contexts))){
  sf_res <- readRDS(paste0(resDir, 
                           contexts[i], "/all/net_", all_tissue_df$nStudies[all_tissue_df$contexts == contexts[i]] ,"/sf_res.rds"))
  lambda[i] <- sf_res$lambda[which.min(abs(sf_res$nEdges - nEdges))]
  nEdges_actual[i] <- sf_res$nEdges[sf_res$lambda == lambda[i]]
  slope[i] <- sf_res$slope.log_k[sf_res$lambda == lambda[i]]
  rsquared[i] <- sf_res$rsquared[sf_res$lambda == lambda[i]]
}

# Read nets 
net <- list()
for(i in c(1:length(contexts))){
  net[[i]] <- readRDS(paste0(resDir, 
                             contexts[i], "/all/net_", all_tissue_df$nStudies[all_tissue_df$contexts == contexts[i]] ,
                             "/lambda_", sprintf("%.2f", lambda[i]), ".rds"))
}

names(net) <- contexts

all_consensus_sf <- readRDS(paste0(resDir, "all_consensus/net_966/sf_res.rds"))
normal_consensus_sf <- readRDS(paste0(resDir, "normal_consensus/net_629/sf_res.rds"))
lambda_all_consensus <-  all_consensus_sf$lambda[which.min(abs(all_consensus_sf$n_edges - nEdges))]
lambda_normal_consensus <- normal_consensus_sf$lambda[which.min(abs(normal_consensus_sf$n_edges - nEdges))]

all_net <- readRDS(paste0(resDir, "all_consensus/net_966/lambda_",
                          sprintf("%.2f", lambda_all_consensus), ".rds"))

normal_net <- readRDS(paste0(resDir, "normal_consensus/net_629/lambda_",
                             sprintf("%.2f", lambda_normal_consensus), ".rds"))

net <- c(net, list("all_consensus" = all_net, "normal_consensus" = normal_net))

geneData <- readRDS(paste0(datDir, "geneData.rds"))
net <- lapply(net, function(inet){
  gene_ids <- rownames(inet)
  geneSymbols <- c()
  for(i in c(1:length(gene_ids))){
    geneSymbols[i] <- geneData$gene_name[geneData$gene_id == gene_ids[i]]
  }
  rownames(inet) <- geneSymbols
  colnames(inet) <- geneSymbols
  inet
})

edgeList <- lapply(net, function(inet){
  ngenes <- dim(inet)[1]
  genes <- rownames(inet)
  genes.ordered <- genes[order(genes)]
  inet <- inet[genes.ordered, genes.ordered]
  inet <- as.matrix(inet)
  inet[lower.tri(inet, diag = T)] <- NA
  inet[inet == 0] <- NA
  iedge <- reshape2::melt(inet, na.rm = T)
  iedge <- paste(iedge[,1], iedge[,2], sep = "_")
})

degreeCentrality_list <- lapply(net, function(inet){
  inet <- as.matrix(inet)
  genes <- rownames(inet)
  inet[lower.tri(inet)] <- 0
  diag(inet) <- 0
  inet[inet!=0] <- 1
  inet_graph <- graph_from_adjacency_matrix(inet, mode = "undirected", diag = F)
  V(inet_graph)$label <- genes
  degree_df <- degree(inet_graph, v = V(inet_graph))
  degree_df <- data.frame(genes = names(degree_df), degree = as.numeric(degree_df))
  degree_df
})

hub_list <- lapply(degreeCentrality_list, function(degree_df){
  threshold <- quantile(degree_df$degree, 0.95)
  degree_df$genes[degree_df$degree >= threshold]
})

uniqueHubGenes <- unique(unlist(hub_list))
cat(length(uniqueHubGenes), "\n")
uniqueEdges <- unique(unlist(edgeList))
cat(length(uniqueEdges), "\n")

in_Network_edges <- lapply(edgeList, function(iedge){
  in_Net <- rep(0, length(uniqueEdges))
  in_Net[uniqueEdges %in% iedge] <- 1
  in_Net
})

in_Network_hubs <- lapply(hub_list, function(ihub){
  in_Net <- rep(0, length(uniqueHubGenes))
  in_Net[uniqueHubGenes %in% ihub] <- 1
  in_Net
})

in_Network_edges <- do.call(cbind, in_Network_edges)
rownames(in_Network_edges) <- uniqueEdges
in_Network_edges <- in_Network_edges[rowSums(in_Network_edges) > 1, ]
in_Network_edges <- as(in_Network_edges, "sparseMatrix")
colnames(in_Network_edges) <- str_to_sentence(colnames(in_Network_edges))
colnames(in_Network_edges) <- gsub("_", " ", colnames(in_Network_edges))
colnames(in_Network_edges)[colnames(in_Network_edges) == "Hescs"] <- "hESCs"
colnames(in_Network_edges)[colnames(in_Network_edges) == "Ipscs"] <- "iPSCs"
colnames(in_Network_edges)[colnames(in_Network_edges) == "Pbmcs t cells"] <- "PBMCs/T cells"

in_Network_hubs <- do.call(cbind, in_Network_hubs)
rownames(in_Network_hubs) <- uniqueHubGenes
in_Network_hubs <- in_Network_hubs[rowSums(in_Network_hubs) > 1, ]
in_Network_hubs <- as(in_Network_hubs, "sparseMatrix")
colnames(in_Network_hubs) <- str_to_sentence(colnames(in_Network_hubs))
colnames(in_Network_hubs) <- gsub("_", " ", colnames(in_Network_hubs))
colnames(in_Network_hubs)[colnames(in_Network_hubs) == "Hescs"] <- "hESCs"
colnames(in_Network_hubs)[colnames(in_Network_hubs) == "Ipscs"] <- "iPSCs"
colnames(in_Network_hubs)[colnames(in_Network_hubs) == "Pbmcs t cells"] <- "PBMCs/T cells"

# Perform NMF
library(pheatmap)
library(tidyverse)
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

model_edges <- nmf(as.matrix(in_Network_edges), k = 8)
rownames(model_edges$h) <- paste0("NMF", c(1:8))
colnames(model_edges$h) <- colnames(in_Network_edges)
rownames(model_edges$w) <- rownames(in_Network_edges)
colnames(model_edges$w) <- paste0("NMF", c(1:8))
tissues <- rownames(t(model_edges$h))
factor_names <- colnames(t(model_edges$h))
pheatmap(t(model_edges$h), color = colorRampPalette(brewer.pal(n = 9, name ="Greens"))(100),
         main = "NMF of Edges\npresent in at least 2 contexts", labels_row = make_bold_names(t(model_edges$h), rownames, tissues), 
         labels_col = make_bold_names(t(model_edges$h), colnames, factor_names))



model_hubs <- RcppML::nmf(as.matrix(in_Network_hubs), k = 8)
rownames(model_hubs$h) <- paste0("NMF", c(1:8))
colnames(model_hubs$h) <- colnames(in_Network_hubs)
rownames(model_hubs$w) <- rownames(in_Network_hubs)
colnames(model_hubs$w) <- paste0("NMF", c(1:8))
tissues <- rownames(t(model_hubs$h))
factor_names <- colnames(t(model_hubs$h))
pheatmap(t(model_hubs$h), color = colorRampPalette(brewer.pal(n = 9, name ="Blues"))(100),
         main = "NMF of Hubs\npresent in at least 2 contexts", labels_row = make_bold_names(t(model_hubs$h), rownames, tissues), 
         labels_col = make_bold_names(t(model_hubs$h), colnames, factor_names))

saveRDS(model_hubs, paste0(homeDir, "results/matrix_factorization/model_hubs.rds"))
saveRDS(model_edges, paste0(homeDir, "results/matrix_factorization/model_edges.rds"))

# Save individual factors
for(i in c(1:dim(model_hubs$w)[2])){
  df <- model_hubs$w
  df <- df[df[ ,i] > quantile(df[ ,i], 0.8), ]
  diff <- c()
  for(j in c(1:dim(df)[1])){
    diff[j] <- df[j, i] - max(df[j, -i])
  }
  genes_factor <- rownames(df)[diff >= 5e-5]
  genes_factor <- data.frame(genes_factor)
  dim(genes_factor)
  write.table(genes_factor, paste0(homeDir, "results/matrix_factorization/hubs/factor_",i,".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}


for(i in c(1:dim(model_edges$w)[2])){
  df <- model_edges$w
  df <- df[df[ ,i] > quantile(df[ ,i], 0.9), ]
  diff <- c()
  for(j in c(1:dim(df)[1])){
    diff[j] <- df[j, i] - max(df[j, -i])
  }
  genes_factor <- rownames(df)[diff >= quantile(diff, 0.9)]
  genes_factor <- data.frame(genes_factor)
  dim(genes_factor)
  write.table(genes_factor, paste0(homeDir, "results/matrix_factorization/edges/factor_",i,".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

