# Description
# This script takes in networks and computes the excess overlap 
# of network central genes with known gene sets
# Output produced includes the centrality for a particular set of values of the penalization parameter lambda
# Excess overlap of the corresponding networks are computed and saved as well. 
rm(list = ls())
.libPaths(c("/home/pravich2/R/x86_64-pc-linux-gnu-library/4.2", 
            "/data/apps/extern/r-packages/4.2.0", 
            "/data/apps/extern/spack_on/gcc/9.3.0/r/4.2.0-whb637mlxrrlrjerioexrx2ayqzq7zot/rlib/R/library"))
library(Matrix)
library(cowplot)
library(reshape2)
library(mvtnorm)
library(igraph)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(caret)
library(scales)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(AnnotationHub)


excess_overlap <- function(gene_set, gene_net_decile){
  genes_deciles <- list()
  P_d <- c()
  P_tot <- length(intersect(gene_set , gene_net_decile[ ,1]))/length(gene_net_decile[,1])
  for(i in c(1:max(gene_net_decile[ ,2]))){
    genes_deciles[[i]] <- gene_net_decile[gene_net_decile[,2] == i, 1]
    P_d[i] <- length(intersect(gene_set, genes_deciles[[i]]))/length(genes_deciles[[i]])
  }
  P_d/P_tot
}

se <- function(gene_set, gene_net_decile){
  genes_deciles <- list()
  P_d <- c()
  se <- c()
  P_tot <- length(intersect(gene_set , gene_net_decile[ ,1]))/length(gene_net_decile[,1])
  for(i in c(1:max(gene_net_decile[ ,2]))){
    genes_deciles[[i]] <- gene_net_decile[gene_net_decile[,2] == i, 1]
    P_d[i] <- length(intersect(gene_set, genes_deciles[[i]]))/length(genes_deciles[[i]])
    se[i] <- sqrt((P_d[i]*(1 - P_d[i])) / length(genes_deciles[[i]]))
  }
  se/P_tot
}

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
compute_enrichment <- function(study, nstudies, lambda){
  agg_type <- "weighted_cov_networks"
  net_list <- lapply(lambda, function(ilambda){
    res.fn <- paste0(homeDir, "results/", agg_type, "/", study, "/", 
                     "net_", nstudies, "/lambda_", ilambda, ".rds")
    net <- readRDS(res.fn)
  })
  
  dat.dir <- paste0(homeDir, "data/")
  protein_kinases <- read.table(paste0(dat.dir, "protein_kinases.txt"))
  gene_sets <- read.csv(paste0(dat.dir, "gene_sets.csv"))
  median_cis <- read.table(paste0(dat.dir, "gene_median_cispval.tsv"), header = T)
  all_tissues_egene <- read.delim(paste0(dat.dir, "all_tissues_egene.txt"), header = T)
  
  all.genes <- mapIds(org.Hs.eg.db, as.character(gene_sets$All.genes[!is.na(gene_sets$All.genes)]), 'SYMBOL', 'ENTREZID')
  MGI.essential <- mapIds(org.Hs.eg.db, as.character(gene_sets$MGI.essential[!is.na(gene_sets$MGI.essential)]), 'SYMBOL', 'ENTREZID')
  autosomal.dominant <-  mapIds(org.Hs.eg.db, as.character(gene_sets$Autosomal.dominant[!is.na(gene_sets$Autosomal.dominant)]), 'SYMBOL', 'ENTREZID')
  Haploinsufficient <- mapIds(org.Hs.eg.db, as.character(gene_sets$Haploinsufficient[!is.na(gene_sets$Haploinsufficient)]), 'SYMBOL', 'ENTREZID')
  High.pLI <- mapIds(org.Hs.eg.db, as.character(gene_sets$High.pLI[!is.na(gene_sets$High.pLI)]), 'SYMBOL', 'ENTREZID')
  High.Shet <- mapIds(org.Hs.eg.db, as.character(gene_sets$High.Shet[!is.na(gene_sets$High.Shet)]), 'SYMBOL', 'ENTREZID')
  High.Phi <- mapIds(org.Hs.eg.db, as.character(gene_sets$High.Phi[!is.na(gene_sets$High.Phi)]), 'SYMBOL', 'ENTREZID')
  High.missense.z <- mapIds(org.Hs.eg.db, as.character(gene_sets$High.missense.z[!is.na(gene_sets$High.missense.z)]), 'SYMBOL', 'ENTREZID')
  ClinVar <- mapIds(org.Hs.eg.db, as.character(gene_sets$ClinVar[!is.na(gene_sets$ClinVar)]), 'SYMBOL', 'ENTREZID')
  OMIM <- mapIds(org.Hs.eg.db, as.character(gene_sets$OMIM[!is.na(gene_sets$OMIM)]), 'SYMBOL', 'ENTREZID')
  GWAS.nearest <- mapIds(org.Hs.eg.db, as.character(gene_sets$GWAS.nearest[!is.na(gene_sets$GWAS.nearest)]), 'SYMBOL', 'ENTREZID')
  TF <- mapIds(org.Hs.eg.db, as.character(gene_sets$TF[!is.na(gene_sets$TF)]), 'SYMBOL', 'ENTREZID')
  DrugBank <- mapIds(org.Hs.eg.db, as.character(gene_sets$DrugBank[!is.na(gene_sets$DrugBank)]), 'SYMBOL', 'ENTREZID')
  High.EDS <- mapIds(org.Hs.eg.db, as.character(gene_sets$High.EDS[!is.na(gene_sets$High.EDS)]), 'SYMBOL', 'ENTREZID')
  eQTL0.deficient <- mapIds(org.Hs.eg.db, as.character(gene_sets$eQTL0.deficient[!is.na(gene_sets$eQTL0.deficient)]), 'SYMBOL', 'ENTREZID')
  Olfactory <- mapIds(org.Hs.eg.db, as.character(gene_sets$Olfactory[!is.na(gene_sets$Olfactory)]), 'SYMBOL', 'ENTREZID')
  
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
  
  res_df_list <- list()
  for(i in c(1:length(net_list))){
    net <- net_list[[i]]
    graph_distance <- get_graph_distance(net)
    graph_cov <- get_graph_cov(net)
    centrality.df <- data.frame(matrix(ncol = 8, nrow = length(V(graph_distance))))
    colnames(centrality.df) <- c("node", "degree", "maximum_weight", 
                                 "strength", "betweeness", "closeness", 
                                 "eigen_centrality", "page_rank")
    centrality.df$node <- V(graph_distance)$label
    centrality.df$degree <- degree(graph_cov, v = V(graph_cov))
    centrality.df$strength <- strength(graph_cov, v = V(graph_cov))
    centrality.df$maximum_weight <- ldply(V(graph_cov), function(v) range(incident(graph_cov, v, mode='total')$weight))$V2
    centrality.df$betweeness <- betweenness(graph_distance, v = V(graph_distance))
    sample.dist.matrix <- distances(graph_distance, v = V(graph_distance), to = V(graph_distance), algorithm = "dijkstra")
    sample.closeness.matrix <- 1/sample.dist.matrix
    diag(sample.closeness.matrix) <- 0
    centrality.df$closeness <- rowSums(sample.closeness.matrix)
    centrality.df$page_rank <- page_rank(graph_cov, v = V(graph_cov))$vector
    centrality.df$eigen_centrality <- eigen_centrality(graph_cov)$vector
    
    for(icol in c(2:8)){
      centrality.df[ ,icol] <- rescale(centrality.df[ ,icol])
      centrality.df[centrality.df[ ,icol] < 0, icol] <- 0
    }
    saveRDS(centrality.df, paste0(homeDir, "results/", agg_type, "/", study, "/", 
                                  "net_", nstudies, "/centrality_", lambda[i], ".rds"))
    type <- "closeness"
    centrality.df <- centrality.df[ ,colnames(centrality.df) %in% c("node", type)]
    colnames(centrality.df) <- c("node", "centrality_measure")
    geneSymbols <- c()
    for(i in c(1:dim(centrality.df)[1])){
      geneSymbols[i] <- geneData$gene_name[geneData$gene_id == centrality.df$node[i]]
    }
    centrality.df$node <- geneSymbols
    # Compute deciles
    nGroups <- 4
    zero_centrality <- centrality.df[centrality.df$centrality_measure == 0, ]
    nonzero_centrality <- centrality.df[centrality.df$centrality_measure > 0, ]
    
    decile.df <- data.frame(matrix(ncol = dim(nonzero_centrality)[2], nrow = dim(nonzero_centrality)[1]))
    colnames(decile.df) <- colnames(nonzero_centrality)
    decile.df$node <- nonzero_centrality$node
    decile.df$centrality_measure <- ntile(nonzero_centrality$centrality_measure, nGroups)
    decile.df <- rbind(decile.df, zero_centrality)
    decile.df$centrality_measure <- decile.df$centrality_measure + 1
    
    excess_overlap.df <- data.frame(matrix(ncol = 15, nrow = max(decile.df$centrality_measure)))
    colnames(excess_overlap.df) <- c("MGI essential", "Autosomal dominant", 
                                     "Haploinsufficient", "ClinVar", "OMIM", 
                                     "GWAS nearest", "Transcription factor", "Drug bank", "eQTL0 deficient", 
                                     "High EDS", "High PLI", "High Shet", "High Phi", "High missense z-score", "All genes")
    excess_overlap.df[ ,1] <- excess_overlap(MGI.essential, decile.df)
    excess_overlap.df[ ,2] <- excess_overlap(autosomal.dominant, decile.df)
    excess_overlap.df[ ,3] <- excess_overlap(Haploinsufficient, decile.df)
    excess_overlap.df[ ,4] <- excess_overlap(ClinVar, decile.df)
    excess_overlap.df[ ,5] <- excess_overlap(OMIM, decile.df)
    excess_overlap.df[ ,6] <- excess_overlap(GWAS.nearest, decile.df)
    excess_overlap.df[ ,7] <- excess_overlap(TF, decile.df)
    excess_overlap.df[ ,8] <- excess_overlap(DrugBank, decile.df)
    excess_overlap.df[ ,9] <- excess_overlap(eQTL0.deficient, decile.df)
    excess_overlap.df[ ,10] <- excess_overlap(High.EDS, decile.df)
    excess_overlap.df[ ,11] <- excess_overlap(High.pLI, decile.df)
    excess_overlap.df[ ,12] <- excess_overlap(High.Shet, decile.df)
    excess_overlap.df[ ,13] <- excess_overlap(High.Phi, decile.df)
    excess_overlap.df[ ,14] <- excess_overlap(High.missense.z, decile.df)
    excess_overlap.df[ ,15] <- excess_overlap(all.genes, decile.df)
    excess_overlap.df$Decile <- c(1:max(decile.df$centrality_measure))
    excess_overlap.df <- melt(excess_overlap.df, id.vars = "Decile")
    
    se.df <- data.frame(matrix(ncol = 15, nrow = max(decile.df$centrality_measure)))
    colnames(se.df) <- c("MGI essential", "Autosomal dominant", 
                         "Haploinsufficient", "ClinVar", "OMIM", 
                         "GWAS nearest", "Transcription factor", "Drug bank", "eQTL0 deficient", 
                         "High EDS", "High PLI", "High Shet", "High Phi", "High missense z-score", "All genes")
    se.df[ ,1] <- se(MGI.essential, decile.df)
    se.df[ ,2] <- se(autosomal.dominant, decile.df)
    se.df[ ,3] <- se(Haploinsufficient, decile.df)
    se.df[ ,4] <- se(ClinVar, decile.df)
    se.df[ ,5] <- se(OMIM, decile.df)
    se.df[ ,6] <- se(GWAS.nearest, decile.df)
    se.df[ ,7] <- se(TF, decile.df)
    se.df[ ,8] <- se(DrugBank, decile.df)
    se.df[ ,9] <- se(eQTL0.deficient, decile.df)
    se.df[ ,10] <- se(High.EDS, decile.df)
    se.df[ ,11] <- se(High.pLI, decile.df)
    se.df[ ,12] <- se(High.Shet, decile.df)
    se.df[ ,13] <- se(High.Phi, decile.df)
    se.df[ ,14] <- se(High.missense.z, decile.df)
    se.df[ ,15] <- se(all.genes, decile.df)
    se.df$Decile <- c(1:max(decile.df$centrality_measure))
    se.df <- melt(se.df, id.vars = "Decile")
    
    colnames(se.df) <- c("Decile", "Gene_set", "SE")
    colnames(excess_overlap.df) <- c("Decile", "Gene_set", "Enrichment")
    res_df <- merge(excess_overlap.df, se.df)
    res_df_list[[i]] <- res_df
  }
  res_df_list }

normal_consensus_res <- compute_enrichment("normal_consensus", "629", c("0.14", "0.16", "0.18", "0.20"))
all_consensus_res <- compute_enrichment("all_consensus", "966", c("0.14", "0.16", "0.18", "0.20"))
cancer_consensus_res <- compute_enrichment("cancer_consensus", "386", c("0.20", "0.22", "0.24", "0.26"))
saveRDS(cancer_consensus_res, paste0(homeDir, "results/weighted_cov_networks/cancer_consensus/net_386/node_closeness_enrichment.rds"))
saveRDS(all_consensus_res, paste0(homeDir, "results/weighted_cov_networks/all_consensus/net_966/node_closeness_enrichment.rds"))
saveRDS(normal_consensus_res, paste0(homeDir, "results/weighted_cov_networks/normal_consensus/net_629/node_closeness_enrichment.rds"))

blood_consensus_res <- compute_enrichment("blood/all", "65", c("0.18", "0.20", "0.22", "0.24", "0.26"))
cns_consensus_res <- compute_enrichment("central_nervous_system/all", "53", c("0.24", "0.26", "0.28", "0.30", "0.32"))
adipose_consensus_res <- compute_enrichment("adipose/all", "11", c("0.20", "0.22", "0.24", "0.26", "0.28"))
liver_consensus_res <- compute_enrichment("liver/all", "28", c("0.20", "0.22", "0.24", "0.26", "0.28"))
lung_consensus_res <- compute_enrichment("lung/all", "10", c( "0.24", "0.26", "0.28", "0.30", "0.32"))
skin_consensus_res <- compute_enrichment("skin/all", "20", c("0.22", "0.24", "0.26", "0.28"))

saveRDS(blood_consensus_res, paste0(homeDir, "results/weighted_cov_networks/blood/all/net_65/node_closeness_enrichment.rds"))
saveRDS(cns_consensus_res, paste0(homeDir, "results/weighted_cov_networks/central_nervous_system/all/net_53/node_closeness_enrichment.rds"))
saveRDS(skin_consensus_res, paste0(homeDir, "results/weighted_cov_networks/skin/all/net_20/node_closeness_enrichment.rds"))
saveRDS(liver_consensus_res, paste0(homeDir, "results/weighted_cov_networks/liver/all/net_28/node_closeness_enrichment.rds"))
saveRDS(adipose_consensus_res, paste0(homeDir, "results/weighted_cov_networks/adipose/all/net_11/node_closeness_enrichment.rds"))
saveRDS(lung_consensus_res, paste0(homeDir, "results/weighted_cov_networks/lung/all/net_10/node_closeness_enrichment.rds"))




