# Description
# This code reads networks of varying densities from a particular aggregation setup
# and held-out test data to compute held out log-likelihood 
# To be used for evaluating the held-out log-likelihood of 
# 1. GTEx
# 2. SRA non-cancer networks
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
lambda <- seq(0.08, 1.00, 0.02)
lambda <- sprintf("%.2f", lambda)

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
res.dir <- paste0(homeDir, "results/")
dat.dir <- paste0(homeDir, "/data")

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

edge.permuted.list <- function(net){
  edge.unformatted <- lapply(net, function(inet){
    if(!is.null(dim(inet))){
      inet <- as.matrix(inet)
      ngenes <- dim(inet)[1]
      #row.names.org <- paste0(rep("gene", ngenes), as.character(1:ngenes))
      row.names.org <- rownames(inet)
      set.seed(123)
      p.row.names <- sample(row.names.org, size = length(row.names.org))
      rownames(inet) <- p.row.names
      colnames(inet) <- p.row.names
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

compute_contingency <- function(A,permA, true_edges){
  contingency_table <- matrix(rep(0,4), nrow = 2, ncol = 2)
  rownames(contingency_table) <- c("In_network", "Out_network")
  colnames(contingency_table) <- c("In_pathway", "Out_pathway")
  contingency_table[1,1] <- length(intersect(A, true_edges))
  contingency_table[1,2] <- length(A) - length(intersect(A, true_edges))
  contingency_table[2,1] <- length(intersect(permA, true_edges))
  contingency_table[2,2] <- length(permA) - length(intersect(permA, true_edges))
  contingency_table
}

# This script analyses networks
inputArgs <- commandArgs(TRUE)

study <- inputArgs[1]
agg_type <- inputArgs[2]
nStudies <- inputArgs[3]

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

# Find edges
inferred.networks <- edge.list(net)
permuted.networks <- edge.permuted.list(net)

print("Inferred and permuted networks computed")

n.edges <- lapply(inferred.networks, function(iedge){
  length(iedge)
})
n.edges <- unlist(n.edges)

# Read in the true networks
true_edges <- readRDS(paste0(dat.dir, "canonical_pathways_edgelist.rds"))
gene1 <- gsub("_.*", "", true_edges)
gene2 <- gsub(".*_", "", true_edges)
true_edges.df <- data.frame(gene1, gene2)
# Restrict the true edges only to those gene pairs where both genes are found in 
# the genes that we do inference over
true_edges.df <- true_edges.df[true_edges.df$gene1 %in% rownames(net[[1]]), ]
true_edges.df <- true_edges.df[true_edges.df$gene2 %in% rownames(net[[1]]), ]
true_edges <- paste(true_edges.df$gene1, true_edges.df$gene2, sep = "_")

# Create a network with all possible edges
all_genes <- unique(c(rownames(net[[1]]), true_edges.df$gene1, true_edges.df$gene2))
ngenes <- length(all_genes)
all.net <- matrix(1, nrow = ngenes, ncol = ngenes)
colnames(all.net) <- all_genes[order(all_genes)]
rownames(all.net) <- all_genes[order(all_genes)]
all.net[lower.tri(all.net, diag = T)] <- NA
all.edges <- reshape2::melt(all.net, na.rm = T)
all.edgelist <- paste(all.edges[ ,1], all.edges[ ,2], sep = "_")
print("All edges computed")

frac.in.true_edges <- lapply(inferred.networks, function(iedge){
  length(intersect(iedge, true_edges))/ length(iedge)
})
frac.in.true_edges <- unlist(frac.in.true_edges)

true_positives <- lapply(inferred.networks, function(iedge){
  length(intersect(iedge, true_edges))
})

false_positives <- lapply(inferred.networks, function(iedge){
  length(iedge) - length(intersect(iedge, true_edges))
})

false_negatives <- lapply(inferred.networks, function(iedge){
  length(true_edges) - length(intersect(iedge, true_edges))
})

true_negatives <- lapply(inferred.networks, function(iedge){
  cat(length(iedge), "\n")
  length(intersect(all.edgelist[!all.edgelist %in% true_edges],all.edgelist[!all.edgelist %in% iedge]))
})

true_positives <- unlist(true_positives)
false_positives <- unlist(false_positives)
false_negatives <- unlist(false_negatives)
true_negatives <- unlist(true_negatives)

recall <- true_positives/(true_positives + false_negatives)
precision <- true_positives/(true_positives + false_positives)

print("Num.  edges and fraction in true_edges computed")
odds_ratio <- c()
lower_odds_ratio <- c()
higher_odds_ratio <- c()
for(i in c(1:length(true_positives))){
  cmat <- matrix(0, nrow = 2, ncol = 2)
  cmat[1, 1] <- true_positives[i]
  cmat[1, 2] <- false_positives[i]
  cmat[2, 1] <- false_negatives[i]
  cmat[2, 2] <- true_negatives[i]
  odds_test <- fisher.test(cmat)
  odds_ratio[i] <- odds_test$estimate
  lower_odds_ratio[i] <- odds_test$conf.int[1]
  higher_odds_ratio[i] <- odds_test$conf.int[2]
}

ct.list <- list()
for(i in c(1:length(inferred.networks))){
  ct.list[[i]] <- compute_contingency(inferred.networks[[i]], permuted.networks[[i]], true_edges)
}

odds.ratio <- lapply(ct.list, function(ct){
  fisher.test(ct)$estimate
})

p.value <- lapply(ct.list, function(ct){
  fisher.test(ct)$p.value
})

odds.ratio <- unlist(odds.ratio)
p.value <- unlist(p.value)


res <- data.frame(lambda, n.edges, frac.in.true_edges,true_positives, false_positives,
                  false_negatives, true_negatives, lower_odds_ratio, odds_ratio, higher_odds_ratio ,precision, recall, odds.ratio, p.value)
colnames(res) <- c("lambda", "n_edges", "frac_in_true", "true_positives", "false_positives", "false_negatives", 
                   "true_negatives", "lower_odds_ratio", "odds_ratio", "higher_odds_ratio", "precision", "recall", "odds_ratio_permuted", "p_value_permuted")
res$F1 <- (2*res$precision*res$recall)/(res$precision + res$recall)

saveRDS(res, paste0(res.dir, agg_type, "/", study, "/net_", nStudies, "/oddsRatio_", nStudies, ".rds"))



