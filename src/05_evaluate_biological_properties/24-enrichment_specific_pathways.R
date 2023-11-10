# Description
# This script evaluates for a number of networks the endrichment of genes present in certain GO terms
# Among central network nodes
rm(list = ls())
.libPaths(c("/data/apps/extern/r-packages/4.2.0", 
            "/home/pravich2/R/x86_64-pc-linux-gnu-library/4.2", 
            "/data/apps/extern/spack_on/gcc/9.3.0/r/4.2.0-whb637mlxrrlrjerioexrx2ayqzq7zot/rlib/R/library"))
# Load required libraries
rm(list = ls())
library(Matrix)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(GOfuncR)

home.dir <- "/data/abattle4/prashanthi/recount3_paper/"
dat.dir <- paste0(home.dir, "data/")
res.dir <- paste0(home.dir, "results/weighted_cov_networks/")
plot.dir <- paste0(home.dir, "plots/")

make_contingency <- function(test, bg, set){
  inset_intest <- length(intersect(test,set))
  inset_inbg <- length(intersect(bg,set))
  notinset_intest <- length(test) - inset_intest
  notinset_inbg <- length(bg) - inset_inbg
  cmat = matrix(c(inset_intest, notinset_intest, inset_inbg, notinset_inbg), nrow = 2)
  cmat
}

kegg_enrichment <- function(genelist, background, pathway){
  cont.mat <- make_contingency(genelist , background , pathway)
  res <- fisher.test(cont.mat, alternative = "greater", conf.int = T)
  res$estimate
}

net_odds_pathway <- function(agg_level, lambda, nstudies, GO_term){
  cat(GO_term, "\n")
  net <- readRDS(paste0(res.dir, agg_level, "/net_", nstudies, "/lambda_", 
                        sprintf("%.2f" , lambda), ".rds"))
  geneData <-  readRDS(paste0(dat.dir, "geneData.rds"))
  pc.genes <- read.delim(paste0(dat.dir, "protein_coding.txt"),
                         header = F, stringsAsFactors = F)
  overlapping_genes <- read.delim(paste0(dat.dir, "ensembl_ids_overlapping_genes.txt"),
                                  stringsAsFactors = F)
  geneData <- geneData[geneData$gene_id %in% pc.genes$V2, ]
  net <- as.matrix(net)
  net[lower.tri(net)] <- 0
  diag(net) <- 0
  net[net!=0] <- 1
  rownames(net) <- geneData$gene_name[match(rownames(net), geneData$gene_id)]
  colnames(net) <- geneData$gene_name[match(colnames(net), geneData$gene_id)]
  net_graph <- graph_from_adjacency_matrix(net, mode = "undirected", diag = F)
  degree_df <- degree(net_graph)
  degree_df <- data.frame(degree_df)
  degree_df <- cbind(rownames(degree_df), degree_df)
  colnames(degree_df) <- c("gene", "degree")
  rownames(degree_df) <- c(1:dim(degree_df)[1])
  pathway <- get_anno_genes(GO_term)
  odds_ratio <- c()
  for(degree_thres in unique(c(seq(5, max(degree_df$degree), 5), max(degree_df$degree)))){
    test <- degree_df$gene[degree_df$degree <= degree_thres]
    odds_ratio <- c(odds_ratio, kegg_enrichment(test, geneData$gene_name, pathway$gene))
  }
  res <- data.frame(unique(c(seq(5, max(degree_df$degree), 5), max(degree_df$degree))),
                    odds_ratio)
  colnames(res) <- c("degree", GO_term)
  res
}

# List of consensus pathways
# mitotic_cell_cycle, chromosome_organization, organelle_organization, microtubule_based_processes, 
# cytoskeleton_dependent_cytokinesis, cellular_response_DNA_damage
consensus_pathways <- c("GO:0000278", "GO:0051276", "GO:0033043", "GO:0007017", "GO:0061640", "GO:0006974")
# Neural specific pathways
neural_pathways <- c("GO:0099175", "GO:0007411", "GO:0007417", 
                     "GO:0007420", "GO:0048854", "GO:0035284", "GO:0051610", 
                     "GO:0090494", "GO:0051932", "GO:0097154", "GO:0099536", 
                     "GO:0001963", "GO:0035249","GO:0051932")
# Blood specific pathways
blood_pathways <- c("GO:0002521", "GO:0030595", "GO:0007596", 
                    "GO:0030168", "GO:0042386", "GO:0043249", 
                    "GO:0034101", "GO:0048821",
                    "GO:0050900", "GO:0045321", "GO:0030595", "GO:0045087")
# Skin specific pathways
skin_pathways <- c("GO:0043589", "GO:0043588", "GO:0098773", "GO:0061436", "GO:1903232", 
                   "GO:0097324", "GO:0032963", "GO:0009411")
# Adipose specific pathways
adipose_pathways <- c("GO:0060612", "GO:0070341", "GO:1904606", "GO:0060191", 
                      "GO:0015916", "GO:0006635", "GO:0030497", "GO:0015908", 
                      "GO:0033762")
# Liver specific pathways 
liver_pathways <- c("GO:0072576", "GO:0001889", "GO:0097421", "GO:0005977", 
                    "GO:0070873", "GO:0006641", "GO:0032868", "GO:0006094", 
                    "GO:0032782", "GO:0035622", "GO:0033762")
# Lung specific pathways
lung_pathways <- c("GO:0060503", "GO:0060510", "GO:0060479", "GO:0015671", 
                   "GO:0060428", "GO:0061145", "GO:0061141", "GO:0014916", 
                   "GO:0070254")


all_consensus_res <- lapply(consensus_pathways, function(GO_term){
  net_odds_pathway("all_consensus", 0.18, 966, GO_term)
})
normal_consensus_res <- lapply(consensus_pathways, function(GO_term){
  net_odds_pathway("normal_consensus", 0.18, 629, GO_term)
})
cancer_consensus_res <- lapply(consensus_pathways, function(GO_term){
  net_odds_pathway("cancer_consensus", 0.24, 386, GO_term)
})
cns_consensus_res <- lapply(consensus_pathways, function(GO_term){
  net_odds_pathway("central_nervous_system/all", 0.28, 53, GO_term)
})
blood_consensus_res <- lapply(consensus_pathways, function(GO_term){
  net_odds_pathway("blood/all", 0.24, 65, GO_term)
})


all_consensus_neural <- lapply(neural_pathways, function(GO_term){
  net_odds_pathway("all_consensus", 0.18, 967, GO_term)
})
normal_consensus_neural <- lapply(neural_pathways, function(GO_term){
  net_odds_pathway("normal_consensus", 0.18, 629, GO_term)
})
cns_GTEx_neural <- lapply(neural_pathways, function(GO_term){
  net_odds_pathway("central_nervous_system/GTEx", 0.32, 13, GO_term)
})
cns_consensus_neural <- lapply(neural_pathways, function(GO_term){
  net_odds_pathway("central_nervous_system/all", 0.28, 53, GO_term)
})
blood_consensus_neural <- lapply(neural_pathways, function(GO_term){
  net_odds_pathway("blood/all", 0.24, 65, GO_term)
})


all_consensus_blood <- lapply(blood_pathways, function(GO_term){
  net_odds_pathway("all_consensus", 0.18, 967, GO_term)
})
normal_consensus_blood <- lapply(blood_pathways, function(GO_term){
  net_odds_pathway("normal_consensus", 0.18, 629, GO_term)
})
blood_GTEx_blood <- lapply(blood_pathways, function(GO_term){
  net_odds_pathway("blood/GTEx", 0.28, 1, GO_term)
})
blood_consensus_blood <- lapply(blood_pathways, function(GO_term){
  net_odds_pathway("blood/all", 0.24, 65, GO_term)
})
cns_consensus_blood <- lapply(blood_pathways, function(GO_term){
  net_odds_pathway("central_nervous_system/all", 0.28, 53, GO_term)
})

all_consensus_skin <- lapply(skin_pathways, function(GO_term){
  net_odds_pathway("all_consensus", 0.18, 967, GO_term)
})
normal_consensus_skin <- lapply(skin_pathways, function(GO_term){
  net_odds_pathway("normal_consensus", 0.18, 629, GO_term)
})
skin_GTEx_skin <- lapply(skin_pathways, function(GO_term){
  net_odds_pathway("skin/GTEx", 0.26, 2, GO_term)
})
skin_consensus_skin <- lapply(skin_pathways, function(GO_term){
  net_odds_pathway("skin/all", 0.26, 20, GO_term)
})
blood_consensus_skin <- lapply(skin_pathways, function(GO_term){
  net_odds_pathway("blood/all", 0.24, 65, GO_term)
})
cns_consensus_skin <- lapply(skin_pathways, function(GO_term){
  net_odds_pathway("central_nervous_system/all", 0.28, 53, GO_term)
})


all_consensus_adipose <- lapply(adipose_pathways, function(GO_term){
  net_odds_pathway("all_consensus", 0.18, 967, GO_term)
})
normal_consensus_adipose <- lapply(adipose_pathways, function(GO_term){
  net_odds_pathway("normal_consensus", 0.18, 629, GO_term)
})
adipose_GTEx_adipose <- lapply(adipose_pathways, function(GO_term){
  net_odds_pathway("adipose/GTEx", 0.28, 2, GO_term)
})
adipose_consensus_adipose <- lapply(adipose_pathways, function(GO_term){
  net_odds_pathway("adipose/all", 0.26, 11, GO_term)
})
blood_consensus_adipose <- lapply(adipose_pathways, function(GO_term){
  net_odds_pathway("blood/all", 0.24, 65, GO_term)
})
cns_consensus_adipose <- lapply(adipose_pathways, function(GO_term){
  net_odds_pathway("central_nervous_system/all", 0.28, 53, GO_term)
})



all_consensus_lung <- lapply(lung_pathways, function(GO_term){
  net_odds_pathway("all_consensus", 0.18, 967, GO_term)
})
normal_consensus_lung <- lapply(lung_pathways, function(GO_term){
  net_odds_pathway("normal_consensus", 0.18, 629, GO_term)
})
lung_GTEx_lung <- lapply(lung_pathways, function(GO_term){
  net_odds_pathway("lung/GTEx", 0.30, 1, GO_term)
})
lung_consensus_lung <- lapply(lung_pathways, function(GO_term){
  net_odds_pathway("lung/all", 0.28, 10, GO_term)
})
blood_consensus_lung <- lapply(lung_pathways, function(GO_term){
  net_odds_pathway("blood/all", 0.24, 65, GO_term)
})
cns_consensus_lung <- lapply(lung_pathways, function(GO_term){
  net_odds_pathway("central_nervous_system/all", 0.28, 53, GO_term)
})


all_consensus_liver <- lapply(liver_pathways, function(GO_term){
  net_odds_pathway("all_consensus", 0.18, 967, GO_term)
})
normal_consensus_liver <- lapply(liver_pathways, function(GO_term){
  net_odds_pathway("normal_consensus", 0.18, 629, GO_term)
})
liver_GTEx_liver <- lapply(liver_pathways, function(GO_term){
  net_odds_pathway("liver/GTEx", 0.38, 1, GO_term)
})
liver_consensus_liver <- lapply(liver_pathways, function(GO_term){
  net_odds_pathway("liver/all", 0.24, 28, GO_term)
})
blood_consensus_liver <- lapply(liver_pathways, function(GO_term){
  net_odds_pathway("blood/all", 0.24, 65, GO_term)
})
cns_consensus_liver <- lapply(liver_pathways, function(GO_term){
  net_odds_pathway("central_nervous_system/all", 0.28, 53, GO_term)
})


convert_list_to_df <- function(ilist, index){
  df <- ilist[[1]]
  for(i in c(2:length(ilist))){
    idf <- ilist[[i]]
    idf <- idf[match(df[ ,index], idf[ ,index]), ]
    idf[ ,index] <- NULL
    df <- cbind(df, idf)
  }
  df
}

all_consensus_res <- convert_list_to_df(all_consensus_res, "degree")
normal_consensus_res <- convert_list_to_df(normal_consensus_res, "degree")
cancer_consensus_res <- convert_list_to_df(cancer_consensus_res, "degree")
cns_consensus_res <- convert_list_to_df(cns_consensus_res, "degree")
blood_consensus_res <- convert_list_to_df(blood_consensus_res, "degree")

all_consensus_neural <- convert_list_to_df(all_consensus_neural, "degree")
normal_consensus_neural <- convert_list_to_df(normal_consensus_neural, "degree")
cns_GTEx_neural <- convert_list_to_df(cns_GTEx_neural, "degree")
cns_consensus_neural <- convert_list_to_df(cns_consensus_neural, "degree")
blood_consensus_neural <- convert_list_to_df(blood_consensus_neural, "degree")

all_consensus_blood <- convert_list_to_df(all_consensus_blood, "degree")
normal_consensus_blood <- convert_list_to_df(normal_consensus_blood, "degree")
blood_GTEx_blood <- convert_list_to_df(blood_GTEx_blood, "degree")
blood_consensus_blood <- convert_list_to_df(blood_consensus_blood, "degree")
cns_consensus_blood <- convert_list_to_df(cns_consensus_blood, "degree")

all_consensus_skin <- convert_list_to_df(all_consensus_skin, "degree")
normal_consensus_skin <- convert_list_to_df(normal_consensus_skin, "degree")
skin_GTEx_skin <- convert_list_to_df(skin_GTEx_skin, "degree")
skin_consensus_skin <- convert_list_to_df(skin_consensus_skin, "degree")
blood_consensus_skin <- convert_list_to_df(blood_consensus_skin, "degree")
cns_consensus_skin <- convert_list_to_df(cns_consensus_skin, "degree")


all_consensus_adipose <- convert_list_to_df(all_consensus_adipose, "degree")
normal_consensus_adipose <- convert_list_to_df(normal_consensus_adipose, "degree")
adipose_GTEx_adipose <- convert_list_to_df(adipose_GTEx_adipose, "degree")
adipose_consensus_adipose <- convert_list_to_df(adipose_consensus_adipose, "degree")
blood_consensus_adipose <- convert_list_to_df(blood_consensus_adipose, "degree")
cns_consensus_adipose <- convert_list_to_df(cns_consensus_adipose, "degree")

all_consensus_liver <- convert_list_to_df(all_consensus_liver, "degree")
normal_consensus_liver <- convert_list_to_df(normal_consensus_liver, "degree")
liver_GTEx_liver <- convert_list_to_df(liver_GTEx_liver, "degree")
liver_consensus_liver <- convert_list_to_df(liver_consensus_liver, "degree")
blood_consensus_liver <- convert_list_to_df(blood_consensus_liver, "degree")
cns_consensus_liver <- convert_list_to_df(cns_consensus_liver, "degree")

all_consensus_lung <- convert_list_to_df(all_consensus_lung, "degree")
normal_consensus_lung <- convert_list_to_df(normal_consensus_lung, "degree")
lung_GTEx_lung <- convert_list_to_df(lung_GTEx_lung, "degree")
lung_consensus_lung <- convert_list_to_df(lung_consensus_lung, "degree")
blood_consensus_lung <- convert_list_to_df(blood_consensus_lung, "degree")
cns_consensus_lung <- convert_list_to_df(cns_consensus_lung, "degree")

all_consensus_res$type <- "All samples"
normal_consensus_res$type <- "Non-cancerous samples"
cancer_consensus_res$type <- "Cancerous samples"
cns_consensus_res$type <- "CNS"
blood_consensus_res$type <- "Blood"


all_consensus_neural$type <- "All samples"
normal_consensus_neural$type <- "Non-cancerous samples"
cns_consensus_neural$type <- "CNS"
cns_GTEx_neural$type <- "CNS (GTEx)"
blood_consensus_neural$type <- "Blood"

all_consensus_skin$type <- "All samples"
normal_consensus_skin$type <- "Non-cancerous samples"
skin_consensus_skin$type <- "Skin"
skin_GTEx_skin$type <- "Skin (GTEx)"
blood_consensus_skin$type <- "Blood"
cns_consensus_skin$type <- "CNS"

all_consensus_blood$type <- "All samples"
normal_consensus_blood$type <- "Non-cancerous samples"
blood_consensus_blood$type <- "Blood"
blood_GTEx_blood$type <- "Blood (GTEx)"
cns_consensus_blood$type <- "CNS"

all_consensus_adipose$type <- "All samples"
normal_consensus_adipose$type <- "Non-cancerous samples"
adipose_consensus_adipose$type <- "Adipose"
adipose_GTEx_adipose$type <- "Adipose (GTEx)"
blood_consensus_adipose$type <- "Blood"
cns_consensus_adipose$type <- "CNS"

all_consensus_liver$type <- "All samples"
normal_consensus_liver$type <- "Non-cancerous samples"
liver_consensus_liver$type <- "Liver"
liver_GTEx_liver$type <- "Liver (GTEx)"
blood_consensus_liver$type <- "Blood"
cns_consensus_liver$type <- "CNS"


all_consensus_lung$type <- "All samples"
normal_consensus_lung$type <- "Non-cancerous samples"
lung_consensus_lung$type <- "Lung"
lung_GTEx_lung$type <- "Lung (GTEx)"
blood_consensus_lung$type <- "Blood"
cns_consensus_lung$type <- "CNS"

consensus_res <- rbind(all_consensus_res, normal_consensus_res, cancer_consensus_res ,cns_consensus_res, blood_consensus_res)
consensus_res$type <- as.factor(consensus_res$type)

neural_res <- rbind(all_consensus_neural, normal_consensus_neural, cns_GTEx_neural, cns_consensus_neural, blood_consensus_neural)
neural_res$type <- as.factor(neural_res$type)

blood_res <- rbind(all_consensus_blood, normal_consensus_blood, blood_GTEx_blood, blood_consensus_blood, cns_consensus_blood)
blood_res$type <- as.factor(blood_res$type)

skin_res <- rbind(all_consensus_skin, normal_consensus_skin, skin_GTEx_skin, skin_consensus_skin, cns_consensus_skin, blood_consensus_skin)
skin_res$type <- as.factor(skin_res$type)

adipose_res <- rbind(all_consensus_adipose, normal_consensus_adipose, adipose_GTEx_adipose, adipose_consensus_adipose, cns_consensus_adipose, blood_consensus_adipose)
adipose_res$type <- as.factor(adipose_res$type)

lung_res <- rbind(all_consensus_lung, normal_consensus_lung, lung_GTEx_lung, lung_consensus_lung, cns_consensus_lung, blood_consensus_lung)
lung_res$type <- as.factor(lung_res$type)

liver_res <- rbind(all_consensus_liver, normal_consensus_liver, liver_GTEx_liver, liver_consensus_liver, cns_consensus_liver, blood_consensus_liver)
liver_res$type <- as.factor(liver_res$type)


saveRDS(consensus_res, paste0(home.dir, "results/enrichment_plots/consensus_res.rds"))
saveRDS(neural_res, paste0(home.dir, "results/enrichment_plots/neural_res.rds"))
saveRDS(skin_res, paste0(home.dir, "results/enrichment_plots/skin_res.rds"))
saveRDS(blood_res, paste0(home.dir, "results/enrichment_plots/blood_res.rds"))
saveRDS(adipose_res, paste0(home.dir, "results/enrichment_plots/adipose_res.rds"))
saveRDS(liver_res, paste0(home.dir, "results/enrichment_plots/liver_res.rds"))
saveRDS(lung_res, paste0(home.dir, "results/enrichment_plots/lung_res.rds"))

pathways <- c(consensus_pathways, neural_pathways, skin_pathways, blood_pathways, adipose_pathways, liver_pathways, lung_pathways)
pathways_df <- get_names(pathways)
saveRDS(pathways_df, paste0(home.dir, "results/enrichment_plots/pathways.rds"))


