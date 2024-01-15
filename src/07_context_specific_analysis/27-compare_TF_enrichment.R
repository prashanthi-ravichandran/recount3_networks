rm(list = ls())
.libPaths(c("/home/pravich2/rlibs/4.0.2/gcc/9.3.0", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-4.0.2-amdvcpog4ugspqwwx3ari7pzkmckelu6/rlib/R/library"))
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(rstatix)
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
# net <- lapply(net, function(inet){
#   gene_ids <- rownames(inet)
#   geneSymbols <- c()
#   for(i in c(1:length(gene_ids))){
#     geneSymbols[i] <- geneData$gene_name[geneData$gene_id == gene_ids[i]]
#   }
#   rownames(inet) <- geneSymbols
#   colnames(inet) <- geneSymbols
#   inet
# })

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

tissue_specific_TFs <- read.csv(paste0(datDir, "tissue_specific_TFs.csv"))
tissue_specific_TFs$X <- NULL

all_TFs <- unique(unlist(str_split(tissue_specific_TFs$Transcription.Factors..Ensembl., ",")))

genes <- rownames(net[[1]])
for(i in c(2:length(net))){
  genes <- intersect(genes, rownames(net[[i]]))
}

genes <- gsub("\\..*","",genes)
non_TFs <- genes[!genes %in% all_TFs]

# Brain TFs
brain_centrality <- degreeCentrality_list[["central_nervous_system"]]
cardiac_centrality <- degreeCentrality_list[["cardiac"]]
skin_centrality <- degreeCentrality_list[["skin"]]
blood_centrality <- degreeCentrality_list[["blood"]]
lung_centrality <- degreeCentrality_list[["lung"]]
muscle_centrality <- degreeCentrality_list[["skeletal_muscle"]]
pancreas_centrality <- degreeCentrality_list[["pancreas"]]
all_centrality <- degreeCentrality_list[["all_consensus"]]
normal_centrality <- degreeCentrality_list[["normal_consensus"]]

brain_TFs <- tissue_specific_TFs$Transcription.Factors..Ensembl.[tissue_specific_TFs$Tissue == "All Brain Tissues"]
brain_TFs <- c(brain_TFs, tissue_specific_TFs$Transcription.Factors..Ensembl.[tissue_specific_TFs$Tissue == "Spinal Cord"])
cardiac_TFs <- tissue_specific_TFs$Transcription.Factors..Ensembl.[tissue_specific_TFs$Tissue == "Heart (Left Ventricle and Atrial Appendage)"]
lung_TFs <- tissue_specific_TFs$Transcription.Factors..Ensembl.[tissue_specific_TFs$Tissue == "Lung"]
muscle_TFs <- tissue_specific_TFs$Transcription.Factors..Ensembl.[tissue_specific_TFs$Tissue == "Muscle - Skeletal"]
pancreas_TFs <- tissue_specific_TFs$Transcription.Factors..Ensembl.[tissue_specific_TFs$Tissue == "Pancreas"]
skin_TFs <- tissue_specific_TFs$Transcription.Factors..Ensembl.[tissue_specific_TFs$Tissue == "Skin (Suprapubic and Lower Leg)"]
blood_TFs <- tissue_specific_TFs$Transcription.Factors..Ensembl.[tissue_specific_TFs$Tissue == "Whole Blood"]
general_TFs <- tissue_specific_TFs$Transcription.Factors..Ensembl.[tissue_specific_TFs$Tissue == "General"]

brain_TFs <- str_split(brain_TFs, ",")
blood_TFs <- str_split(blood_TFs, ",")
cardiac_TFs <- str_split(cardiac_TFs, ",")
lung_TFs <- str_split(lung_TFs, ",")
muscle_TFs <- str_split(muscle_TFs, ",")
pancreas_TFs <- str_split(pancreas_TFs, ",")
skin_TFs <- str_split(skin_TFs, ",")
general_TFs <- str_split(general_TFs, ",")

brain_TFs <- unlist(brain_TFs)
blood_TFs <- unlist(blood_TFs)
cardiac_TFs <- unlist(cardiac_TFs)
lung_TFs <- unlist(lung_TFs)
muscle_TFs <- unlist(muscle_TFs)
pancreas_TFs <- unlist(pancreas_TFs)
skin_TFs <- unlist(skin_TFs)
general_TFs <- unlist(general_TFs)
all_tissue_spec_TFs <- unique(c(brain_TFs, blood_TFs, cardiac_TFs, lung_TFs, muscle_TFs, 
                                pancreas_TFs, skin_TFs))

brain_centrality$genes <- gsub("\\..*","",brain_centrality$genes)
blood_centrality$genes <- gsub("\\..*","",blood_centrality$genes)
cardiac_centrality$genes <- gsub("\\..*","",cardiac_centrality$genes)
lung_centrality$genes <- gsub("\\..*","",lung_centrality$genes)
muscle_centrality$genes <- gsub("\\..*","",muscle_centrality$genes)
pancreas_centrality$genes <- gsub("\\..*","",pancreas_centrality$genes)
skin_centrality$genes <- gsub("\\..*","",skin_centrality$genes)
all_centrality$genes <- gsub("\\..*","",all_centrality$genes)
normal_centrality$genes <- gsub("\\..*","",normal_centrality$genes)


centrality_of_tissue_TFs <- list()
centrality_of_tissue_TFs[[1]] <- brain_centrality$degree[brain_centrality$genes %in% brain_TFs]
centrality_of_tissue_TFs[[2]] <- blood_centrality$degree[blood_centrality$genes %in% blood_TFs]
centrality_of_tissue_TFs[[3]] <- cardiac_centrality$degree[cardiac_centrality$genes %in% cardiac_TFs]
centrality_of_tissue_TFs[[4]] <- lung_centrality$degree[lung_centrality$genes %in% lung_TFs]
centrality_of_tissue_TFs[[5]] <- muscle_centrality$degree[muscle_centrality$genes %in% muscle_TFs]
centrality_of_tissue_TFs[[6]] <- pancreas_centrality$degree[pancreas_centrality$genes %in% pancreas_TFs]
centrality_of_tissue_TFs[[7]] <- skin_centrality$degree[skin_centrality$genes %in% skin_TFs]
centrality_of_tissue_TFs[[8]] <- all_centrality$degree[all_centrality$genes %in% all_tissue_spec_TFs]
centrality_of_tissue_TFs[[9]] <- normal_centrality$degree[normal_centrality$genes %in% all_tissue_spec_TFs]

names(centrality_of_tissue_TFs) <- c("Brain", "Blood", "Cardiac", "Lung", "Muscle", "Pancreas", "Skin", "Universal consensus", "Normal consensus")

centrality_of_tissue_TFs_df <- data.frame(unlist(centrality_of_tissue_TFs))
centrality_of_tissue_TFs_df$Tissue <- row.names(centrality_of_tissue_TFs_df)
row.names(centrality_of_tissue_TFs_df) <- c(1:dim(centrality_of_tissue_TFs_df)[1])
colnames(centrality_of_tissue_TFs_df) <- c("Degree", "Tissue")
centrality_of_tissue_TFs_df$Category <- "Tissue specific TFs"
centrality_of_tissue_TFs_df$Tissue[grep("Brain", centrality_of_tissue_TFs_df$Tissue)] <- "Brain"
centrality_of_tissue_TFs_df$Tissue[grep("Blood", centrality_of_tissue_TFs_df$Tissue)] <- "Blood"
centrality_of_tissue_TFs_df$Tissue[grep("Cardiac", centrality_of_tissue_TFs_df$Tissue)] <- "Cardiac"
centrality_of_tissue_TFs_df$Tissue[grep("Lung", centrality_of_tissue_TFs_df$Tissue)] <- "Lung"
centrality_of_tissue_TFs_df$Tissue[grep("Muscle", centrality_of_tissue_TFs_df$Tissue)] <- "Muscle"
centrality_of_tissue_TFs_df$Tissue[grep("Pancreas", centrality_of_tissue_TFs_df$Tissue)] <- "Pancreas"
centrality_of_tissue_TFs_df$Tissue[grep("Skin", centrality_of_tissue_TFs_df$Tissue)] <- "Skin"
centrality_of_tissue_TFs_df$Tissue[grep("Normal consensus", centrality_of_tissue_TFs_df$Tissue)] <- "Non-cancerous consensus"
centrality_of_tissue_TFs_df$Tissue[grep("Universal consensus", centrality_of_tissue_TFs_df$Tissue)] <- "Universal consensus"

centrality_of_general_TFs <- list()
centrality_of_general_TFs[[1]] <- brain_centrality$degree[brain_centrality$genes %in% general_TFs]
centrality_of_general_TFs[[2]] <- blood_centrality$degree[blood_centrality$genes %in% general_TFs]
centrality_of_general_TFs[[3]] <- cardiac_centrality$degree[cardiac_centrality$genes %in% general_TFs]
centrality_of_general_TFs[[4]] <- lung_centrality$degree[lung_centrality$genes %in% general_TFs]
centrality_of_general_TFs[[5]] <- muscle_centrality$degree[muscle_centrality$genes %in% general_TFs]
centrality_of_general_TFs[[6]] <- pancreas_centrality$degree[pancreas_centrality$genes %in% general_TFs]
centrality_of_general_TFs[[7]] <- skin_centrality$degree[skin_centrality$genes %in% general_TFs]
centrality_of_general_TFs[[8]] <- all_centrality$degree[all_centrality$genes %in% general_TFs]
centrality_of_general_TFs[[9]] <- normal_centrality$degree[normal_centrality$genes %in% general_TFs]

names(centrality_of_general_TFs) <- c("Brain", "Blood", "Cardiac", "Lung", "Muscle", "Pancreas", "Skin", "Universal consensus", "Non-cancerous consensus")
centrality_of_general_TFs_df <- data.frame(unlist(centrality_of_general_TFs))
centrality_of_general_TFs_df$Tissue <- row.names(centrality_of_general_TFs_df)
row.names(centrality_of_general_TFs_df) <- c(1:dim(centrality_of_general_TFs_df)[1])
colnames(centrality_of_general_TFs_df) <- c("Degree", "Tissue")
centrality_of_general_TFs_df$Category <- "General TFs"
centrality_of_general_TFs_df$Tissue[grep("Brain", centrality_of_general_TFs_df$Tissue)] <- "Brain"
centrality_of_general_TFs_df$Tissue[grep("Blood", centrality_of_general_TFs_df$Tissue)] <- "Blood"
centrality_of_general_TFs_df$Tissue[grep("Cardiac", centrality_of_general_TFs_df$Tissue)] <- "Cardiac"
centrality_of_general_TFs_df$Tissue[grep("Lung", centrality_of_general_TFs_df$Tissue)] <- "Lung"
centrality_of_general_TFs_df$Tissue[grep("Muscle", centrality_of_general_TFs_df$Tissue)] <- "Muscle"
centrality_of_general_TFs_df$Tissue[grep("Pancreas", centrality_of_general_TFs_df$Tissue)] <- "Pancreas"
centrality_of_general_TFs_df$Tissue[grep("Skin", centrality_of_general_TFs_df$Tissue)] <- "Skin"
centrality_of_general_TFs_df$Tissue[grep("Universal consensus", centrality_of_general_TFs_df$Tissue)] <- "Universal consensus"
centrality_of_general_TFs_df$Tissue[grep("Non-cancerous consensus", centrality_of_general_TFs_df$Tissue)] <- "Non-cancerous consensus"

set.seed(0)
non_TFs_select <- sample(non_TFs, 100)
centrality_of_non_TFs <- list()
centrality_of_non_TFs[[1]] <- brain_centrality$degree[brain_centrality$genes %in% non_TFs_select]
centrality_of_non_TFs[[2]] <- blood_centrality$degree[blood_centrality$genes %in% non_TFs_select]
centrality_of_non_TFs[[3]] <- cardiac_centrality$degree[cardiac_centrality$genes %in% non_TFs_select]
centrality_of_non_TFs[[4]] <- lung_centrality$degree[lung_centrality$genes %in% non_TFs_select]
centrality_of_non_TFs[[5]] <- muscle_centrality$degree[muscle_centrality$genes %in% non_TFs_select]
centrality_of_non_TFs[[6]] <- pancreas_centrality$degree[pancreas_centrality$genes %in% non_TFs_select]
centrality_of_non_TFs[[7]] <- skin_centrality$degree[skin_centrality$genes %in% non_TFs_select]
centrality_of_non_TFs[[8]] <- all_centrality$degree[all_centrality$genes %in% non_TFs_select]
centrality_of_non_TFs[[9]] <- normal_centrality$degree[normal_centrality$genes %in% non_TFs_select]

names(centrality_of_non_TFs) <- c("Brain", "Blood", "Cardiac", "Lung", "Muscle", "Pancreas", "Skin", "Universal consensus", "Non-cancerous consensus")
centrality_of_non_TFs_df <- data.frame(unlist(centrality_of_non_TFs))
centrality_of_non_TFs_df$Tissue <- row.names(centrality_of_non_TFs_df)
row.names(centrality_of_non_TFs_df) <- c(1:dim(centrality_of_non_TFs_df)[1])
colnames(centrality_of_non_TFs_df) <- c("Degree", "Tissue")
centrality_of_non_TFs_df$Category <- "Not a TF"
centrality_of_non_TFs_df$Tissue[grep("Brain", centrality_of_non_TFs_df$Tissue)] <- "Brain"
centrality_of_non_TFs_df$Tissue[grep("Blood", centrality_of_non_TFs_df$Tissue)] <- "Blood"
centrality_of_non_TFs_df$Tissue[grep("Cardiac", centrality_of_non_TFs_df$Tissue)] <- "Cardiac"
centrality_of_non_TFs_df$Tissue[grep("Lung", centrality_of_non_TFs_df$Tissue)] <- "Lung"
centrality_of_non_TFs_df$Tissue[grep("Muscle", centrality_of_non_TFs_df$Tissue)] <- "Muscle"
centrality_of_non_TFs_df$Tissue[grep("Pancreas", centrality_of_non_TFs_df$Tissue)] <- "Pancreas"
centrality_of_non_TFs_df$Tissue[grep("Skin", centrality_of_non_TFs_df$Tissue)] <- "Skin"
centrality_of_non_TFs_df$Tissue[grep("Universal consensus", centrality_of_non_TFs_df$Tissue)] <- "Universal consensus"
centrality_of_non_TFs_df$Tissue[grep("Non-cancerous consensus", centrality_of_non_TFs_df$Tissue)] <- "Non-cancerous consensus"

TFs_res <- rbind(centrality_of_tissue_TFs_df, centrality_of_general_TFs_df, centrality_of_non_TFs_df)
TFs_res$Category <- factor(TFs_res$Category, levels = c("Tissue specific TFs","General TFs","Not a TF"))
TFs_res$Tissue[TFs_res$Tissue == "Non-cancerous consensus"] <- "Non-cancerous\nconsensus"
TFs_res$Tissue[TFs_res$Tissue == "Universal consensus"] <- "Universal\nconsensus"
TFs_res$Tissue[TFs_res$Tissue == "Blood"] <- "Blood\nconsensus"
TFs_res$Tissue[TFs_res$Tissue == "Lung"] <- "Lung\nconsensus"
TFs_res$Tissue[TFs_res$Tissue == "Skin"] <- "Skin\nconsensus"
TFs_res$Tissue[TFs_res$Tissue == "Pancreas"] <- "Pancreas\nconsensus"
TFs_res$Tissue[TFs_res$Tissue == "Cardiac"] <- "Cardiac\nconsensus"
TFs_res$Tissue[TFs_res$Tissue == "Muscle"] <- "Muscle\nconsensus"
TFs_res$Tissue[TFs_res$Tissue == "Brain"] <- "CNS\nconsensus"
TFs_res$Tissue <- factor(TFs_res$Tissue, levels = c("Blood\nconsensus", "Lung\nconsensus", "Skin\nconsensus", "Pancreas\nconsensus", "Cardiac\nconsensus", "Muscle\nconsensus", "CNS\nconsensus", "Universal\nconsensus", "Non-cancerous\nconsensus"))

saveRDS(TFs_res, paste0(resDir, "TF_enrichment.rds"))


TFs_res$Category <- as.character(TFs_res$Category)
TFs_res$Category <- factor(TFs_res$Category, levels = c("Not a TF", "General TFs", "Tissue specific TFs"))

TFs_res$Tissue <- as.character(TFs_res$Tissue)
TFs_res$Tissue[grep("CNS", TFs_res$Tissue)] <- "CNS"
TFs_res$Tissue[grep("Blood", TFs_res$Tissue)] <- "Blood"
TFs_res$Tissue[grep("Lung", TFs_res$Tissue)] <- "Lung"
TFs_res$Tissue[grep("Cardiac", TFs_res$Tissue)] <- "Cardiac"
TFs_res$Tissue[grep("Pancreas", TFs_res$Tissue)] <- "Pancreas"
TFs_res$Tissue[grep("Skin", TFs_res$Tissue)] <- "Skin"
TFs_res$Tissue[grep("Muscle", TFs_res$Tissue)] <- "Skeletal Muscle"
TFs_res$Tissue[grep("Non-cancerous", TFs_res$Tissue)] <- "Non-cancer\nconsensus"
TFs_res$Tissue <- factor(TFs_res$Tissue, levels = c("Blood", "Cardiac", "CNS", "Lung", 
                                                    "Pancreas", "Skeletal Muscle", "Skin", "Non-cancer\nconsensus", "Universal\nconsensus"))

stat.test <- TFs_res %>%
  group_by(Tissue) %>%
  wilcox_test(Degree ~ Category)
stat.test <- stat.test %>%
  add_xy_position(x = "Tissue", dodge = 0.8)
stat.test$p.adj.signif[stat.test$p.adj.signif == "****"] <- "p < 1e-4"
stat.test$p.adj.signif[stat.test$p.adj.signif == "***"] <- "p < 0.001"
stat.test$p.adj.signif[stat.test$p.adj.signif == "**"] <- "p < 0.01"
stat.test$p.adj.signif[stat.test$p.adj.signif == "*"] <- "p < 0.1"

ggplot(TFs_res, aes(x = Tissue, y = Degree)) + 
  geom_boxplot(aes(fill = Category), outlier.colour = "darkgrey", outlier.size =  0.5) + theme_classic() +
  ggtitle("Tissue-specific network\nnode degree distributions") + coord_flip() +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01,
                     bracket.nudge.y = -1, hide.ns = TRUE, size = 3.5, fontface = "bold") + 
  xlab("") + 
  theme(legend.title = element_text(size = 12, face = "bold"), plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

