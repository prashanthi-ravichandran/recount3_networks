# Description
# This script evaluates for a number of networks the enrichment of genes present in certain GO terms
# Among central network nodes
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
library(dplyr)


home.dir <- "/data/abattle4/prashanthi/recount3/"
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
  res <- fisher.test(cont.mat, conf.int = T)
  res
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
  p_value <- c()
  odds_ratio_lower <- c()
  odds_ratio_upper <- c()
  zero_degree <- degree_df[degree_df$degree == 0, ]
  nonzero_degree <- degree_df[degree_df$degree > 0, ]
  zero_degree$decile <- 0
  nonzero_degree$decile <- ntile(nonzero_degree$degree, 9)
  degree_df <- rbind(zero_degree, nonzero_degree)
  degree_df$decile <- degree_df$decile + 1
  
  for(idegree in c(1:max(degree_df$decile))){
    test <- degree_df$gene[degree_df$decile >= idegree]
    enrich_res <- kegg_enrichment(test, rownames(net), pathway$gene)
    odds_ratio <- c(odds_ratio, enrich_res$estimate)
    p_value <- c(p_value, enrich_res$p.value)
    odds_ratio_lower <- c(odds_ratio_lower, enrich_res$conf.int[1])
    odds_ratio_upper <- c(odds_ratio_upper, enrich_res$conf.int[2])
  }
  res_odds <- data.frame(c(1:max(degree_df$decile)), odds_ratio)
  colnames(res_odds) <- c("degree", GO_term)
  
  res_pvalue <- data.frame(c(1:max(degree_df$decile)),p_value)
  colnames(res_pvalue) <- c("degree", GO_term)
  
  res_odds_lower <- data.frame(c(1:max(degree_df$decile)), odds_ratio_lower)
  colnames(res_odds_lower) <- c("degree", GO_term)
  
  res_odds_upper <- data.frame(c(1:max(degree_df$decile)),odds_ratio_upper)
  colnames(res_odds_upper) <- c("degree", GO_term)
  
  list("odds_ratio" = res_odds,
       "pvalue" = res_pvalue, 
       "odds_ratio_lower" = res_odds_lower,
       "odds_ratio_upper" = res_odds_upper)
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
                    "GO:0030168", 
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
  net_odds_pathway("all_consensus", 0.18, 966, GO_term)
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
  net_odds_pathway("all_consensus", 0.18, 966, GO_term)
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
  net_odds_pathway("all_consensus", 0.18, 966, GO_term)
})
normal_consensus_skin <- lapply(skin_pathways, function(GO_term){
  net_odds_pathway("normal_consensus", 0.18, 629, GO_term)
})
skin_GTEx_skin <- lapply(skin_pathways, function(GO_term){
  net_odds_pathway("skin/GTEx", 0.28, 2, GO_term)
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
  net_odds_pathway("all_consensus", 0.18, 966, GO_term)
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
  net_odds_pathway("all_consensus", 0.18, 966, GO_term)
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
  net_odds_pathway("all_consensus", 0.18, 966, GO_term)
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


convert_list_to_df <- function(ilist, net){
  odds <- ilist[[1]][["odds_ratio"]]
  odds_lower <- ilist[[1]][["odds_ratio_lower"]]
  odds_upper <- ilist[[1]][["odds_ratio_upper"]]
  pvalue <- ilist[[1]][["pvalue"]]
  for(i in c(2:length(ilist))){
    odds <- odds %>% left_join(ilist[[i]][["odds_ratio"]], by = "degree")
    odds_lower <- odds_lower %>% left_join(ilist[[i]][["odds_ratio_lower"]], by = "degree")
    odds_upper <- odds_upper %>% left_join(ilist[[i]][["odds_ratio_upper"]], by = "degree")
    pvalue <- pvalue %>% left_join(ilist[[i]][["pvalue"]], by = "degree")
  }
  odds <- reshape2::melt(odds, id.vars = "degree")
  odds_lower <- reshape2::melt(odds_lower, id.vars = "degree")
  odds_upper <- reshape2::melt(odds_upper, id.vars = "degree")
  pvalue <- reshape2::melt(pvalue, id.vars = "degree")
  
  colnames(odds) <- c("degree", "pathway", "odds_ratio")
  colnames(odds_lower) <- c("degree", "pathway", "odds_ratio_lower")
  colnames(odds_upper) <- c("degree", "pathway", "odds_ratio_upper")
  colnames(pvalue) <- c("degree", "pathway", "pvalue")
  
  res <- odds %>% left_join(odds_lower, by = c("degree", "pathway")) %>%
         left_join(odds_upper, by = c("degree", "pathway")) %>%
         left_join(pvalue, by = c("degree", "pathway"))
  res$net <- net
  res
}

all_consensus_res <- convert_list_to_df(all_consensus_res, "Universal consensus")
normal_consensus_res <- convert_list_to_df(normal_consensus_res, "Non-cancerous consensus")
cancer_consensus_res <- convert_list_to_df(cancer_consensus_res, "Cancerous consensus")
cns_consensus_res <- convert_list_to_df(cns_consensus_res, "CNS")
blood_consensus_res <- convert_list_to_df(blood_consensus_res, "Blood")

all_consensus_neural <- convert_list_to_df(all_consensus_neural, "Universal consensus")
normal_consensus_neural <- convert_list_to_df(normal_consensus_neural, "Non-cancerous consensus")
cns_GTEx_neural <- convert_list_to_df(cns_GTEx_neural, "CNS (GTEx)")
cns_consensus_neural <- convert_list_to_df(cns_consensus_neural, "CNS")
blood_consensus_neural <- convert_list_to_df(blood_consensus_neural, "Blood")

all_consensus_blood <- convert_list_to_df(all_consensus_blood, "Universal consensus")
normal_consensus_blood <- convert_list_to_df(normal_consensus_blood, "Non-cancerous consensus")
blood_GTEx_blood <- convert_list_to_df(blood_GTEx_blood, "Blood (GTEx)")
blood_consensus_blood <- convert_list_to_df(blood_consensus_blood, "Blood")
cns_consensus_blood <- convert_list_to_df(cns_consensus_blood, "CNS")

all_consensus_skin <- convert_list_to_df(all_consensus_skin, "Universal consensus")
normal_consensus_skin <- convert_list_to_df(normal_consensus_skin, "Non-cancerous consensus")
skin_GTEx_skin <- convert_list_to_df(skin_GTEx_skin, "Skin (GTEx)")
skin_consensus_skin <- convert_list_to_df(skin_consensus_skin, "Skin")
blood_consensus_skin <- convert_list_to_df(blood_consensus_skin, "Blood")
cns_consensus_skin <- convert_list_to_df(cns_consensus_skin, "CNS")

all_consensus_adipose <- convert_list_to_df(all_consensus_adipose, "Universal consensus")
normal_consensus_adipose <- convert_list_to_df(normal_consensus_adipose, "Non-cancerous consensus")
adipose_GTEx_adipose <- convert_list_to_df(adipose_GTEx_adipose, "Adipose (GTEx)")
adipose_consensus_adipose <- convert_list_to_df(adipose_consensus_adipose, "Adipose")
blood_consensus_adipose <- convert_list_to_df(blood_consensus_adipose, "Blood")
cns_consensus_adipose <- convert_list_to_df(cns_consensus_adipose, "CNS")

all_consensus_liver <- convert_list_to_df(all_consensus_liver, "Universal consensus")
normal_consensus_liver <- convert_list_to_df(normal_consensus_liver, "Non-cancerous consensus")
liver_GTEx_liver <- convert_list_to_df(liver_GTEx_liver, "Liver (GTEx)")
liver_consensus_liver <- convert_list_to_df(liver_consensus_liver, "Liver")
blood_consensus_liver <- convert_list_to_df(blood_consensus_liver, "Blood")
cns_consensus_liver <- convert_list_to_df(cns_consensus_liver, "CNS")

all_consensus_lung <- convert_list_to_df(all_consensus_lung, "Universal consensus")
normal_consensus_lung <- convert_list_to_df(normal_consensus_lung, "Non-cancerous consensus")
lung_GTEx_lung <- convert_list_to_df(lung_GTEx_lung, "Lung (GTEx")
lung_consensus_lung <- convert_list_to_df(lung_consensus_lung, "Lung")
blood_consensus_lung <- convert_list_to_df(blood_consensus_lung, "Blood")
cns_consensus_lung <- convert_list_to_df(cns_consensus_lung, "CNS")

consensus_res <- rbind(all_consensus_res, normal_consensus_res, cancer_consensus_res ,cns_consensus_res, blood_consensus_res)
consensus_res$net <- as.factor(consensus_res$net)

neural_res <- rbind(all_consensus_neural, normal_consensus_neural, cns_GTEx_neural, cns_consensus_neural, blood_consensus_neural)
neural_res$net <- as.factor(neural_res$net)

blood_res <- rbind(all_consensus_blood, normal_consensus_blood, blood_GTEx_blood, blood_consensus_blood, cns_consensus_blood)
blood_res$net <- as.factor(blood_res$net)

skin_res <- rbind(all_consensus_skin, normal_consensus_skin, skin_GTEx_skin, skin_consensus_skin, cns_consensus_skin, blood_consensus_skin)
skin_res$net <- as.factor(skin_res$net)

adipose_res <- rbind(all_consensus_adipose, normal_consensus_adipose, adipose_GTEx_adipose, adipose_consensus_adipose, cns_consensus_adipose, blood_consensus_adipose)
adipose_res$net <- as.factor(adipose_res$net)

lung_res <- rbind(all_consensus_lung, normal_consensus_lung, lung_GTEx_lung, lung_consensus_lung, cns_consensus_lung, blood_consensus_lung)
lung_res$net <- as.factor(lung_res$net)

liver_res <- rbind(all_consensus_liver, normal_consensus_liver, liver_GTEx_liver, liver_consensus_liver, cns_consensus_liver, blood_consensus_liver)
liver_res$net <- as.factor(liver_res$net)


saveRDS(consensus_res, paste0(home.dir, "results/enrichment_plots/consensus_res.rds"))
saveRDS(neural_res, paste0(home.dir, "results/enrichment_plots/neural_res.rds"))
saveRDS(skin_res, paste0(home.dir, "results/enrichment_plots/skin_res.rds"))
saveRDS(blood_res, paste0(home.dir, "results/enrichment_plots/blood_res.rds"))
saveRDS(adipose_res, paste0(home.dir, "results/enrichment_plots/adipose_res.rds"))
saveRDS(liver_res, paste0(home.dir, "results/enrichment_plots/liver_res.rds"))
saveRDS(lung_res, paste0(home.dir, "results/enrichment_plots/lung_res.rds"))

pathways <- c(consensus_pathways, neural_pathways, skin_pathways, blood_pathways, adipose_pathways, liver_pathways, lung_pathways)
pathways_df <- get_names(pathways)
list_genes <- get_anno_genes(pathways)
ngenes <- c()
for(i in c(1:dim(pathways_df)[1])){
  ngenes[i] <- length(unique(list_genes$gene[list_genes$go_id == pathways_df$go_id[i]]))
}
pathways_df$ngenes <- ngenes
saveRDS(pathways_df, paste0(home.dir, "results/enrichment_plots/pathways.rds"))

consensus_res$net <- as.character(consensus_res$net)
consensus_res$net <- factor(consensus_res$net, levels = c("CNS", "Blood", "Cancerous consensus", "Non-cancerous consensus", "Universal consensus"))
p1 <- ggplot(consensus_res[consensus_res$pathway == "GO:0000278", ], aes(x = degree, y = odds_ratio, fill = net, color = net))+
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)+ scale_fill_brewer(palette="Set1") + 
  scale_color_brewer(palette="Set1") + ggtitle("Mitotic cell cycle") +
  theme_classic() + geom_hline(yintercept = 1, color="black", linetype="dashed") + 
  theme(legend.title = element_blank(), plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + 
  geom_errorbar(aes(ymin=odds_ratio_lower, ymax=odds_ratio_upper), width=.2, position=position_dodge(.9)) + xlab("Degree decile") + ylab("Odds ratio") + 
  scale_x_continuous(breaks=seq(1, 10, 1))

p2 <- ggplot(consensus_res[consensus_res$pathway == "GO:0051276", ], aes(x = degree, y = odds_ratio, fill = net, color = net))+
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)+ scale_fill_brewer(palette="Set1") + 
  scale_color_brewer(palette="Set1") + ggtitle("Chromosome organization") +
  theme_classic() + geom_hline(yintercept = 1, color="black", linetype="dashed") + 
  theme(legend.title = element_blank(), plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + 
  geom_errorbar(aes(ymin=odds_ratio_lower, ymax=odds_ratio_upper), width=.2, position=position_dodge(.9)) + xlab("Degree decile") + ylab("Odds ratio") + 
  scale_x_continuous(breaks=seq(1, 10, 1))

p3 <- ggplot(consensus_res[consensus_res$pathway == "GO:0033043", ], aes(x = degree, y = odds_ratio, fill = net, color = net))+
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)+ scale_fill_brewer(palette="Set1") + 
  scale_color_brewer(palette="Set1") + ggtitle("Regulation of organelle organization") +
  theme_classic() + geom_hline(yintercept = 1, color="black", linetype="dashed") + 
  theme(legend.title = element_blank(), plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + 
  geom_errorbar(aes(ymin=odds_ratio_lower, ymax=odds_ratio_upper), width=.2, position=position_dodge(.9)) + xlab("Degree decile") + ylab("Odds ratio") + 
  scale_x_continuous(breaks=seq(1, 10, 1))

p4 <- ggplot(consensus_res[consensus_res$pathway == "GO:0007017", ], aes(x = degree, y = odds_ratio, fill = net, color = net))+
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)+ scale_fill_brewer(palette="Set1") + 
  scale_color_brewer(palette="Set1") + ggtitle("Microtubule-based process") +
  theme_classic() + geom_hline(yintercept = 1, color="black", linetype="dashed") + 
  theme(legend.title = element_blank(), plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + 
  geom_errorbar(aes(ymin=odds_ratio_lower, ymax=odds_ratio_upper), width=.2, position=position_dodge(.9)) + xlab("Degree decile") + ylab("Odds ratio") + 
  scale_x_continuous(breaks=seq(1, 10, 1))

p5 <- ggplot(consensus_res[consensus_res$pathway == "GO:0061640", ], aes(x = degree, y = odds_ratio, fill = net, color = net))+
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)+ scale_fill_brewer(palette="Set1") + 
  scale_color_brewer(palette="Set1") + ggtitle("Cytoskeleton-dependent cytokinesis") +
  theme_classic() + geom_hline(yintercept = 1, color="black", linetype="dashed") + 
  theme(legend.title = element_blank(), plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + 
  geom_errorbar(aes(ymin=odds_ratio_lower, ymax=odds_ratio_upper), width=.2, position=position_dodge(.9)) + xlab("Degree decile") + ylab("Odds ratio") + 
  scale_x_continuous(breaks=seq(1, 10, 1))

p6 <- ggplot(consensus_res[consensus_res$pathway == "GO:0006974", ], aes(x = degree, y = odds_ratio, fill = net, color = net))+
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)+ scale_fill_brewer(palette="Set1") + 
  scale_color_brewer(palette="Set1") + ggtitle("Cellular response to DNA damage stimulus") +
  theme_classic() + geom_hline(yintercept = 1, color="black", linetype="dashed") + 
  theme(legend.title = element_blank(), plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + 
  geom_errorbar(aes(ymin=odds_ratio_lower, ymax=odds_ratio_upper), width=.2, position=position_dodge(.9)) + xlab("Degree decile") + ylab("Odds ratio") + 
  scale_x_continuous(breaks=seq(1, 10, 1))