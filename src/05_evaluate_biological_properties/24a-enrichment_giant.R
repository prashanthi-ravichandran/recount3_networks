# Load required libraries
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
            "/home/pravich2/rlibs/4.0.2/gcc/9.3.0", 
            "/data/apps/extern/r-sva/3.42.0/rlib/R/library"))
library(Matrix)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(GOfuncR)
library(dplyr)
library(org.Hs.eg.db)

make_contingency <- function(test, bg, set){
  inset_intest <- length(intersect(test,set))
  inset_inbg <- length(intersect(bg,set))
  notinset_intest <- length(test) - inset_intest
  notinset_inbg <- length(bg) - inset_inbg
  cmat = matrix(c(inset_intest, notinset_intest, inset_inbg, notinset_inbg), nrow = 2)
  cmat
}

compute_enrichment <- function(genelist, background, pathway){
  cont.mat <- make_contingency(genelist , background , pathway)
  res <- fisher.test(cont.mat, conf.int = T)
  res
}

home.dir <- "/data/abattle4/prashanthi/recount3_paper/"
dat.dir <- paste0(home.dir, "data/")
res.dir <- paste0(home.dir, "results/weighted_cov_networks/")
plot.dir <- paste0(home.dir, "plots/")

compute_odds <- function(agg_level, nstudies, lambda, ref_name){
  net <- readRDS(paste0(res.dir, agg_level, "/net_", nstudies, "/lambda_", 
                        sprintf("%.2f" , lambda), ".rds"))
  geneData <-  readRDS(paste0(dat.dir, "geneData.rds"))
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
  zero_degree <- degree_df[degree_df$degree == 0, ]
  nonzero_degree <- degree_df[degree_df$degree > 0, ]
  zero_degree$decile <- 0
  nonzero_degree$decile <- ntile(nonzero_degree$degree, 4)
  degree_df <- rbind(zero_degree, nonzero_degree)
  degree_df$decile <- degree_df$decile + 1
  # Read and process reference
  reference <- read.table(paste0(dat.dir, "human_base/", ref_name))
  reference$V1 <- as.character(reference$V1)
  reference$V2 <- as.character(reference$V2)
  reference$protein1_symbol <- mapIds(org.Hs.eg.db, keys = reference$V1, column="SYMBOL",
                                      keytype="ENTREZID",multiVals="first")
  reference$protein2_symbol <- mapIds(org.Hs.eg.db, keys = reference$V2, column="SYMBOL",
                                      keytype="ENTREZID",multiVals="first")
  reference <- reference[reference$V3 == 1, ]
  reference <- reference[ ,c("protein1_symbol", "protein2_symbol")]
  reference <- reference[!is.na(reference$protein1_symbol), ]
  reference <- reference[!is.na(reference$protein2_symbol), ]
  reference <- reference[reference$protein1_symbol %in% rownames(net), ]
  reference <- reference[reference$protein2_symbol %in% rownames(net), ]
  
  reference <- as.matrix(reference)
  reference_graph <- graph_from_edgelist(reference)
  reference_degree <- degree(reference_graph)
  reference_degree <- data.frame("genes" = names(reference_degree), 
                                 "degree" = as.numeric(reference_degree))
  central_genes <- reference_degree[reference_degree$degree >= quantile(reference_degree$degree, 0.8), ]
  
  cat(length(central_genes$genes), "\n")
  cat(length(intersect(central_genes$genes, degree_df$gene)))
  if(length(intersect(central_genes$genes, degree_df$gene)) < 100){
    stop("Insufficient overlap between reference and networks")
  }
  # Compute enrichment
  odds_ratio <- c()
  p_value <- c()
  odds_ratio_lower <- c()
  odds_ratio_upper <- c()
  for(idegree in c(1:max(degree_df$decile))){
    test <- degree_df$gene[degree_df$decile >= idegree]
    enrich_res <- compute_enrichment(test, rownames(net), central_genes$genes)
    odds_ratio <- c(odds_ratio, enrich_res$estimate)
    p_value <- c(p_value, enrich_res$p.value)
    odds_ratio_lower <- c(odds_ratio_lower, enrich_res$conf.int[1])
    odds_ratio_upper <- c(odds_ratio_upper, enrich_res$conf.int[2])
  }
  res_odds <- data.frame(c(1:max(degree_df$decile)), odds_ratio, p_value, 
                         odds_ratio_lower, odds_ratio_upper)
  colnames(res_odds) <- c("degree", "odds", "pvalue", "odds_ratio_lower", "odds_ratio_upper")
  res_odds
}


######## Blood
blood_net_blood_ref <- compute_odds("blood/all", 65, 0.24, "blood.dat")
consensus_net_blood_ref <- compute_odds("all_consensus", 966, 0.18, "blood.dat")
blood_GTEx_blood_ref <- compute_odds("blood/GTEx", 1, 0.28, "blood.dat")
fibroblast_net_blood_ref <- compute_odds("fibroblasts/all", 13, 0.30, "blood.dat")

blood_net_blood_ref$net <- "Blood"
consensus_net_blood_ref$net <- "Universal Consensus"
blood_GTEx_blood_ref$net <- "Blood (GTEx)"
fibroblast_net_blood_ref$net <- "Fibroblasts"
blood_res <- rbind(blood_net_blood_ref, blood_GTEx_blood_ref, consensus_net_blood_ref, fibroblast_net_blood_ref)


p1 <- ggplot(blood_res, aes(x = degree, y = odds, fill = net)) + geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 1, color = "black", linetype = 2) + ggtitle("Blood") + geom_errorbar(mapping = aes(x = degree, ymin = odds_ratio_lower, ymax = odds_ratio_upper, fill = net), width = 0.2, position=position_dodge(0.9)) + theme_classic()+ theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + ylab("Odds ratio") + xlab("Degree quintile")


######## CNS
cns_net_cns_ref <- compute_odds("central_nervous_system/all", 53, 0.30, "brain.dat")
consensus_net_cns_ref <- compute_odds("all_consensus", 966, 0.18, "brain.dat")
cns_GTEx_cns_ref <- compute_odds("central_nervous_system/GTEx", 13, 0.32, "brain.dat")
blood_net_cns_ref <- compute_odds("blood/all", 65, 0.24, "brain.dat")

cns_net_cns_ref$net <- "Central nervous system"
consensus_net_cns_ref$net <- "Universal consensus"
cns_GTEx_cns_ref$net <- "Central nervous system (GTEx)"
blood_net_cns_ref$net <- "Blood"

cns_res <- rbind(cns_net_cns_ref, consensus_net_cns_ref, cns_GTEx_cns_ref, blood_net_cns_ref)

cns_res$net <- factor(cns_res$net, levels = c("Central nervous system", "Central nervous system (GTEx)", "Blood", "Universal consensus"))
p2 <- ggplot(cns_res, aes(x = degree, y = odds, fill = net)) + geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 1, color = "black", linetype = 2) + theme_classic() + ggtitle("Central nervous system") + geom_errorbar(mapping = aes(x = degree, ymin = odds_ratio_lower, ymax = odds_ratio_upper, fill = net), width = 0.2, position=position_dodge(0.9)) + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + ylab("Odds ratio") + xlab("Degree quintile") + guides(fill=guide_legend(nrow=2,byrow=TRUE))

######## Skin
skin_net_skin_ref <- compute_odds("skin/all", 20, 0.26, "skin.dat")
consensus_net_skin_ref <- compute_odds("all_consensus", 966, 0.18, "skin.dat")
skin_GTEx_skin_ref <- compute_odds("skin/GTEx", 2, 0.28, "skin.dat")
blood_net_skin_ref <- compute_odds("blood/all", 65, 0.24, "skin.dat")

skin_net_skin_ref$net <- "Skin"
consensus_net_skin_ref$net <- "Universal consensus"
skin_GTEx_skin_ref$net <- "Skin (GTEx)"
blood_net_skin_ref$net <- "Blood"

skin_res <- rbind(skin_net_skin_ref, consensus_net_skin_ref, skin_GTEx_skin_ref, blood_net_skin_ref)

skin_res$net <- factor(skin_res$net, levels = c("Skin", "Skin (GTEx)", "Blood", "Universal consensus"))

p3 <- ggplot(skin_res, aes(x = degree, y = odds, fill = net)) + geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 1, color = "black", linetype = 2) + theme_classic() + ggtitle("Skin") + geom_errorbar(mapping = aes(x = degree, ymin = odds_ratio_lower, ymax = odds_ratio_upper, fill = net), width = 0.2, position=position_dodge(0.9)) + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + ylab("Odds ratio") + xlab("Degree quintile")

###### Liver
liver_net_liver_ref <- compute_odds("liver/all", 28, 0.26, "liver.dat")
consensus_net_liver_ref <- compute_odds("all_consensus", 966, 0.18, "liver.dat")
liver_GTEx_liver_ref <- compute_odds("liver/GTEx", 1, 0.40, "liver.dat")
blood_net_liver_ref <- compute_odds("blood/all", 65, 0.24, "liver.dat")
liver_net_liver_ref$net <- "Liver"
consensus_net_liver_ref$net <- "Universal consensus"
liver_GTEx_liver_ref$net <- "Liver (GTEx)"
blood_net_liver_ref$net <- "Blood"
liver_res <- rbind(liver_net_liver_ref, consensus_net_liver_ref, liver_GTEx_liver_ref, blood_net_liver_ref)

liver_res$net <- factor(liver_res$net, levels = c("Liver", "Liver (GTEx)", "Blood", "Universal consensus"))

p4 <- ggplot(liver_res, aes(x = degree, y = odds, fill = net)) + geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 1, color = "black", linetype = 2) + theme_classic() + ggtitle("Liver") + geom_errorbar(mapping = aes(x = degree, ymin = odds_ratio_lower, ymax = odds_ratio_upper, fill = net), width = 0.2, position=position_dodge(0.9)) + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + ylab("Odds ratio") + xlab("Degree quintile")

### Lung
lung_net_lung_ref <- compute_odds("lung/all", 10, 0.30, "lung.dat")
consensus_net_lung_ref <- compute_odds("all_consensus", 966, 0.18, "lung.dat")
lung_GTEx_lung_ref <- compute_odds("lung/GTEx", 1, 0.30, "lung.dat")
blood_net_lung_ref <- compute_odds("blood/all", 65, 0.24, "lung.dat")

lung_net_lung_ref$net <- "Lung"
consensus_net_lung_ref$net <- "Universal consensus"
lung_GTEx_lung_ref$net <- "Lung (GTEx)"
blood_net_lung_ref$net <- "Blood"

lung_res <- rbind(lung_net_lung_ref, consensus_net_lung_ref, lung_GTEx_lung_ref, blood_net_lung_ref)

lung_res$net <- factor(lung_res$net, levels = c("Lung", "Lung (GTEx)", "Blood", "Universal consensus"))

p5 <- ggplot(lung_res, aes(x = degree, y = odds, fill = net)) + geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 1, color = "black", linetype = 2) + theme_classic() + ggtitle("Lung") + geom_errorbar(mapping = aes(x = degree, ymin = odds_ratio_lower, ymax = odds_ratio_upper, fill = net), width = 0.2, position=position_dodge(0.9)) + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + ylab("Odds ratio") + xlab("Degree quintile")

##### Adipose
adipose_net_adipose_ref <- compute_odds("adipose/all", 11, 0.26, "adipose_tissue.dat")
consensus_net_adipose_ref <- compute_odds("all_consensus", 966, 0.18, "adipose_tissue.dat")
adipose_GTEx_adipose_ref <- compute_odds("adipose/GTEx", 2, 0.28, "adipose_tissue.dat")
blood_net_adipose_ref <- compute_odds("blood/all", 65, 0.24, "adipose_tissue.dat")

adipose_net_adipose_ref$net <- "Adipose"
consensus_net_adipose_ref$net <- "Universal consensus"
adipose_GTEx_adipose_ref$net <- "Adipose (GTEx)"
blood_net_adipose_ref$net <- "Blood"

adipose_res <- rbind(adipose_net_adipose_ref, consensus_net_adipose_ref, adipose_GTEx_adipose_ref, blood_net_adipose_ref)

adipose_res$net <- factor(adipose_res$net, levels = c("Adipose", "Adipose (GTEx)", "Blood", "Universal consensus"))

p6 <- ggplot(adipose_res, aes(x = degree, y = odds, fill = net)) + geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept = 1, color = "black", linetype = 2) + theme_classic() + ggtitle("Adipose") + geom_errorbar(mapping = aes(x = degree, ymin = odds_ratio_lower, ymax = odds_ratio_upper, fill = net), width = 0.2, position=position_dodge(0.9)) + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + ylab("Odds ratio") + xlab("Degree quintile")

