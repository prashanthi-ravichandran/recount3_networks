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
    #cat(length(iedge), "\n")
    iedge
  })
  edge.unformatted
}

inputArgs <- commandArgs(TRUE)
tissue <- inputArgs[1]
agg_level <- inputArgs[2]
nStudies <- inputArgs[3]
#tissue <- "blood"
#agg_level <- "all"
#nStudies <- "65"
lambda <- seq(0.18, 0.40, 0.02)
lambda <- sprintf("%.2f", lambda)
homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
res.dir <- paste0(homeDir, "results/weighted_cov_networks/")
dat.dir <- paste0(homeDir, "data/")

geneData <- readRDS(paste0(dat.dir, "geneData.rds"))
# Read in networks
index <- 1
net <- list()
for(ilambda in lambda){
  print(index)
  print(ilambda)
  inet <- readRDS(paste0(res.dir, tissue, "/", agg_level, "/net_", nStudies, "/lambda_", ilambda, ".rds"))
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

edgeList <- edge.list(net)

# Read in reference
snap_pathways <- read.delim(paste0(dat.dir, "PPI/PPT-Ohmnet_tissues-combined.edgelist"))

select_tissues <- switch (tissue,
                          "adipose" = c("adipose_tissue"), 
                          "airway" = c("bronchial_epithelial_cell","bronchus","trachea"), 
                          "B_cells" = c("b_lymphocyte"), 
                          "blood" = c("blood_plasma", "blood", "blood_platelet"), 
                          "breast" = c("mammary_epithelium", "mammary_gland"), 
                          "cardiac" = c("cardiac_muscle", "heart"), 
                          "central_nervous_system" = c("amygdala", "brain", "caudate_nucleus", "caudate_putamen", "central_nervous_system",
                                                       "cerebellar_cortex", "cerebellum", "cerebral_cortex", "corpus_callosum", 
                                                       "dentate_gyrus", "diencephalon", "forebrain", "frontal_lobe", "hippocampus",
                                                       "hypothalamus",  "locus_ceruleus", "medulla_oblongata", "midbrain", 
                                                       "nucleus_accumbens" , "occipital_lobe", "occipital_pole", "parietal_lobe", "pons", 
                                                       "spinal_cord", "substantia_nigra", "subthalamic_nucleus" , "telencephalon" , "temporal_lobe", 
                                                       "thalamus"),
                          "colon" = c("large_intestine", "colon", "cecum"), 
                          "esophagus" = c("esophagus"), 
                          "eye" = c("lens", "eye", "retina", "tear_gland", "cornea", "choroid"), 
                          "fibroblasts" = c("skin_fibroblast"),
                          "intestine" = c("intestine", "duodenum", "ileum", "jejunum", "small_intestine"), 
                          "ipscs" = c("embryo"), 
                          "hescs" = c("embryo"), 
                          "kidney" = c("renal_glomerulus", "nephron", "renal_tubule", "kidney"), 
                          "liver" = c("liver", "hepatocyte"), 
                          "lung" = c("lung"), 
                          "multipotent_cells" = c("hematopoietic_stem_cell"), 
                          "myeloid_cells" = c("granulocyte", "dendritic_cell", "macrophage", "monocyte"), 
                          "nervous_system" = c("peripheral_nervous_system", "neuron", "nervous_system" ), 
                          "pancreas" = c("pancreas"), 
                          "pbmcs_t_cells" = c("culture_condition_cd8_cell", "lymphocyte", "megakaryocyte" , "mononuclear_phagocyte", 
                                              "natural_killer_cell", "neutrophil", "t_lymphocyte", "basophil"), 
                          "prostate" = c("prostate_gland" ), 
                          "skeletal_muscle" = c("muscle", "skeletal_muscle"), 
                          "skin" = c("skin_fibroblast","skin", "epidermis", "keratinocyte", "hair_follicle"), 
                          "stomach" = c("stomach"), 
                          "vascular" =c("vascular_endothelial_cell", "vascular_endothelium", "artery", "aorta", 
                                        "blood_vessel", "umbilical_vein_endothelial_cell")
)

snap_pathways <- snap_pathways[snap_pathways$tissue %in% select_tissues,  ]
snap_pathways$protein1 <- as.character(snap_pathways$protein1)
snap_pathways$protein2 <- as.character(snap_pathways$protein2)
snap_pathways$protein1_symbol <- mapIds(org.Hs.eg.db, keys = snap_pathways$protein1, column="SYMBOL",
                                        keytype="ENTREZID",multiVals="first")
snap_pathways$protein2_symbol <- mapIds(org.Hs.eg.db, keys = snap_pathways$protein2, column="SYMBOL",
                                        keytype="ENTREZID",multiVals="first")
snap_pathways <- snap_pathways[snap_pathways$protein1_symbol %in% rownames(net[[1]]), ]
snap_pathways <- snap_pathways[snap_pathways$protein2_symbol %in% rownames(net[[1]]), ]

true_tissue_pathways <-c()
for(i in c(1:dim(snap_pathways)[1])){
  genes <- c(snap_pathways$protein1_symbol[i], snap_pathways$protein2_symbol[i])
  genes <- sort(genes)
  true_tissue_pathways[i] <- paste(genes[1], genes[2], sep = "_")
}

true_tissue_pathways <- unique(true_tissue_pathways)

# Create full network
unique_genes <- unique(c(rownames(net[[1]]), snap_pathways$protein1_symbol,
                         snap_pathways$protein2_symbol))
unique_genes <- unique_genes[order(unique_genes)]
ngenes <- length(unique_genes)
all.net <- matrix(1, nrow = ngenes, ncol = ngenes)
colnames(all.net) <- unique_genes
rownames(all.net) <- unique_genes
all.net[lower.tri(all.net, diag = T)] <- NA
all.edges <- reshape2::melt(all.net, na.rm = T)
all.edgelist <- paste(all.edges[ ,1], all.edges[ ,2], sep = "_")

# Compute the actual odds ratio
compute_odds <- function(all_edges, net_edges, true_edges){
  contingency_table <- matrix(rep(0,4), nrow = 2, ncol = 2)
  rownames(contingency_table) <- c("In_network", "Out_network")
  colnames(contingency_table) <- c("In_pathway", "Out_pathway")
  contingency_table[1,1] <- length(intersect(net_edges, true_edges))
  contingency_table[1,2] <- length(net_edges) - length(intersect(net_edges, true_edges))
  contingency_table[2,1] <- length(true_edges) - length(intersect(net_edges, true_edges))
  contingency_table[2,2] <- length(intersect(all_edges[!all_edges %in% net_edges], all_edges[!all_edges %in% true_edges]))
  odds <- (contingency_table[1,1]/contingency_table[2,1])/(contingency_table[1,2]/contingency_table[2,2])
  odds
}

odds_actual <- lapply(edgeList, function(iedge){
  compute_odds(all.edgelist, iedge, true_tissue_pathways)
})

odds_actual <- unlist(odds_actual)

# Permuted odds
nPerm <- 100
odds_perm <- list()
for(p in c(1:nPerm)){
  cat(p, "\n")
  net_perm <- lapply(net, function(inet){
    nodes_perm <- sample(rownames(inet))
    rownames(inet) <- nodes_perm
    colnames(inet) <- nodes_perm
    inet
  })
  edgeList_perm <- edge.list(net_perm)
  odds <- lapply(edgeList_perm, function(iedge){
    compute_odds(all.edgelist, iedge, true_tissue_pathways)
  })
  odds <- unlist(odds)
  odds_perm[[p]] <- odds
}

odds_perm <- do.call(cbind, odds_perm)
colnames(odds_perm) <- paste0("perm_", c(1:nPerm))

res_odds <- list("lambda" = lambda, "odds_actual" = odds_actual, "odds_perm" = odds_perm)
saveRDS(res_odds, paste0(res.dir, tissue, "/", agg_level, "/net_", nStudies, "/res_odds_perm.rds"))
