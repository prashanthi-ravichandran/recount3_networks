# Description
# This script reads in the processed count and metadata and performs xCell deconvolution
# We then subset to the non-cancerous samples and compute t-SNE dimensionality reduction 
# which is used to identify tissue-context clusters and outliers
rm(list = ls())
# Author: Prashanthi Ravichandran
# Packages
library(irlba)
library(reshape2)
library(dplyr)
library(RSpectra)
library(scales)
library(rsvd)
library(stats)
library(ggplot2)
library(parallel)
library(matrixStats)
library(umap)
library(Rtsne)
library(forcats)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ruler)

# Required files
# Outputs from 08-compute_all_counts_all_metadata.R
# data/all_counts.rds 
# data/all_metadata.rds
# data/manual_annotations.tsv
# data/geneData.rds downloaded from the RSE objects obtained from Recount3
# xcell deconvolution functions from src/xcell_sort.R

# Output files
# 1. data/xcellScores.rds
# 2. data/tsne_normal_xcell.rds
# 3. tissue_xcell_tsne/ : A directory that contains the tSNE reductions for samples belonging to each of the 27 contexts

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
source(paste0(homeDir, "src/xcell_sort.R"))
metadata <- readRDS(paste0(homeDir, "data/all_metadata.rds"))

metadata$tissue.category[metadata$tissue.category == "Whole_Blood"] <- "blood"
metadata$tissue.category <- gsub("_", " ", metadata$tissue.category)
metadata$tissue.category <- tolower(metadata$tissue.category)
metadata$tissue.category[metadata$tissue.category == "brain cerebellum"] <- "cerebellum"
metadata$tissue.category[metadata$tissue.category == "skeletal muscle "] <- "skeletal muscle"
metadata$tissue.category[metadata$tissue.category == "muscle skeletal"] <- "skeletal muscle"
metadata$tissue.category[metadata$tissue.category == "ipsc"] <- "ipscs"
metadata$tissue.category[metadata$tissue.category == "prostate"] <- "prostrate"
metadata$source <- metadata$study
metadata$source[grep("DRP", metadata$source)] <- "sra"
metadata$source[grep("SRP", metadata$source)] <- "sra"
metadata$source[grep("ERP", metadata$source)] <- "sra"
metadata$source[metadata$source != "sra" & metadata$source != "GTEx"] <- "TCGA"

metadata$cancer[metadata$cancer == "cancer "] <- "cancer"
metadata$cancer[metadata$cancer == ""] <- "none"
metadata$cancer[metadata$tissue.category == "cells leukemia cell line cml"] <- "cancer"
metadata <- metadata[!metadata$tissue.category == "cells leukemia cell line cml", ]


metadata_cancer <- metadata[metadata$cancer == "cancer", ]
metadata_normal <- metadata[metadata$cancer == "none", ]
metadata_normal <- metadata_normal[!metadata_normal$tissue.category == "unknown", ]
metadata_cancer <- metadata_cancer[!metadata_cancer$tissue.category == "unknown", ]

counts <- readRDS(paste0(homeDir, "/data/all_counts.rds"))
counts <- t(counts)
gtex_normal <- counts[rownames(counts) %in% metadata_normal$sample[metadata_normal$source == "GTEx"], ]
gtex_cancer <- counts[rownames(counts) %in% metadata_cancer$sample[metadata_cancer$source == "GTEx"], ]
tcga_normal <- counts[rownames(counts) %in% metadata_normal$sample[metadata_normal$source == "TCGA"], ]
tcga_cancer <- counts[rownames(counts) %in% metadata_cancer$sample[metadata_cancer$source == "TCGA"], ]
sra_normal <- counts[rownames(counts) %in% metadata_normal$sample[metadata_normal$source == "sra"], ]
sra_cancer <- counts[rownames(counts) %in% metadata_cancer$sample[metadata_cancer$source == "sra"], ]

all_normal <- rbind(gtex_normal, tcga_normal, sra_normal)
all_cancer <- rbind(gtex_cancer, tcga_cancer, sra_cancer)

metadata_normal <- metadata_normal[match(rownames(all_normal), metadata_normal$sample), ]
metadata_cancer <- metadata_cancer[match(rownames(all_cancer), metadata_cancer$sample), ]

metadata_normal$tissue.category[metadata_normal$tissue.category == "smooth vessel"] <- "smooth muscle"
metadata_normal$tissue.category[metadata_normal$tissue.category == "trophoblast"] <- "trophoblasts"
metadata_normal$tissue.category[metadata_normal$tissue.category == "adipocyte"] <- "adipocytes"
tissue_palette <- read.csv(paste0(homeDir, "data/colors_clustering.csv"), header = FALSE)
plt.source.cols <- metadata_normal$source
plt.source.cols[plt.source.cols == "sra"] <- "green"
plt.source.cols[plt.source.cols == "GTEx"] <- "blue"
plt.source.cols[plt.source.cols == "TCGA"] <- "red"
plt.colors.tissues <- c()
tissue <- c()
for(i in c(1:dim(metadata_normal)[1])){
  #print(metadata_normal$tissue.category[i])
  plt.colors.tissues[i] <- tissue_palette$V2[tissue_palette$V1 == metadata_normal$tissue.category[i]]
  tissue[i] <- tissue_palette$V3[tissue_palette$V1 == metadata_normal$tissue.category[i]]
}
tissue_palette_category <- tissue_palette
tissue_palette_category$V1 <- NULL
tissue_palette_category <- distinct(tissue_palette_category, V3, .keep_all = TRUE)
rownames(tissue_palette_category) <- tissue_palette_category$V3

metadata_normal$colors <- plt.colors.tissues
metadata_normal$tissue <- tissue

plt.colors.category <- tissue_palette_category[metadata_normal$tissue, ]
plt.colors.category <- plt.colors.category$V2
plt.colors.category[metadata_normal$tissue == "Central nervous system"] <- "#EEEE00"

xcell_flag <- TRUE
if(xcell_flag){
  counts.xcellScore <- readRDS(paste0(homeDir,"data/xcellScores.rds"))
} else{
  source(paste0(homeDir, "src/xcell_sort.R"))
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  geneData <- readRDS(paste0(homeDir, "data/geneData.rds"))
  geneData <- geneData[geneData$gene_id %in% rownames(counts), ]
  geneData <- geneData[match(rownames(counts), geneData$gene_id), ]
  rownames(counts) <- geneData$gene_name
  counts.xcellScore <- xCellAnalysis(counts)
  saveRDS(counts.xcellScore, paste0(homeDir, "data/xcellScores.rds"))
}

normal_xcell_scores <- counts.xcellScore[ ,metadata_normal$sample]
normal_xcell_scores <- t(normal_xcell_scores)
set.seed(100)

tsne_computed <- 1
if(tsne_computed == 1){
  tsne_normal_xcell <- readRDS(paste0(homeDir, "data/tsne_normal_xcell.rds"))
}else{
  tsne_normal_xcell <- Rtsne(normal_xcell_scores, verbose = TRUE, check_duplicates = FALSE)
  saveRDS(tsne_normal_xcell, paste0(homeDir, "data/tsne_normal_xcell.rds"))
}

normal_xcell_tsne_df <- tsne_normal_xcell$Y
normal_xcell_tsne_df <- data.frame(normal_xcell_tsne_df)
colnames(normal_xcell_tsne_df) <- c("D1", "D2")

normal_xcell_tsne_df$tissue <- metadata_normal$tissue
normal_xcell_tsne_df$tissue_subtype <- metadata_normal$tissue.category
normal_xcell_tsne_df$col_tissue <- plt.colors.category
normal_xcell_tsne_df$col_subtype <- metadata_normal$colors
normal_xcell_tsne_df$study <- metadata_normal$study
normal_xcell_tsne_df$sample_id <- metadata_normal$sample
normal_xcell_tsne_df$tissue[normal_xcell_tsne_df$tissue == "salivary gland"] <- "Salivary gland"
normal_xcell_tsne_df$tissue[normal_xcell_tsne_df$tissue == "head and neck"] <- "Head and neck"
colnames(normal_xcell_tsne_df) <- c("D1", "D2", "Tissue", "Tissue subtype", "col_tissue", "col_subtype", "Study")
normal_xcell_tsne_df$Tissue[normal_xcell_tsne_df$Tissue == "Prostate"] <- "Prostrate"
normal_xcell_tsne_df$Tissue <- as.factor(normal_xcell_tsne_df$Tissue)
tissue_palette_category <- tissue_palette_category[order(tissue_palette_category$V3), ]

p1 <- ggplot(normal_xcell_tsne_df, aes(x = D1, y = D2, color = Tissue)) + 
  geom_point(alpha = 0.1) + 
  scale_colour_manual(values = tissue_palette_category$V2) + 
  theme_classic() + ggtitle("t-SNE Visualization of Cell Type Enrichment Scores\nof 63193 Non-cancerous Samples representing 48 Tissue Contexts") + 
  guides(colour = guide_legend(ncol = 3, override.aes = list(alpha = 1, size = 3))) + xlim(c(-50, 50)) + ylim(c(-50, 50)) + 
  theme(plot.title = element_text(face="bold", size = 14), axis.text = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"))
p1

blood_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Blood", ]
immune_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Immune system" , ]
CNS_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tssue == "Central nervous system" , ]
skin_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Skin" , ]
musculoskeletal_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Muscoloskeletal" , ]
adipose_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Adipose" , ]
liver_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Liver" , ]
nervous_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Nervous system" , ]
intestine_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Intestine" , ]
cardiac_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Cardiac" , ]
vascular_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Vascular" , ]
iPSCs_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "iPSCs" , ]
esophagus_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Esophagus" , ]
lung_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Lung" , ]
colon_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Colon" , ]
kidney_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Kidney" , ]
fibroblasts_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Fibroblasts" , ]
breast_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Breast" , ]
hESCs_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "hESCs" , ]
prostrate_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Prostrate" , ]
stomach_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Stomach" , ]
pancreas_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Pancreas" , ]
multipotent_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Multipotent cells" , ]
airway_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Airway" , ]
eye_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue == "Eye" , ]
skeletalMuscle_xcell_tsne <- normal_xcell_tsne_df[normal_xcell_tsne_df$Tissue_subtype == "skeletal muscle" , ]


saveDir <- paste0(homeDir, "data/tissue_xcell_tsne/")
saveRDS(blood_xcell_tsne, paste0(saveDir, "blood.rds"))
saveRDS(immune_xcell_tsne, paste0(saveDir, "immune.rds"))
saveRDS(CNS_xcell_tsne, paste0(saveDir, "CNS.rds"))
saveRDS(skin_xcell_tsne, paste0(saveDir, "skin.rds"))
saveRDS(musculoskeletal_xcell_tsne, paste0(saveDir, "musculoskeletal.rds"))
saveRDS(adipose_xcell_tsne, paste0(saveDir, "adipose.rds"))
saveRDS(liver_xcell_tsne, paste0(saveDir, "liver.rds"))
saveRDS(nervous_xcell_tsne, paste0(saveDir, "nervous.rds"))
saveRDS(intestine_xcell_tsne, paste0(saveDir, "intestine.rds"))
saveRDS(cardiac_xcell_tsne , paste0(saveDir, "cardiac.rds"))
saveRDS(vascular_xcell_tsne, paste0(saveDir, "vascular.rds"))
saveRDS(iPSCs_xcell_tsne , paste0(saveDir, "iPSCs.rds"))
saveRDS(esophagus_xcell_tsne , paste0(saveDir, "esophagus.rds"))
saveRDS(lung_xcell_tsne , paste0(saveDir, "lung.rds"))
saveRDS(colon_xcell_tsne , paste0(saveDir, "colon.rds")) 
saveRDS(kidney_xcell_tsne, paste0(saveDir, "kidney.rds"))
saveRDS(fibroblasts_xcell_tsne, paste0(saveDir, "fibroblasts.rds"))
saveRDS(breast_xcell_tsne , paste0(saveDir, "breast.rds"))
saveRDS(hESCs_xcell_tsne , paste0(saveDir, "hESCs.rds"))
saveRDS(prostrate_xcell_tsne , paste0(saveDir, "prostate.rds"))
saveRDS(stomach_xcell_tsne , paste0(saveDir, "stomach.rds"))
saveRDS(pancreas_xcell_tsne, paste0(saveDir, "pancreas.rds"))
saveRDS(multipotent_xcell_tsne, paste0(saveDir, "multipotent.rds"))
saveRDS(airway_xcell_tsne, paste0(saveDir, "airway.rds"))
saveRDS(eye_xcell_tsne , paste0(saveDir, "eye.rds"))
saveRDS(skeletalMuscle_xcell_tsne, paste0(saveDir, "skeletalMuscle.rds"))
