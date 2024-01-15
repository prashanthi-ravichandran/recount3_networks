# Description
# This script does prepares uncorrected data used in inference of the universal consensus network
# 1. For GTEx projects exclude those with no applicable study, or CML controls
# 2. For GTEx projects split the samples by tissue of origin and save them 
# 3. For SRA projects, obtain replicate merged data and exclude microRNA and single cell data
# 4. Log transforms RPKM gene expression
# 5. Saves the processed expression data

# Input: RPKM data with gene filters obtained from data/rpkm/gene_filters (for GTEx and TCGA)
# Input: Replicate merged data from SRA obtained from data/replicates_merged


rm(list = ls())
library(recount3)
library(dplyr)
library(tidyverse)
# Script to pre-process all consensus
homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
rpkmDir <- paste0(homeDir, "data/rpkm/gene_filters/")
dir.create(paste0(homeDir, "data/all_consensus/"))
dir.create(paste0(homeDir, "data/all_consensus/uncorrected_expression"))
uncorrectedDir <- paste0(homeDir, "data/all_consensus/uncorrected_expression/")

studies <- c("GTEx", "sra", "tcga")

# Start with GTEx
gtex.projects <- list.files(paste0(rpkmDir, studies[1], "/"))

gtex_list <- lapply(gtex.projects, function(iproj){
  readRDS(paste0(rpkmDir, studies[1], "/", iproj))
})


names(gtex_list) <- gsub(".rds", "", gtex.projects)
gtex_list[["STUDY_NA"]] <- NULL

lapply(gtex_list, function(rse){
  rse$gtex.smtsd <- gsub("\\(|\\)", "", rse$gtex.smtsd)
  rse$gtex.smtsd <- gsub(" - | ", "_", rse$gtex.smtsd)
  tissues <- unique(rse$gtex.smtsd)
  for(i in c(1:length(tissues))){
    print(tissues[i])
    rse_subset <- rse[ ,rse$gtex.smtsd == tissues[i]]
    print(dim(rse_subset))
    saveRDS(rse_subset, paste0(uncorrectedDir, studies[1], "/", tissues[i], ".rds"))
  }
})

# For SRA we need to apply microRNA and other filters to the data with replicates merged
replicatesMerged <- "/data/abattle4/prashanthi/recount3/data/replicates_merged/gene_filters/"
sra.projects <- list.files(replicatesMerged)
sra_list <- lapply(sra.projects, function(iproj){
  cat(iproj, sep="\n")
  readRDS(paste0(replicatesMerged, iproj))
})

names(sra_list) <- gsub(".rds", "", sra.projects)

sra_list <- lapply(sra_list, function(iexpr){
  cat(dim(iexpr), "\n")
  num_genes <- dim(iexpr)[1]
  counts <- assays(iexpr)$RPKM
  num_zero_expression <- colSums(counts <= 0)
  prop_zero_expression <- num_zero_expression/num_genes
  iexpr <- iexpr[ , prop_zero_expression < 0.5]
  cat(dim(iexpr), "\n")
  iexpr
})

# Remove any potential single cell samples 
n.sra <- lapply(sra_list, function(iexpr){
  dim(iexpr)[2]
})

n.sra <- unlist(n.sra)
sum(n.sra)

# Remove potential single cell samples
sra_list <- lapply(sra_list, function(iexpr){
  cat(dim(iexpr), "\n")
  iexpr[ ,colData(iexpr)$recount_pred.pred.type == "rna-seq" | is.na(colData(iexpr)$recount_pred.pred.type)]
})

n.sra <- lapply(sra_list, function(iexpr){
  dim(iexpr)[2]
})

n.sra <- unlist(n.sra)
sum(n.sra)

sra_list <- sra_list[n.sra >= 20]
sra_list[["SRP096986"]] <- NULL # This is a single cell study
sra_list[["SRP135684"]] <- NULL # this is a single cell study
sra_list[["SRP166966"]] <- NULL # this is a single cell study
sra_list[["SRP200058"]] <- NULL # this is a single cell study
sra_list[["SRP063998"]] <- NULL # this is a single cell study

tissue_df_annotations <- readRDS("/data/abattle4/prashanthi/recount3/data/tissue_df.rds")
sra_list <- lapply(sra_list, function(iexpr){
  cat(dim(iexpr), "\n")
  iexpr <- iexpr[ ,colnames(iexpr) %in% tissue_df_annotations$sample]
  cat(dim(iexpr), "\n")
  iexpr
})

n.sra <- lapply(sra_list, function(iexpr){
  dim(iexpr)[2]
})

n.sra <- unlist(n.sra)
sum(n.sra)
sra_list_filtered <- sra_list[n.sra > 0]

for(i in c(1:length(sra_list_filtered))){
  print(i)
  saveRDS(sra_list_filtered[[i]], paste0(uncorrectedDir, studies[2], "/", names(sra_list_filtered)[i], ".rds"))
}

# Process TCGA
tcga.projects <- list.files(paste0(rpkmDir, studies[3], "/"))
tcga_list <- lapply(tcga.projects, function(iproj){
  readRDS(paste0(rpkmDir, studies[3], "/", iproj))
})

names(tcga_list) <- gsub(".rds", "", tcga.projects)
tcga_list <- lapply(tcga_list, function(irse){
  irse[ ,!is.na(irse$tcga.cgc_sample_sample_type)]
})

tcga_list <- lapply(tcga_list, function(irse){
  irse[ ,!irse$tcga.cgc_sample_sample_type == ""]
})

tcga_list <- lapply(tcga_list, function(irse){
  irse[ ,!is.na(irse$tcga.xml_patient_id)]
})

tcga_nSamples <- lapply(tcga_list, function(irse){
  dim(irse)[2]
})

tcga_nInd <- lapply(tcga_list, function(irse){
  length(unique(irse$tcga.xml_patient_id))
})

tcga_nSamples <- unlist(tcga_nSamples)
tcga_nInd <- unlist(tcga_nInd)

for(i in c(1:length(tcga_list))){
  print(i)
  rse <- tcga_list[[i]]
  counts <- as.data.frame(t(rse@assays@data$RPKM))
  counts$patient <- paste0(rse$tcga.xml_patient_id, "_", rse$tcga.cgc_sample_sample_type)
  print(dim(rse))
  if(length(unique(counts$patient)) < dim(counts)[1]){
    counts <- counts %>% group_by(patient) %>% summarise_all(median) %>% as.data.frame
    rownames(counts) <- counts$patient
    counts$patient <- NULL
    counts <- t(counts)
    metaData <- as.data.frame(colData(rse))
    metaData <- metaData %>% group_by(tcga.xml_patient_id, tcga.cgc_sample_sample_type) %>% 
      dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE), 
                       across(where(is.character), dplyr::first),
                       across(where(is.logical), all)) %>% as.data.frame
    rse_merged <- SummarizedExperiment(assays = list(RPKM = counts), colData = metaData, rowData = rowData(rse))
  }else{
    rse_merged <- rse
  }
  print(dim(rse_merged))
  saveRDS(rse_merged, paste0(uncorrectedDir, studies[3], "/", names(tcga_list)[i], ".rds"))
}

