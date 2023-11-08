# Description
# This script reads data and organizes it into the individual
# tissue contexts that the samples belong to 
# excluding immune tissues

# Author: Prashanthi Ravichandran
rm(list = ls())
# Packages
library(stringr)

# Inputs
# 1. data/filter_two_fail.rds
# 2. data/all_consensus/GTEx/uncorrected_expression/
# 3. data/all_consensus/sra/uncorrected_expression/
# 4. data/all_consensus/tcga/uncorrected_expression/

# Outputs
# 1. *tissue*/all/uncorrected_expression
# 2. *tissue*/GTEx/uncorrected_expression


homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
filtered_tissue_specific <- readRDS(paste0(homeDir, "data/filter_two_fail.rds"))
filtered_tissue_specific$tissue <- tolower(filtered_tissue_specific$tissue)
filtered_tissue_specific$tissue <- gsub(" ", "_", filtered_tissue_specific$tissue)
filtered_tissue_specific$source <- filtered_tissue_specific$study
filtered_tissue_specific$source[grep("SRP", filtered_tissue_specific$source)] <- "sra"
filtered_tissue_specific$source[grep("ERP", filtered_tissue_specific$source)] <- "sra"
filtered_tissue_specific$source[grep("DRP", filtered_tissue_specific$source)] <- "sra"
filtered_tissue_specific$source[!filtered_tissue_specific$source %in% c("GTEx", "sra")] <- "tcga"

datDir <- paste0(homeDir, "data/")
for(i in unique(filtered_tissue_specific$tissue)){
  if(!dir.exists(paste0(datDir, i))){
    dir.create(paste0(datDir, i))
  }
  if(!dir.exists(paste0(datDir, i, "/all"))){
    dir.create(paste0(datDir, i, "/all"))
  }
  if(!dir.exists(paste0(datDir, i, "/GTEx"))){
    dir.create(paste0(datDir, i, "/GTEx"))
  }
  if(!dir.exists(paste0(datDir, i, "/all/uncorrected_expression"))){
    dir.create(paste0(datDir, i, "/all/uncorrected_expression"))
  }
  if(!dir.exists(paste0(datDir, i, "/all/corrected_expression"))){
    dir.create(paste0(datDir, i, "/all/corrected_expression"))
  }
  if(!dir.exists(paste0(datDir, i, "/all/single_covariances"))){
    dir.create(paste0(datDir, i, "/all/single_covariances"))
  }
  if(!dir.exists(paste0(datDir, i, "/all/weighted_covariances"))){
    dir.create(paste0(datDir, i, "/all/weighted_covariances"))
  }
  if(!dir.exists(paste0(datDir, i, "/GTEx/uncorrected_expression"))){
    dir.create(paste0(datDir, i, "/GTEx/uncorrected_expression"))
  }
  if(!dir.exists(paste0(datDir, i, "/GTEx/corrected_expression"))){
    dir.create(paste0(datDir, i, "/GTEx/corrected_expression"))
  }
  if(!dir.exists(paste0(datDir, i, "/GTEx/single_covariances"))){
    dir.create(paste0(datDir, i, "/GTEx/single_covariances"))
  }
  if(!dir.exists(paste0(datDir, i, "/GTEx/weighted_covariances"))){
    dir.create(paste0(datDir, i, "/GTEx/weighted_covariances"))
  }
}

for(i in unique(filtered_tissue_specific$tissue)){
  print(i)
  df <- filtered_tissue_specific[filtered_tissue_specific$tissue == i, ]
  df_GTEx <- df[df$source == "GTEx", ]
  df_sra <-  df[df$source == "sra", ]
  df_tcga <-  df[df$source == "tcga", ]
  
  df_GTEx$tissue_subtype <- str_to_title(df_GTEx$tissue_subtype)
  df_GTEx$tissue_subtype <- gsub(" ", "_", df_GTEx$tissue_subtype)
  unique_GTEx <- unique(df_GTEx$tissue_subtype)
  for(j in unique_GTEx){
    if(j == "Blood"){
      j <- "Whole_Blood"
    }
    if(j == "Brain_Spinal_Cord_Cervical_C-1" ){
      j <- "Brain_Spinal_cord_cervical_c-1" 
    }
    if(j == "Brain_Frontal_Cortex_Ba9" ){
      j <- "Brain_Frontal_Cortex_BA9" 
    }
    if(j == "Brain_Anterior_Cingulate_Cortex_Ba24" ){
      j <- "Brain_Anterior_cingulate_cortex_BA24" 
    }
    if(j == "Brain_Caudate_Basal_Ganglia" ){
      j <- "Brain_Caudate_basal_ganglia"
    }
    if(j == "Cerebellum" ){
      j <- "Brain_Cerebellum"
    }
    if(j == "Brain_Nucleus_Accumbens_Basal_Ganglia" ){
      j <- "Brain_Nucleus_accumbens_basal_ganglia"
    }
    if(j == "Brain_Putamen_Basal_Ganglia" ){
      j <- "Brain_Putamen_basal_ganglia"
    }
    if(j == "Brain_Substantia_Nigra" ){
      j <- "Brain_Substantia_nigra"
    }
    if(j == "Cells_Cultured_Fibroblasts" ){
      j <- "Cells_Cultured_fibroblasts"
    }
    if(j == "Skeletal_Muscle" ){
      j <- "Muscle_Skeletal"
    }
    if(j == "Skin_Sun_Exposed_Lower_Leg" ){
      j <- "Skin_Sun_Exposed_Lower_leg"
    }
    print(j)
    jstudy <- readRDS(paste0(datDir, "all_consensus/uncorrected_expression/GTEx/", j, ".rds"))
    saveRDS(jstudy, paste0(datDir, i, "/GTEx/uncorrected_expression/", j, ".rds"))
    saveRDS(jstudy, paste0(datDir, i, "/all/uncorrected_expression/", j, ".rds"))
  }
  unique_sra <- unique(df_sra$study)
  for(j in unique_sra){
    print(j)
    jstudy <- readRDS(paste0(datDir, "all_consensus/uncorrected_expression/sra/", j, ".rds"))
    jstudy <- jstudy[ ,colnames(jstudy) %in% df_sra$sample_id]
    saveRDS(jstudy, paste0(datDir, i, "/all/uncorrected_expression/", j, ".rds"))
  }
  
  unique_tcga <- unique(df_tcga$study)
  for(j in unique_tcga){
    print(j)
    jstudy <- readRDS(paste0(datDir, "all_consensus/uncorrected_expression/tcga/", j, ".rds"))
    jstudy <- jstudy[ ,colnames(jstudy) %in% df_tcga$sample_id]
    saveRDS(jstudy, paste0(datDir, i, "/all/uncorrected_expression/", j, ".rds"))
  }
}

