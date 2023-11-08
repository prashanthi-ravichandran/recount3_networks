rm(list = ls())
library(recount3)
# Script to pre-process cancer consensus study
# The inputs to this includes 
# 1. RPKM transformed TCGA counts found in data/rpkm/gene_filters
# 2. SRA data from data/replicates_merged
# 3. Manual annotations of cancer status and tissue type from data/manual_annotations.tsvs


homeDir <- "/data/abattle4/prashanthi/recount3_paper"
rpkmDir <- paste0(homeDir, "data/rpkm/gene_filters/")
dir.create(paste0(homeDir, "data/cancer_consensus/"))
dir.create(paste0(homeDir, "data/cancer_consensus/uncorrected_expression"))
uncorrectedDir <- paste0(homeDir, "data/cancer_consensus/uncorrected_expression/")

studies <- c("sra", "tcga")

# For SRA we need to apply microRNA and other filters to the data with replicates merged
replicatesMerged <- paste0(homeDir, "data/replicates_merged/gene_filters/")
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

tissue_df_annotations <- readRDS(paste0(homeDir, "data/manual_annotations.tsv"))
tissue_df_annotations <- tissue_df_annotations[tissue_df_annotations$cancer == "cancer", ]

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
sra_list <- sra_list[n.sra > 0]

for(i in c(1:length(sra_list))){
  print(i)
  saveRDS(sra_list[[i]], paste0(uncorrectedDir, studies[1], "/", names(sra_list)[i], ".rds"))
}

# Finally TCGA
tcga.projects <- list.files(paste0(rpkmDir, studies[2], "/"))
tcga_list <- lapply(tcga.projects, function(iproj){
  readRDS(paste0(rpkmDir, studies[2], "/", iproj))
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


tcga_cancer_samples <- lapply(tcga_list, function(irse){
  irse[ ,!irse$tcga.cgc_sample_sample_type == "Solid Tissue Normal"]
})

tcga_nSamples <- lapply(tcga_cancer_samples, function(irse){
  dim(irse)[2]
})

tcga_nInd <- lapply(tcga_cancer_samples, function(irse){
  length(unique(irse$tcga.xml_patient_id))
})

tcga_nSamples <- unlist(tcga_nSamples)
tcga_nInd <- unlist(tcga_nInd)

tcga_cancer_samples<- tcga_cancer_samples[tcga_nSamples > 0]
for(i in c(1:length(tcga_cancer_samples))){
  print(i)
  rse <- tcga_cancer_samples[[i]]
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
  saveRDS(rse_merged, paste0(uncorrectedDir, studies[2], "/", names(tcga_cancer_samples)[i], ".rds"))
}


