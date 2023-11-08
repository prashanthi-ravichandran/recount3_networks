# Description
# This script does the following tasks for the SRA log2(RPKM + 1) data
# 1. Exclude samples that are obtained through size-fractionation library prep
# 2. Exclude samples which have small RNA in the experiment title
# 3. Merge the replicates
# 4. Saves the processed expression data

# Input: log2(RPKM + 1) gene expression data

rm(list = ls())
library(plyr)
library(dplyr)
library(data.table)
library(recount3)
library(parallel)

nsamples <- function(rse_list){
  n <- lapply(rse_list, function(rse){
    dim(rse)[2]
  })
  n <- unlist(n)
  n
}

print(sessionInfo())

gene_filters <- "gene_filters"
homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
processedDir <- paste0(homeDir, "data/replicates_merged/", gene_filters, "/")
rpkmDir <- paste0(homeDir, "data/rpkm/", gene_filters, "/")

human_projects <- available_projects()
human_projects <- human_projects[human_projects$n_samples >= 30, ]

sra_info <- subset(human_projects, file_source == "sra")


# Read in the data
rse_list <- lapply(sra_info$project, function(iproj){
   cat(iproj, "\n")
   readRDS(paste0(rpkmDir, "sra/", iproj, ".rds"))
})

names(rse_list) <- sra_info$project

n1 <- nsamples(rse_list)

# Remove the potential cell lines 
#rse_list <- lapply(rse_list, function(rse){
#  cat(dim(rse), "\n")
#  potential_cell_lines <- c(grep("cell line", rse$sra.sample_attributes, ignore.case = TRUE), 
#                            grep("cell lines", rse$sra.sample_attributes, ignore.case = TRUE))
#  rse <- rse[ ,!c(1:dim(rse)[2]) %in% potential_cell_lines]
#  cat(dim(rse), "\n")
#  rse
#})

#n2 <- nsamples(rse_list)
#rse_list <- rse_list[n2 > 0]


rse_list <- lapply(rse_list, function(rse){
  cat(dim(rse), "\n")
  smallRNA_samples <- c(grep("size fractionation", rse$sra.library_selection, ignore.case = TRUE), 
                        grep("small RNA", rse$sra.experiment_title, ignore.case = TRUE))
  rse <- rse[ ,!c(1:dim(rse)[2]) %in% smallRNA_samples]
  cat(dim(rse), "\n")
  rse
})

n3 <- nsamples(rse_list)
rse_list <- rse_list[n3 > 0]


# Merge replicates
for(i in c(1:length(rse_list))){
  print(i)
  rse <- rse_list[[i]]
  counts <- as.data.frame(t(rse@assays@data$RPKM))
  counts$experiment <- rse$sra.experiment_acc
  print(dim(rse))
  if(length(unique(counts$experiment)) < dim(counts)[1]){
    counts <- counts %>% group_by(experiment) %>% summarise_all(median) %>% as.data.frame
    rownames(counts) <- counts$experiment
    counts$experiment <- NULL
    counts <- t(counts)
    metaData <- as.data.frame(colData(rse))
    metaData <- metaData %>% group_by(sra.experiment_acc) %>% 
      dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE), 
                       across(where(is.character), dplyr::first),
                       across(where(is.logical), all)) %>% as.data.frame
    rse_merged <- SummarizedExperiment(assays = list(RPKM = counts), colData = metaData, rowData = rowData(rse))
  }else{
    rse_merged <- rse
  }
  print(dim(rse_merged))
  saveRDS(rse_merged, paste0(processedDir, names(rse_list)[i], ".rds"))
}


