# Description
# This script downloads all human projects deposited in Recount3 
# with 30 or more samples
rm(list = ls())
# Author: Prashanthi Ravichandran
# Packages
library(tidyverse)
library(recount3)
# Required files
# data/protein_coding.txt: List of protein coding genes
# data/ensembl_ids_overlapping_genes.txt: List of genes with overlapping sequences that are excluded prior to network inference
#===============================================================================

# Code

select.genes <- function(rse.object, threshold, ...){
  counts <- SummarizedExperiment::assay(rse.object, 1)
  min.samples <- round(0.25*dim(rse.object)[2]) # default min.samples expression in at least 1/4 of samples
  keep <- apply(counts, 1, function(x, n = min.samples){
    t = sum(x >= threshold) >= n
    t
  })
  rse.object <- rse.object[keep,]
}
# Change to the location of the cloned directory
homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
dir.create(paste0(homeDir, "data/raw/"))
dir.create(paste0(homeDir, "data/rpkm/"))

rawDir <- paste0(homeDir, "data/raw/")
rpkmDir <- paste0(homeDir, "data/rpkm/")
plotsDir <- paste0(homeDir, "plots/")

pc.genes <- read.delim(paste0(homeDir, "data/protein_coding.txt"), 
                       header = F, stringsAsFactors = F)
pc.genes <- pc.genes$V2[!pc.genes$V1 %in% c("chrM","chrY")]

overlapping_genes <- read.delim(paste0(homeDir, "data/ensembl_ids_overlapping_genes.txt"),
                                stringsAsFactors = F)

human_projects <- available_projects()

pdf(paste0(plotsDir, "Distribution_of_sample_sizes_Recount3.pdf"))
hist((human_projects$n_samples), main = "Distribution of sample sizes (Recount3: 8742 studies)", xlab = "Number samples", ylab = "density")
dev.off()


human_projects <- human_projects[human_projects$n_samples >= 30, ]

GTEx_info <- subset(human_projects, file_source == "gtex")
sra_info <- subset(human_projects, file_source == "sra")
tcga_info <- subset(human_projects, file_source == "tcga")

pdf(paste0(plotsDir, "Distribution_of_sample_sizes_GTEx.pdf"))
hist((GTEx_info$n_samples), main = "Distribution of sample sizes (GTEx: 29 studies)", xlab = "Number samples", ylab = "density")
dev.off()

pdf(paste0(plotsDir, "Distribution_of_sample_sizes_sra.pdf"))
hist(sra_info$n_samples, main = "Distribution of sample sizes (sra: 1685 studies)", xlab = "Number samples", ylab = "density")
dev.off()

pdf(paste0(plotsDir, "Distribution_of_sample_sizes_tcga.pdf"))
hist(tcga_info$n_samples, main = "Distribution of sample sizes (TCGA: 33 studies)", xlab = "Number samples", ylab = "density")
dev.off()

lapply(1:dim(GTEx_info)[1], function(i, dirName){
  rse <- tryCatch(create_rse(GTEx_info[i, ]), error=function(e) NULL)
  if(!is.null(rse)){
  print(dim(rse))
  saveRDS(rse, file = paste(dirName, "GTEx/", GTEx_info$project[i],".rds", sep = ""))}
}, rawDir)


lapply(1:dim(sra_info)[1], function(i, dirName){
  print(i)
  rse <- tryCatch(create_rse(sra_info[i, ]), error=function(e) NULL)
  if(!is.null(rse)){
  print(dim(rse))
  saveRDS(rse, file = paste(dirName, "sra/", sra_info$project[i],".rds", sep = ""))}
}, rawDir)

lapply(1:dim(tcga_info)[1], function(i, dirName){
  rse <- tryCatch(create_rse(tcga_info[i, ]), error=function(e) NULL)
  if(!is.null(rse)){
  print(dim(rse))
  saveRDS(rse, file = paste(dirName, "tcga/", tcga_info$project[i],".rds", sep = ""))}
}, rawDir)



