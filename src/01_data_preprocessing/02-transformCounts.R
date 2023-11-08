# Description
# This script does the following tasks
# 1. Selects samples that have expression in atleast 1 gene
# 2. Transform gene expression counts to RPKM
# 3. Subsets the genes to protein coding autosomal genes that are non-overlapping
# 4. Log transforms RPKM gene expression
# 5. Saves the processed expression data

# Input: Raw downloaded data from 01-downloadData.R
rm(list = ls())
library(recount3)

select.genes <- function(rse.object){
  counts <- assays(rse.object)$RPKM
  sum_genes <- rowSums(counts)
  rse.object <- rse.object[!sum_genes == 0,]
  rse.object
}

select.samples <- function(rse.object, threshold, ...){
  counts <- assays(rse.object)$raw_counts
  keep <- apply(counts, 2, function(x){
    t = sum(x == 0) < dim(counts)[1]
    t
  })
  rse.object <- rse.object[ ,keep]
}

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
rawDir <- paste0(homeDir, "data/raw/")
rpkmDir <- paste0(homeDir, "data/rpkm/gene_filters/")
plotsDir <- paste0(homeDir, "plots/")

pc.genes <- read.delim(paste0(homeDir, "data/protein_coding.txt"), 
                       header = F, stringsAsFactors = F)
pc.genes <- pc.genes$V2[!pc.genes$V1 %in% c("chrM","chrY")]

overlapping_genes <- read.delim(paste0(homeDir, "data/ensembl_ids_overlapping_genes.txt"),
                                stringsAsFactors = F)

human_projects <- available_projects()
human_projects <- human_projects[human_projects$n_samples >= 30, ]

GTEx_info <- subset(human_projects, file_source == "gtex")
sra_info <- subset(human_projects, file_source == "sra")
tcga_info <- subset(human_projects, file_source == "tcga")

# Compute RPKM data
lapply(1:dim(GTEx_info)[1], function(i, pc_genes, ov_genes, dirName){
  print(i)
  rse <- readRDS(paste0(rawDir, "GTEx/", GTEx_info$project[i], ".rds"))
  rse <- select.samples(rse)
  assays(rse)$counts <- transform_counts(rse, by = "auc")
  assays(rse)$RPKM <- recount::getRPKM(rse, length_var = "score")
  rse <- select.genes(rse)
  rse <- rse[rownames(rse) %in% pc_genes,]
  rse <- rse[!rownames(rse) %in% ov_genes,]
  counts <- assays(rse)$RPKM
  assays(rse)$RPKM<- log2(counts+1)
  print(dim(rse))
  saveRDS(rse, file = paste(dirName, "/GTEx/", GTEx_info$project[i],".rds", sep = ""))
}, pc.genes, overlapping_genes$x, rpkmDir)

lapply(1:dim(tcga_info)[1], function(i, pc_genes, ov_genes, dirName){
  print(i)
  rse <- readRDS(paste0(rawDir, "tcga/", tcga_info$project[i], ".rds"))
  rse <- select.samples(rse)
  assays(rse)$counts <- transform_counts(rse, by = "auc")
  assays(rse)$RPKM <- recount::getRPKM(rse, length_var = "score")
  rse <- select.genes(rse)
  rse <- rse[rownames(rse) %in% pc_genes,]
  rse <- rse[!rownames(rse) %in% ov_genes,]
  counts <- assays(rse)$RPKM
  assays(rse)$RPKM<- log2(counts+1)
  print(dim(rse))
  saveRDS(rse, file = paste(dirName, "/tcga/", tcga_info$project[i],".rds", sep = ""))
}, pc.genes, overlapping_genes$x, rpkmDir)

studies_downloaded_separately <- c("SRP107565", "SRP104120", "SRP108320", "SRP108120",
       "SRP079357", "SRP108292", "SRP107198", "SRP107901",
       "SRP102239", "SRP102433", "SRP106195", "DRP003950",
       "SRP106050", "SRP107927", "SRP105227", "SRP106627",
       "SRP102186", "SRP106011", "SRP105266", "SRP103588",
       "SRP106875", "SRP106008", "SRP106630", "SRP107072",
       "SRP103772", "SRP103878", "SRP101959", "SRP104124",
       "SRP102685", "SRP105816", "SRP108251", "SRP105369",
       "SRP105769", "SRP108135", "SRP079342", "SRP102542",
       "SRP104148", "SRP102999", "SRP103200", "SRP102077",
       "SRP106077", "SRP108393", "SRP108032", "SRP102119",
       "SRP103819", "SRP105756", "SRP108121", "SRP106621",
       "SRP108321", "SRP103821", "SRP102952", "SRP107025",
       "SRP104125", "SRP106817", "SRP079684", "SRP107036",
       "SRP102104", "SRP103204", "SRP107867", "SRP102683",
       "SRP102483", "SRP102922", "SRP108128")

lapply(1:dim(sra_info)[1], function(i, pc_genes, ov_genes, dirName){
  print(i)
  rse <- readRDS(paste0(rawDir, "sra/", sra_info$project[i], ".rds"))
  rse <- select.samples(rse)
  rse <- rse[ ,!rse$recount_qc.bc_auc.all_reads_all_bases == 0]
  assays(rse)$counts <- transform_counts(rse, by = "auc")
  if(sra_info$project[i] %in% studies_downloaded_separately){
    assays(rse)$RPKM <- recount::getRPKM(rse)
  }else{
    assays(rse)$RPKM <- recount::getRPKM(rse, length_var = "score")
  }
  rse <- select.genes(rse)
  rse <- rse[rownames(rse) %in% pc_genes,]
  rse <- rse[!rownames(rse) %in% ov_genes,]
  counts <- assays(rse)$RPKM
  assays(rse)$RPKM<- log2(counts+1)
  print(dim(rse))
  saveRDS(rse, file = paste(dirName, "/sra/", sra_info$project[i],".rds", sep = ""))
}, pc.genes, overlapping_genes$x, rpkmDir)
