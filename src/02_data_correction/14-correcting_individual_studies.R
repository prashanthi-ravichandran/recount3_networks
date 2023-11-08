# Description
# This script is used to correct each individual study/tissue for each of the following networks
# 1. All consensus
# 2. Normal consensus
# 3. Cancer consensus
# 4. GTEx
# 5. SRA non-cancerous
# We create both corrected expression matrices for each study and a study metadata file.
rm(list = ls())
# Author: Prashanthi Ravichandran
# Packages
library(rsvd)
library(stats)
library(quadprog)
library(recount3)
library(pracma)
library(sva)
library(parallel)
library(matrixStats)

pc.estimate = function(dat){
  mod = matrix(1, nrow = dim(dat)[1], ncol = 1)
  colnames(mod) = "Intercept"
  n.pc <- num.sv(t(dat), mod, method="be")
  n.pc
}

pc_correct_mat <- function(dat){
  dat.expr <- t(dat)
  print(dim(dat.expr))
  est.pc.rm <- pc.estimate(dat.expr)
  usv <- svd(dat.expr)
  loadings <- usv$u
  n.pc <- c(1:est.pc.rm)
  frac_variance <- sum(usv$d[n.pc])/sum(usv$d)
  print(paste("removing", est.pc.rm, "PCs", nrow(dat.expr)))
  print(paste("Variance regressed is", frac_variance, sep = " "))
  ## use residuals from top n.pc principal components
  dat.expr.adjusted <- lm(dat.expr ~ loadings[,n.pc])$residuals
  list(t(dat.expr.adjusted), est.pc.rm)
}

remove.zero.genes <- function(dat.expr){
  # Remove genes for which scaled expression is NA
  # This just means that the variance was zero and all the samples had the same value
  # In this case, since we have log transformed the raw expression data to give log2(expr + 1)
  # All the samples had 0 raw expression in these genes. We exclude these genes and proceed with the analysis
  gene_variances <- apply(dat.expr, 1, var)
  genes.zero.variance <- names(gene_variances)[gene_variances <= 1e-6]
  dat.expr[!rownames(dat.expr) %in% genes.zero.variance, ]
}

inputArgs <- commandArgs(TRUE)
num_cores <- detectCores() - 1
agg_option <- inputArgs[1]

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
datDir <- paste0(homeDir, "data/")
if(agg_option %in% c("normal_consensus", "all_consensus")){
  studies <- c("GTEx", "sra", "tcga")
  projects <- lapply(studies, function(study){
    list.files(paste0(datDir, agg_option, "/uncorrected_expression/", study, "/"))
  })
}else{
  if(agg_option == "cancer_consensus"){
    studies <- c("sra", "tcga")
    projects <- lapply(studies, function(study){
      list.files(paste0(datDir, agg_option, "/uncorrected_expression/", study, "/"))
    })
  }else{
    projects <- list.files(paste0(datDir, agg_option, "/uncorrected_expression/"))
  }
}

if(agg_option %in% c("normal_consensus", "all_consensus")){
  expr_list <- list()
  studies <- c("GTEx", "sra", "tcga")
  expr_list[[1]] <- lapply(projects[[1]], function(iproj){
    cat(iproj, sep = "\n")
    readRDS(paste0(datDir, agg_option, "/uncorrected_expression/", studies[1], "/", iproj))
  })
  
  expr_list[[2]] <- lapply(projects[[2]], function(iproj){
    cat(iproj, sep = "\n")
    readRDS(paste0(datDir, agg_option, "/uncorrected_expression/", studies[2], "/", iproj))
  })
  
  
  names(expr_list[[1]]) <- gsub(".rds", "", projects[[1]])
  names(expr_list[[2]]) <- gsub(".rds", "", projects[[2]])
  
  all_expr <- c(expr_list[[1]], expr_list[[2]])
  previous_processed <- 1
  if(previous_processed == 1){
    completed_files <- list.files(paste0(datDir, agg_option,"/corrected_expression/"))
    completed_files <- gsub(".rds", "", completed_files)
    all_expr <- all_expr[!names(all_expr) %in% completed_files]
  }
}else{
  if(agg_option == "cancer_consensus"){
    expr_list <- list()
    studies <- c("sra", "tcga")
    expr_list[[1]] <- lapply(projects[[1]], function(iproj){
      cat(iproj, sep = "\n")
      readRDS(paste0(datDir, agg_option, "/uncorrected_expression/", studies[1], "/", iproj))
    })
    
    expr_list[[2]] <- lapply(projects[[2]], function(iproj){
      cat(iproj, sep = "\n")
      readRDS(paste0(datDir, agg_option, "/uncorrected_expression/", studies[2], "/", iproj))
    })
    
    expr_list[[3]] <- lapply(projects[[3]], function(iproj){
      cat(iproj, sep = "\n")
      readRDS(paste0(datDir, agg_option, "/uncorrected_expression/", studies[3], "/", iproj))
    })
    
    names(expr_list[[1]]) <- gsub(".rds", "", projects[[1]])
    names(expr_list[[2]]) <- gsub(".rds", "", projects[[2]])
    names(expr_list[[3]]) <- gsub(".rds", "", projects[[3]])
    
    all_expr <- c(expr_list[[1]], expr_list[[2]], expr_list[[3]])
    previous_processed <- 1
    if(previous_processed == 1){
      completed_files <- list.files(paste0(datDir, agg_option,"/corrected_expression/"))
      completed_files <- gsub(".rds", "", completed_files)
      all_expr <- all_expr[!names(all_expr) %in% completed_files]
    }
    
  }else{
  all_expr <- lapply(projects, function(iproj){
    cat(iproj, sep = "\n")
    readRDS(paste0(datDir, agg_option, "/uncorrected_expression/", iproj))
  })
  names(all_expr) <- gsub(".rds", "", projects)
  previous_processed <- 1
  if(previous_processed == 1){
    completed_files <- list.files(paste0(datDir, agg_option,"/corrected_expression/"))
    completed_files <- gsub(".rds", "", completed_files)
    all_expr <- all_expr[!names(all_expr) %in% completed_files]
  }
  }
}

# Remove zero variance genes
# For each group remove zero variances genes 
counts_list <- lapply(all_expr , function(iexpr){
  counts <- assays(iexpr)$RPKM
  cat(dim(counts), "\n")
  counts <- remove.zero.genes(counts)
  cat(dim(counts), "\n")
  as.matrix(counts)
})

nsamples <- lapply(counts_list, function(icount){
  dim(icount)[2]
})

nsamples <- unlist(nsamples)
min_sample_size <- 15

counts_list <- counts_list[nsamples >= min_sample_size]

# Quantile normalize the count matrix
counts_list <- lapply(counts_list , function(icount){
  cat(dim(icount), "\n")
  limma::normalizeQuantiles(icount)
})

counts_list <- lapply(counts_list , function(icount){
  cat(dim(icount), "\n")
  icount <- remove.zero.genes(icount)
  cat(dim(icount), "\n")
  icount
})

counts_list <- lapply(counts_list, function(icount){
  t(scale(t(icount)))
})

projects <- names(counts_list)
nsamples <- lapply(counts_list, function(icount){
  dim(icount)[2]
})
nsamples <- unlist(nsamples)
nPCs <- c()
saveDir <- paste0(datDir, agg_option,"/corrected_expression/")
for(i in c(1:length(counts_list))){
  print(i)
  res <- pc_correct_mat(counts_list[[i]])
  corr_counts <- res[[1]]
  nPCs[i] <- res[[2]]
  saveRDS(corr_counts, paste0(saveDir, names(counts_list)[i], ".rds"))
}

study_metaData <- data.frame("study" = projects, 
                             "nsamples" = nsamples, 
                             "nPCs" = nPCs)
rownames(study_metaData) <- c(1:dim(study_metaData)[1])
saveRDS(study_metaData, paste0(datDir, agg_option, "/study_metaData.rds"))

