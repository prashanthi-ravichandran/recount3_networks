# Description
# This script is used to sequentially aggregate and estimate corrected data 
# for GTEx and SRA non-cancerous networks
rm(list = ls())
# Author: Prashanthi Ravichandran
# Packages
library(dplyr)
library(tidyverse)
library(coop)
library(parallel)
library(rsvd)
library(stats)
library(quadprog)
library(recount3)
library(pracma)
library(sva)
library(matrixStats)

# Input
# 1. If the input flag is SRA then use data/sra_normal/uncorrected_expression
# 2. If the input flag is GTEx then use data/GTEx/uncorrected_expression

common.genes <- function(expr_list){
  genes <- rownames(expr_list[[1]])
  for(i in c(1:length(expr_list))){
    genes <- intersect(genes, rownames(expr_list[[i]]))
  }
  genes
}

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
  #print(est.pc.rm)
  usv <- svd(dat.expr)
  #loadings <- compute.pc.loadings(dat.expr.filt)
  loadings <- usv$u
  n.pc <- c(1:est.pc.rm)
  frac_variance <- sum(usv$d[n.pc])/sum(usv$d)
  print(paste("removing", est.pc.rm, "PCs", nrow(dat.expr)))
  print(paste("Variance regressed is", frac_variance, sep = " "))
  ## use residuals from top n.pc principal components
  dat.expr.adjusted <- lm(dat.expr ~ loadings[,n.pc])$residuals
  t(dat.expr.adjusted)
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

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
sraDir <- paste0(homeDir, "data/sra_normal/")
sra_uncorr <- paste0(sraDir, "uncorrected_expression/")
sra_corr <- paste0(sraDir, "correcting_after_pooling/")
GTExDir <- paste0(homeDir, "data/GTEx/")
GTEx_uncorr <- paste0(GTExDir, "uncorrected_expression/")
GTEx_corr <- paste0(GTExDir, "correcting_after_pooling/")

inputArgs <- commandArgs(TRUE)
if(inputArgs[1] == "GTEx"){
  expr_files <- list.files(GTEx_uncorr)
  expr_files <- expr_files[!expr_files == "Cells_Leukemia_cell_line_CML.rds"]
  expr.list <- lapply(expr_files, function(ifile){
    cat(ifile, "\n")
    readRDS(paste0(GTEx_uncorr, ifile))
  })
  names(expr.list) <- gsub(".rds", "", expr_files)
  study_metaData <- readRDS(paste0(GTExDir, "/study_metaData.rds"))
  study_metaData <- study_metaData[order(study_metaData$n), ]
  expr.list <- expr.list[names(expr.list) %in% study_metaData$study]
  expr.list <- expr.list[match(study_metaData$study, names(expr.list))]
  
  counts_list <- lapply(expr.list , function(iexpr){
    counts <- assays(iexpr)$RPKM
    cat(dim(counts), "\n")
    as.matrix(counts)
  })
  nsamples <- lapply(counts_list, function(icount){
    dim(icount)[2]
  })
  
  nsamples <- unlist(nsamples)
  min_sample_size <- 15
  
  counts_list <- counts_list[nsamples >= min_sample_size]
  nstudies <- length(counts_list)
  stepSize <- 10
  aggregate_levels <- seq(0, nstudies, stepSize)[-1]
  aggregate_levels  <- unique(c(1, aggregate_levels, nstudies))
  
  cov_files <- list.files(paste0(GTExDir, "/single_covariances/"))
  icov <- readRDS(paste0(GTExDir, "/single_covariances/", cov_files[1]))
  select.genes <- rownames(icov)
  
  for(i in aggregate_levels){
  print(i)
  counts_list_subset <- counts_list[1:i]
  genes_common <- common.genes(counts_list_subset)
  counts_list_subset <- lapply(counts_list_subset, function(icount){
    icount <- icount[rownames(icount) %in% genes_common, ]
    icount[match(genes_common, rownames(icount)), ]
  })
  iexpr <- do.call(cbind, counts_list_subset)
  iexpr <- remove.zero.genes(iexpr)
  iexpr <- limma::normalizeQuantiles(iexpr)
  iexpr <- remove.zero.genes(iexpr)
  iexpr <- t(scale(t(iexpr)))
  iexpr_corrected <- pc_correct_mat(iexpr)
  iexpr_corrected <- iexpr_corrected[rownames(iexpr_corrected) %in% select.genes, ]
  iexpr_corrected <- limma::normalizeQuantiles(iexpr_corrected)
  cov <- covar(scale(t(iexpr_corrected)))
  print(paste0(GTEx_corr, "cov_",as.character(i),".rds"))
  saveRDS(cov, paste0(GTEx_corr, "cov_",as.character(i),".rds"))
  rm(cov)
  }
}else{
  expr_files <- list.files(sra_uncorr)
  expr.list <- lapply(expr_files, function(ifile){
    cat(ifile, "\n")
    readRDS(paste0(sra_uncorr, ifile))
  })
  names(expr.list) <- gsub(".rds", "", expr_files)
  study_metaData <- readRDS(paste0(sraDir, "/study_metaData.rds"))
  study_metaData <- study_metaData[order(study_metaData$n), ]
  expr.list <- expr.list[names(expr.list) %in% study_metaData$study]
  expr.list <- expr.list[match(study_metaData$study, names(expr.list))]
  
  counts_list <- lapply(expr.list , function(iexpr){
    counts <- assays(iexpr)$RPKM
    cat(dim(counts), "\n")
    as.matrix(counts)
  })
  nsamples <- lapply(counts_list, function(icount){
    dim(icount)[2]
  })
  
  nsamples <- unlist(nsamples)
  min_sample_size <- 15
  
  counts_list <- counts_list[nsamples >= min_sample_size]
  nstudies <- length(counts_list)
  stepSize <- 10
  aggregate_levels <- seq(0, nstudies, stepSize)[-1]
  aggregate_levels  <- unique(c(1, aggregate_levels, nstudies))
  
  cov_files <- list.files(paste0(sraDir, "/single_covariances/"))
  icov <- readRDS(paste0(sraDir, "/single_covariances/", cov_files[1]))
  select.genes <- rownames(icov)
  
  for(i in aggregate_levels){
    print(i)
    counts_list_subset <- counts_list[1:i]
    genes_common <- common.genes(counts_list_subset)
    counts_list_subset <- lapply(counts_list_subset, function(icount){
      icount <- icount[rownames(icount) %in% genes_common, ]
      icount[match(genes_common, rownames(icount)), ]
    })
    iexpr <- do.call(cbind, counts_list_subset)
    iexpr <- remove.zero.genes(iexpr)
    iexpr <- limma::normalizeQuantiles(iexpr)
    iexpr <- remove.zero.genes(iexpr)
    iexpr <- t(scale(t(iexpr)))
    iexpr_corrected <- pc_correct_mat(iexpr)
    iexpr_corrected <- iexpr_corrected[rownames(iexpr_corrected) %in% select.genes, ]
    iexpr_corrected <- limma::normalizeQuantiles(iexpr_corrected)
    cov <- covar(scale(t(iexpr_corrected)))
    print(paste0(sra_corr, "cov_",as.character(i),".rds"))
    saveRDS(cov, paste0(sra_corr, "cov_",as.character(i),".rds"))
    rm(cov)
  }
}