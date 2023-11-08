# Description
# This script reads in individual covariance matrices and computes the 
# unweighted and weighted aggregate of the covariance matrices 
rm(list = ls())
homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
datDir <- "/data/abattle4/prashanthi/recount3_paper/data/"
inputArgs <- commandArgs(TRUE)
agg_level <- inputArgs[1]
singleCovDir <- paste0(datDir, agg_level, "/single_covariances/")
unweightedCovDir <- paste0(datDir, agg_level, "/unweighted_covariances/")
weightedCovDir <- paste0(datDir, agg_level, "/weighted_covariances/")

nSum <- function(df){
  sum_n <- c()
  for(i in c(1:dim(df)[1])){
    sum_n[i] <- sum(df$n[1:i])
  }
  sum_n
}

unweighted_aggregate_step <- function(cov_list, df, stepSize){
  nstudies <- dim(df)[1]
  aggregate_levels <- seq(0, nstudies, stepSize)[-1]
  aggregate_levels  <- unique(c(1, aggregate_levels, nstudies))
  for(i in aggregate_levels){
    w <- rep(1/i, i)
    print(i)
    cov_average <- matrix(0, dim(cov_list[[1]])[1], dim(cov_list[[1]])[2])
    for(j in c(1:length(w))){
      cov_average <- cov_average + cov_list[[j]]*w[j]
    }
    dim(cov_average)
    saveRDS(cov_average, paste0(unweightedCovDir , "cov_",as.character(i),".rds"))
    rm(cov_average)}
}

weighted_aggregate_step <- function(cov_list, df, stepSize){
  nstudies <- dim(df)[1]
  aggregate_levels <- seq(0, nstudies, stepSize)[-1]
  aggregate_levels  <- unique(c(1, aggregate_levels, nstudies))
  for(i in aggregate_levels){
    w <- df$n[1:i]/sum(df$n[1:i])
    print(i)
    cov_average <- matrix(0, dim(cov_list[[1]])[1], dim(cov_list[[1]])[2])
    for(j in c(1:length(w))){
      cov_average <- cov_average + cov_list[[j]]*w[j]
    }
    dim(cov_average)
    saveRDS(cov_average, paste0(weightedCovDir, "cov_",as.character(i),".rds"))
    rm(cov_average)}
}

weighted_aggregate <- function(cov_list, df){
  nstudies <- dim(df)[1]
  w <- df$n/sum(df$n)
  cov_average <- matrix(0, dim(cov_list[[1]])[1], dim(cov_list[[1]])[2])
  for(j in c(1:length(w))){
    cov_average <- cov_average + cov_list[[j]]*w[j]
  }
  dim(cov_average)
  saveRDS(cov_average, 
          paste0(weightedCovDir ,"cov_",nstudies,".rds"))
}

unweighted_aggregate <- function(cov_list, df){
  nstudies <- dim(df)[1]
  cov_average <- matrix(0, dim(cov_list[[1]])[1], dim(cov_list[[1]])[2])
  for(j in c(1:nstudies)){
    cov_average <- cov_average + cov_list[[j]]*(1/nstudies)
  }
  dim(cov_average)
  saveRDS(cov_average, 
          paste0(unweightedCovDir ,"cov_",nstudies,".rds"))
}

study_metaData <- readRDS(paste0(datDir, agg_level, "/study_metaData.rds"))
cov_files <- list.files(paste0(singleCovDir))

cov.list <- lapply(cov_files, function(ifile){
  cat(ifile, "\n")
  readRDS(paste0(singleCovDir, ifile))
})

names(cov.list) <- gsub(".rds", "", cov_files)
study_metaData <- study_metaData[study_metaData$study %in% names(cov.list), ]
study_metaData <- study_metaData[order(study_metaData$n), ]
cov.list <- cov.list[match(study_metaData$study, names(cov.list))]
study_metaData$nSum <- nSum(study_metaData)

if(agg_level == "GTEx"){
  stepSize <- 10
  unweighted_aggregate_step(cov.list, study_metaData, stepSize)
  weighted_aggregate_step(cov.list, study_metadata, stepSize)
}else{
  if(agg_level == "sra_normal"){
    stepSize <- 100
    unweighted_aggregate_step(cov.list, study_metaData, stepSize)
    weighted_aggregate_step(cov.list, study_metadata, stepSize)
  }else{
    weighted_aggregate(cov.list, study_metaData)
  }
}

