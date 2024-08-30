# Description
# This code reads in a user specified covariance matrix and penalization parameter lambda
# to infer a network using graphical lasso
rm(list = ls())
.libPaths(c("/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-markdown-1.1-g65guwxov2k2maedod2wedfr5jon6egu/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-knitr-1.28-cn7dhiz6mwl53op4gpecto35sljr4muz/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-yaml-2.2.0-ttbd4ipwa5jccbl733saf7c2zijmdixr/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-stringr-1.4.0-qbr2amu2xnxkicncfgvu7klzll4dg46v/rlib/R/library",
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-stringi-1.4.3-n22ruwgbot2i3becz3comesic75s47r6/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-magrittr-1.5-wy6q2ditqph62m4dib33mds3wsbluj7g/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-glue-1.4.1-5ejojwbzonzsvkmwa7o2h57gpbsutbl4/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-mime-0.7-n4hlpgd2g5fdbt3idl374eftb2utep3r/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-highr-0.8-rcw4hovy72yo4wj6mfrgtgbkovlqb3fg/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-evaluate-0.14-laasv6wsxntd3ht34ydl4ott7zzfthko/rlib/R/library", 
            "/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-4.0.2-amdvcpog4ugspqwwx3ari7pzkmckelu6/rlib/R/library", 
            "/home/pravich2/rlibs/4.0.2/gcc/9.3.0"))
# Load required libraries
#library(huge)
library("QUIC") # package for graphical lasso
library("parallel")
library(Matrix)
library(glassoFast)


## graphical lasso
# range of penalty parameter

tol = 1e-04
maxIter = 1000

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
datDir <- paste0(homeDir, "data/")
resDir <- paste0(homeDir, "results/")
## Input arguments
inputargs <- commandArgs(TRUE)
lambda <- as.numeric(inputargs[1])
print(inputargs)

dat.fn <- paste0(datDir, "blood/all/single_covariances/Whole_Blood.rds")
res.fn <- paste0(resDir, "single_study_networks/Whole_Blood/lambda_", sprintf("%.2f", lambda), ".rds")

# Read in the covariance matrix
cov.mat    <- readRDS(dat.fn)

print("Data read complete")

# learn co-expression network with graphical lasso
if(lambda == 0){
  dat.net <- solve(cov.mat, tol = 1e-100)
}else{
  eachRho.mat <- matrix(lambda, nrow = dim(cov.mat)[1], ncol = dim(cov.mat)[2])
  diag(eachRho.mat) <- 0
  dat.net <- QUIC(cov.mat, rho = eachRho.mat, tol = tol, maxIter = maxIter, msg = 4)$X
  colnames(dat.net) <- colnames(cov.mat)
  rownames(dat.net) <- rownames(cov.mat)
}


saveRDS(dat.net , res.fn)