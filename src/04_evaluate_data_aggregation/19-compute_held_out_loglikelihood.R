# Description
# This code reads networks of varying densities from a particular aggregation setup
# and held-out test data to compute held out log-likelihood 
# To be used for evaluating the held-out log-likelihood of 
# 1. GTEx
# 2. SRA non-cancer networks
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
rm(list = ls())
library(Matrix)
library(matrixStats)
library(reshape2)
library(mvtnorm)
library(igraph)

edge.list <- function(net){
  edge.unformatted <- lapply(net, function(inet){
    if(!is.null(dim(inet))){
      ngenes <- dim(inet)[1]
      genes <- rownames(inet)
      genes.ordered <- genes[order(genes)]
      inet <- inet[genes.ordered, genes.ordered]
      inet <- as.matrix(inet)
      inet[lower.tri(inet, diag = T)] <- NA
      inet[inet == 0] <- NA
      iedge <- reshape2::melt(inet, na.rm = T)
      iedge <- paste(iedge[,1], iedge[,2], sep = "_")
    }
    else{
      iedge <- NULL
    }
    iedge
  })
  edge.unformatted
}

# This script analyses networks
lambda <- seq(0.08, 1.00, 0.02)
lambda <- sprintf("%.2f", lambda)

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
res.dir <- paste0(homeDir, "results/")
dat.dir <- paste0(homeDir, "/data")

inputArgs <- commandArgs(TRUE)

study <- inputArgs[1]
agg_type <- inputArgs[2]
nStudies <- inputArgs[3]

res.dir <- paste0(res.dir, agg_type, "/")

if(study == "GTEx"){
  test.study <- "sra_normal"
}else{
  if(study == "sra_normal"){
    test.study <- "GTEx"
  }
}

# Read in the networks
index <- 1
net <- list()
for(ilambda in lambda){
  print(ilambda)
  inet <- readRDS(paste0(res.dir, study, "/net_", nStudies, "/lambda_", ilambda, ".rds"))
  if(is.list(inet)){
    net[[index]] <- Matrix(inet[[1]], sparse = T)
  }else{
    net[[index]] <- Matrix(inet, sparse = T)
  }
  index <- index + 1
}

# Find edges
inferred.networks <- edge.list(net)

n.edges <- lapply(inferred.networks, function(iedge){
  length(iedge)
})
n.edges <- unlist(n.edges)

# Read in the test data
gtex_metaData <- readRDS(paste0(dat.dir, "GTEx/study_metaData.rds"))
sra_metaData  <- readRDS(paste0(dat.dir, "sra_normal/study_metaData.rds"))
select.GTEx <- c("Adipose_Subcutaneous", "Cells_Cultured_fibroblasts", "Heart_Left_Ventricle", "Artery_Aorta", "Colon_Sigmoid")
select.sra <- c("SRP116272", "SRP151763", "SRP150552", "SRP174638", "SRP187978")
test.covariances <- list()

if(test.study == "sra_normal"){
  test.covariances <- lapply(select.sra, function(istudy){
    readRDS(paste0(dat.dir, test.study, "/single_covariances/", istudy, ".rds"))
  })
  names(test.covariances) <- select.sra
}
if(test.study == "GTEx"){
  test.covariances <- lapply(select.GTEx, function(istudy){
    readRDS(paste0(dat.dir, test.study, "/single_covariances/", istudy, ".rds"))
  })
  names(test.covariances) <- select.GTEx
}

num.studies <- length(test.covariances)
if(test.study == "sra_normal"){
  ntest <- sra_metaData$n[match(names(test.covariances), 
                                sra_metaData$study)]
}
if(test.study == "GTEx"){
  ntest <- gtex_metaData$n[match(names(test.covariances),gtex_metaData$study)]
}

common.genes <- intersect(rownames(net[[1]]), rownames(test.covariances[[1]]))
net <- lapply(net, function(inet){
  inet <- inet[rownames(inet) %in% common.genes, ]
  inet <- inet[ ,colnames(inet) %in% common.genes]
  inet <- inet[match(common.genes, rownames(inet)), ]
  inet <- inet[ ,match(common.genes, colnames(inet))]
  inet
})

# Obtain test covariances
test.covariances <- lapply(test.covariances, function(icov){
  icov <- icov[rownames(icov) %in% common.genes, ]
  icov <- icov[ ,colnames(icov) %in% common.genes]
  icov <- icov[match(common.genes, rownames(icov)), ]
  icov <- icov[ ,match(common.genes, colnames(icov))]
  icov
})

# Compute held-out log-lieklihood for each test study
index <- 1
loglikelihood <- list()
for(itest in test.covariances){
  print(index)
  iloglikelihood <- lapply(net, function(theta){
    ll <-  ntest[index]*(Matrix::determinant(theta, logarithm=T)$modulus[1] - sum(itest * theta))
    ll
  })
  loglikelihood[[index]]<- unlist(iloglikelihood)
  index <- index + 1
}

ll_df <- cbind(loglikelihood[[1]], loglikelihood[[2]], 
               loglikelihood[[3]], loglikelihood[[4]], 
               loglikelihood[[5]])
ll_mean <- rowMeans(ll_df)

res <- data.frame(lambda, n.edges, ll_mean)

saveRDS(res, paste0(res.dir, study, "/net_", nStudies, "/ll_", study.id, ".rds"))


