# Description
# This script evaluates whether a particular network is scale or not
# For the range of lambda values over which we infer networks in a particular aggregation setting
# we determine if the degree distribution follows scale-free properties
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
library(igraph)
library(ggplot2)
library(maptools)

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
dat.dir <- paste0(homeDir, "data/")
res.dir <- paste0(homeDir, "results/weighted_cov_networks/")
plot.dir <- paste0(homeDir, "plots/")

lambda <- seq(0.10,1.00,0.02)
inputArgs <- commandArgs(TRUE)
agg_setting <- inputArgs[1]

nstudies <- switch(  
agg_setting,  
"all_consensus"= 966,  
"normal_consensus"= 629,  
"cancer_consensus"= 386,  
"sra_normal"= 566,
"GTEx"= 49, 
"adipose/all" = 11, 
"airway/all" = 8, 
"blood/all" = 65,
"B_cells/all" = 17, 
"breast/all" = 12, 
"cardiac/all" = 13, 
"central_nervous_system/all" = 53, 
"colon/all" = 10, 
"eye/all" = 2, 
"esophagus/all" = 3, 
"fibroblasts/all" = 13, 
"hescs/all" = 21, 
"intestine/all" = 9, 
"ipscs/all" = 14, 
"kidney/all" = 26, 
"liver/all" = 28, 
"lung/all" = 10, 
"multipotent_cells/all" = 6, 
"myeloid_cells/all" = 29, 
"pancreas/all" = 3, 
"pbmcs_t_cells/all" = 58, 
"prostate/all" = 6, 
"skeletal_muscle/all" = 8, 
"stomach/all" = 4, 
"skin/all" = 20, 
"vascular/all" = 4, 
"adipose/GTEx" = 2, 
"blood/GTEx" = 1,
"B_cells/GTEx" = 1, 
"breast/GTEx" = 1, 
"cardiac/GTEx" = 2, 
"central_nervous_system/GTEx" = 13, 
"colon/GTEx" = 2, 
"esophagus/GTEx" = 3, 
"fibroblasts/GTEx" = 1, 
"intestine/GTEx" = 1, 
"kidney/GTEx" = 1, 
"liver/GTEx" = 1, 
"lung/GTEx" = 1, 
"pancreas/GTEx" = 1, 
"prostate/GTEx" = 1, 
"skeletal_muscle/GTEx" = 1, 
"stomach/GTEx" = 1, 
"skin/GTEx" = 2, 
"vascular/GTEx" = 3, 
)  

net <- list()
for(i in c(1:length(lambda))){
  ifile <- paste0(res.dir, agg_type ,"/", agg_level, "/net_", nstudies, "/lambda_", 
                  sprintf("%.2f" , lambda[i]), ".rds")
  print(ifile)
  net[[i]] <- readRDS(ifile)
}

n_edges <- lapply(net, function(inet){
  inet <- as.matrix(inet)
  inet[lower.tri(inet)] <- 0
  diag(inet) <- 0
  inet[inet!=0] <- 1
  net_igraph <- graph_from_adjacency_matrix(inet, mode = "undirected", diag = F)
  gsize(net_igraph)
})

mean_degree <- lapply(net, function(inet){
  inet <- as.matrix(inet)
  inet[lower.tri(inet)] <- 0
  diag(inet) <- 0
  inet[inet != 0] <- 1
  net_igraph <- graph_from_adjacency_matrix(inet, mode = "undirected", diag = F)
  degree <- degree(net_igraph)
  mean(degree)
})

median_degree <- lapply(net, function(inet){
  inet <- as.matrix(inet)
  inet[lower.tri(inet)] <- 0
  diag(inet) <- 0
  inet[inet != 0] <- 1
  net_igraph <- graph_from_adjacency_matrix(inet, mode = "undirected", diag = F)
  degree <- degree(net_igraph)
  median(degree)
})

max_degree <- lapply(net, function(inet){
  inet <- as.matrix(inet)
  inet[lower.tri(inet)] <- 0
  diag(inet) <- 0
  inet[inet != 0] <- 1
  net_igraph <- graph_from_adjacency_matrix(inet, mode = "undirected", diag = F)
  degree <- degree(net_igraph)
  max(degree)
})

n_edges <- unlist(n_edges)
mean_degree <- unlist(mean_degree)
median_degree <- unlist(median_degree)
max_degree <- unlist(max_degree)

degree_distribution_list <- lapply(net, function(inet){
  inet <- as.matrix(inet)
  inet[lower.tri(inet)] <- 0
  diag(inet) <- 0
  inet[inet != 0] <- 1
  net_igraph <- graph_from_adjacency_matrix(inet, mode = "undirected", diag = F)
  df <- data.frame(degree_distribution(net_igraph))
  df$degree <- c(1:dim(df)[1])
  colnames(df) <- c("p_k", "k")
  df <- df[!df$p_k == 0, ]
  df$log_k <- log10(df$k)
  df$log_p_k <- log10(df$p_k)
  df
})

lm_models <- lapply(degree_distribution_list, function(df){
  lm(log_p_k ~ log_k, data = df)
})

slope <- lapply(lm_models, function(fit){
  fit$coefficients[2]
})

r_squared <- lapply(lm_models, function(fit){
  summary(fit)$r.squared
})

lower_limit <- lapply(lm_models, function(fit){
  confint(fit)[2, 1]
})

upper_limit <- lapply(lm_models, function(fit){
  confint(fit)[2, 2]
})

slope <- unlist(slope)
r_squared <- unlist(r_squared)
lower_limit <- unlist(lower_limit)
upper_limit <- unlist(upper_limit)

sf_res <- data.frame(lambda, n_edges, mean_degree, median_degree, max_degree, 
                     lower_limit, slope, upper_limit, r_squared)
print("SF res computed")
saveRDS(sf_res, paste0(res.dir, "/", agg_setting, "/net_", nstudies, "/sf_res.rds"))


