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
resDir <- paste0(homeDir, "results/weighted_cov_networks/")

normal_centrality <- readRDS(paste0(resDir, "normal_consensus/net_629/centrality_0.16.rds"))
cancer_centrality <- readRDS(paste0(resDir, "cancer_consensus/net_386/centrality_0.22.rds"))

geneData <- readRDS(paste0(homeDir, "/data/geneData.rds"))

convert_to_Symbol <- function(a){
  b <- c()
  for(i in c(1:length(a))){
    b[i] <- geneData$gene_name[geneData$gene_id == a[i]]
  }
  b
}

normal_symbol <- convert_to_Symbol(normal_centrality$node)
cancer_symbol <- convert_to_Symbol(cancer_centrality$node)

normal_centrality$symbol <- normal_symbol
cancer_centrality$symbol <- cancer_symbol

normal_hubs <- normal_centrality$symbol[normal_centrality$closeness >= quantile(normal_centrality$closeness, 0.9)]
cancer_hubs <- cancer_centrality$symbol[cancer_centrality$closeness >= quantile(cancer_centrality$closeness, 0.9)]

common_hubs <- intersect(normal_hubs, cancer_hubs)
normal_hubs <- normal_hubs[!normal_hubs %in% common_hubs]
cancer_hubs <- cancer_hubs[!cancer_hubs %in% common_hubs]

common_hubs <- intersect(normal_hubs, cancer_hubs)
normal_hubs <- normal_hubs[!normal_hubs %in% common_hubs]
cancer_hubs <- cancer_hubs[!cancer_hubs %in% common_hubs]

saveRDS(cancer_hubs, paste0(resDir, "cancer_hubs.rds"))
saveRDS(normal_hubs, paste0(resDir, "normal_hubs.rds"))
saveRDS(common_hubs, paste0(resDir, "common_hubs.rds"))



