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
library(Rtsne)
library(ggplot2)
library(cowplot)
library(dplyr)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(irlba)
library(stringr)
library(SummarizedExperiment)

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
source(paste0(homeDir, "src/plot_utils.R"))

tissue <- "skin"
resDir <- "/data/abattle4/prashanthi/recount3_paper/results/weighted_cov_networks/"
all_samples <- readRDS(paste0(resDir, tissue, "/all/net_20/oddsRatio.rds"))
GTEx_samples <- readRDS(paste0(resDir, tissue, "/GTEx/net_2/oddsRatio.rds"))

all_samples$Samples_agg <- "All samples"
GTEx_samples$Samples_agg <- "GTEx"

res_odds <- rbind(all_samples, GTEx_samples)
res_odds$TPR <- res_odds$true_positives/(res_odds$true_positives + res_odds$false_negatives)

ggplot(res_odds[res_odds$nEdges > 100 & res_odds$nEdges < 100000, ],
       aes(x = nEdges, y = TPR, color = Samples_agg)) + geom_point() + geom_line() + theme_bw() +
  labs(x=expression(log[10]~"|"~italic(E)~"|")) + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                                      legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
                                                      axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") +
  facet_zoom2(xlim = c(100, 10000), ylim = c(0, 0.01), zoom.size = 1, show.area = TRUE) + ggtitle("Skin")

