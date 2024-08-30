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
library(ggplot2)
library(cowplot)

blood_GTEx_res <- readRDS("/data/abattle4/prashanthi/recount3_paper/results/weighted_cov_networks/blood/GTEx/net_1/oddsRatio.rds")
blood_all_res <- readRDS("/data/abattle4/prashanthi/recount3_paper/results/weighted_cov_networks/blood/all/net_65/oddsRatio.rds")

blood_GTEx_res$Study <- "GTEx"
blood_all_res$Study <- "All samples"

blood_res <- rbind(blood_GTEx_res, blood_all_res)


scientific_10 <- function(x) {   parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) }

p1 <- ggplot(blood_res, aes(x = Lambda, y = F1, color = Study)) + geom_point(aes(size = nEdges)) + geom_line() + theme_classic() +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) + xlab("Lambda") + ylab("F1-score") + ggtitle("F1 score vs.\nPenalization parameter") + 
    scale_size_binned(breaks=c(1000,10000,100000,200000,300000),  labels=function(x) scientific_10(x)) + 
  guides(size = guide_bins(title = "Number of edges", theme = theme(legend.axis.line = element_blank())), 
         color = guide_legend(title = "Input Data"))

p2 <- ggplot(blood_res, aes(x = recall, y = precision, color = Study)) + geom_point(aes(size = nEdges)) + geom_line() + theme_classic() +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) + xlab("Recall") + ylab("Precision") + ggtitle("Precision vs. Recall") + 
  scale_size_binned(breaks=c(1000,10000,100000,200000,300000),  labels=function(x) scientific_10(x)) + 
  guides(size = guide_bins(title = "Number of edges", theme = theme(legend.axis.line = element_blank())), 
         color = guide_legend(title = "Input Data")) 

max_all_F1 <- max(blood_all_res$F1)
max_GTEx_F1 <- max(blood_GTEx_res$F1)
plot_df <- data.frame("Study" = c("GTEx (optimal F1)", "All samples (optimal F1)", "GTEx (Lambda = 0.12)"), 
                      "F1" = c(max_GTEx_F1, max_all_F1, blood_GTEx_res$F1[blood_GTEx_res$Lambda == "0.12"]))
plot_df$Study <- factor(plot_df$Study, levels = c("All samples (optimal F1)", "GTEx (optimal F1)", "GTEx (Lambda = 0.12)"))
p3 <- ggplot(plot_df, aes(x = Study, y = F1, fill = Study)) + geom_bar(stat = "identity") + theme_classic()+
  xlab("Study") + ylab("F1 score") + ggtitle("F1 score") + theme(plot.title = element_text(face="bold", size = 16), axis.text.y = element_text(size = 12, face = "bold"),
                                                                 axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                                                         legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"),
                                                                         axis.title = element_text(size = 12, face = "bold"))


plot_grid(plot_grid(p2, p1, nrow = 1), p3, nrow = 2, rel_heights = c(1, 0.7))

