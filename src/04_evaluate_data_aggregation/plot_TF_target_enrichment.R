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
rm(list = ls())
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

get_res_df <- function(study, nStudies){
  resDir <- "/data/abattle4/prashanthi/recount3_paper/results/weighted_cov_networks/"
    res_df <- lapply(nStudies, function(n){
      readRDS(paste0(resDir, study, "/net_", n, "/TF_target_enrichment.rds"))
    })
    for(i in c(1:length(res_df))){
      res_df[[i]]$net <- study
      res_df[[i]]$nStudies <- nStudies[i]
    }
    res_df <- do.call(rbind, res_df)
    res_df
}

pal.colors <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Greens"))(11)[-c(1:4)]
sra_normal_res <- get_res_df("sra_normal", c(1, 100, 200, 300, 400, 500, 566))
sra_normal_res$nStudies <- factor(sra_normal_res$nStudies, c(1, 100, 200, 300, 400, 500, 566))

ggplot(sra_normal_res[sra_normal_res$n.edges >= 500 & sra_normal_res$n.edges <= 100000, ], aes(x = recall, y = precision, color = nStudies)) + 
  geom_point()+geom_line() + scale_colour_manual(name = "Number of\nStudies\naggregated", values = pal.colors) + theme_classic() + xlab("Recall") +
  ylab("Precision") + ggtitle("SRA") + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                            legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
                            axis.title = element_text(size = 12, face = "bold"), legend.position = "right")

resDir <- "/data/abattle4/prashanthi/recount3_paper/results/weighted_cov_networks/"
all_consensus_res <- readRDS(paste0(resDir, "all_consensus/net_966/TF_target_enrichment.rds"))
normal_consensus_res <- readRDS(paste0(resDir, "normal_consensus/net_629/TF_target_enrichment.rds"))
cancer_consensus_res <- readRDS(paste0(resDir, "cancer_consensus/net_386/TF_target_enrichment.rds"))

all_consensus_res$net <- "Universal consensus"
normal_consensus_res$net <- "Non-cancer consensus"
cancer_consensus_res$net <- "Cancer consensus"
consensus_res <- rbind(all_consensus_res, normal_consensus_res, cancer_consensus_res)
consensus_res$net <- factor(consensus_res$net, levels = c("Universal consensus", "Non-cancer consensus", "Cancer consensus"))
ggplot(consensus_res[consensus_res$n.edges >= 500 & consensus_res$n.edges <= 100000, ], aes(x = recall, y = precision, color = net)) + 
  geom_point()+geom_line() + scale_colour_brewer(palette = "Dark2") + theme_classic() + xlab("Recall") +
  ylab("Precision") + ggtitle("Precision-recall curve of consensus networks") + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                            legend.text = element_text(size = 12, face = "bold"), legend.title = element_blank(), 
                                            axis.title = element_text(size = 12, face = "bold"), legend.position = "right")

blood_res <- readRDS(paste0(resDir, "blood/all/net_65/TF_target_enrichment.rds"))
cns_res <- readRDS(paste0(resDir, "central_nervous_system/all/net_53/TF_target_enrichment.rds"))
adipose_res <- readRDS(paste0(resDir, "adipose/all/net_11/TF_target_enrichment.rds"))
skin_res <- readRDS(paste0(resDir, "skin/all/net_20/TF_target_enrichment.rds"))
liver_res <- readRDS(paste0(resDir, "liver/all/net_28/TF_target_enrichment.rds"))
lung_res <- readRDS(paste0(resDir, "lung/all/net_10/TF_target_enrichment.rds"))

blood_res$Tissue <- "Blood"
cns_res$Tissue <- "Central nervous system"
adipose_res$Tissue <- "Adipose"
skin_res$Tissue <- "Skin"
liver_res$Tissue <- "Liver"
lung_res$Tissue <- "Lung"

tissue_res <- rbind(blood_res, cns_res, adipose_res, skin_res, liver_res, lung_res)

ggplot(tissue_res[consensus_res$n.edges >= 500 & consensus_res$n.edges <= 100000, ], aes(x = recall, y = precision, color = Tissue)) + 
  geom_point()+geom_line() + scale_colour_brewer(palette = "Dark2") + theme_classic() + xlab("Recall") +
  ylab("Precision") + ggtitle("Precision-recall curve of\ntissue-context specific networks") + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                                                                      legend.text = element_text(size = 12, face = "bold"), legend.title = element_blank(), 
                                                                                      axis.title = element_text(size = 12, face = "bold"), legend.position = "right")


