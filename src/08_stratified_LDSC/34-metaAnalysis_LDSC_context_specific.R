# Description
# This script reads in the output from sLDSC for context-specific networks
# and compute enrichment and coefficients by meta-analyzing across traits

rm(list = ls())
.libPaths(c("/data/apps/linux-centos8-cascadelake/gcc-9.3.0/r-4.0.2-amdvcpog4ugspqwwx3ari7pzkmckelu6/rlib/R/library", 
            "/home/pravich2/rlibs/4.0.2/gcc/9.3.0" ))
library(ggpubr)
library(rmeta)
library(ggplot2)
library(cowplot)
library(dplyr)
library(rstatix)

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
resDir <- "/data/abattle4/prashanthi/recount3_paper/results/s_LDSC/"

read_res <- function(study, lambda, centrality, control, trait_set){
  cat(control, "\n")
  cat(centrality, "\n")
  res_dir <- paste0(resDir, 
                    study, "_", lambda, "/", centrality, "/results/", trait_set, "/")
  res_files <- list.files(res_dir)
  res_files <- res_files[grep(control, res_files)]
  res_files <- res_files[grep(".results", res_files)]
  res_files <- gsub(paste0(control, ".results"), "", res_files)
  my_sd <- read.table(paste0(resDir, 
                             study, "_", lambda, "/", centrality, "/ldscore/", study, "_", lambda, control ,".sd"), header = FALSE)$V1
  
  #cat("Start traits")
  cat_enrich <- lapply(res_files, function(trait){
    cat(trait, "\n")
    data <- read.table(paste(res_dir , trait, control ,".results",sep=""),h=T)  
    log <- read.table(paste(res_dir , trait, control, ".log",sep=""),h=F,fill=T) 
    h2g <- as.numeric(as.character(log[which(log$V4=="h2:"),5]))
    M <- 5961159
    data$tau <- M*my_sd*data$Coefficient/h2g
    data$tau_sd <- M*my_sd*data$Coefficient_std_error/h2g
    myenrstat    <- (h2g/M)*((data$Prop._h2/data$Prop._SNPs)-(1-data$Prop._h2)/(1-data$Prop._SNPs))
    myenrstat_z  <- qnorm(data$Enrichment_p/2) #step2
    myenrstat_sd <- myenrstat/myenrstat_z #step3
    data$enrstat <- myenrstat
    data$enrstat_sd   <- myenrstat_sd
    data$trait <- trait
    data
    #data[data$Category == "L2_1", ]
  })
  cat_enrich <- do.call("rbind", cat_enrich)
  cat_enrich$Centrality <- centrality
  #cat_enrich$Trait <- res_files
  cat_enrich
}

get_summ_res <- function(study, lambda, control, trait_set){
  if(trait_set == "UKBB"){
    nTraits <- switch(study, 
                      "all_consensus" = 219, 
                      "normal_consensus" = 219, 
                      "blood_consensus" = 0, 
                      "CNS_consensus" = 43, 
                      "blood_GTEx" = 0,
                      "CNS_GTEx" = 43)
  }else{
    if(trait_set == "independent_sumstats"){
      nTraits <- switch(study, 
                        "all_consensus" = 47, 
                        "normal_consensus" = 47, 
                        "blood_consensus" = 9, 
                        "CNS_consensus" = 8, 
                        "blood_GTEx" = 9,
                        "CNS_GTEx" = 8)
    }else{
      nTraits <- switch(study, 
                        "all_consensus" = 176,
                        "CNS_consensus" = 47, 
                        "CNS_GTEx" = 47)
    }
  }
  
  degree_enrich <- read_res(study, lambda, "degree", control, trait_set)
  closeness_enrich <- read_res(study, lambda, "closeness", control, trait_set)
  betweeness_enrich <- read_res(study, lambda, "betweeness",control, trait_set)
  eigen_centrality_enrich <- read_res(study, lambda, "eigen_centrality", control, trait_set)
  max_weight_enrich <- read_res(study, lambda, "maximum_weight", control, trait_set)
  page_rank_enrich <- read_res(study, lambda, "page_rank", control, trait_set)
  strength_enrich <- read_res(study, lambda, "strength", control, trait_set)
  
  if(!length(unique(degree_enrich$trait)) == nTraits){
    stop("Check degree")
  }
  if(!length(unique(closeness_enrich$trait)) == nTraits){
    stop("Check closeness")
  }
  if(!length(unique(betweeness_enrich$trait)) == nTraits){
    stop("Check betweeness")
  }
  if(!length(unique(eigen_centrality_enrich$trait)) == nTraits){
    stop("Check eigen centrality")
  }
  if(!length(unique(max_weight_enrich$trait)) == nTraits){
    stop("Check maximum weight")
  }
  if(!length(unique(page_rank_enrich$trait)) == nTraits){
    stop("Check page rank")
  }
  if(!length(unique(strength_enrich$trait)) == nTraits){
    stop("Check strength")
  }
  
  concat_res <- rbind(degree_enrich, closeness_enrich, betweeness_enrich, 
                      max_weight_enrich, page_rank_enrich, strength_enrich)
  concat_res$Centrality[concat_res$Centrality == "betweeness"] <- "Betweeness"
  concat_res$Centrality[concat_res$Centrality == "closeness"] <- "Closeness"
  concat_res$Centrality[concat_res$Centrality == "degree"] <- "Degree"
  concat_res$Centrality[concat_res$Centrality == "maximum_weight"] <- "Maximum weight"
  concat_res$Centrality[concat_res$Centrality == "page_rank"] <- "Page rank"
  concat_res$Centrality[concat_res$Centrality == "strength"] <- "Strength"
  
  concat_res
}

metaAnalysis <- function(res, centrality, control){
  res <- res[res$Control == control, ]
  res <- res[res$Centrality == centrality, ]
  res_list <- split.data.frame(res, res$Category)
  trait <- res_list[[1]]$trait
  res_list <- lapply(res_list, function(ires){
    ires[match(trait, ires$trait), ]
  })
  
  enr <- lapply(res_list, function(ires){
    ires$Enrichment
  })
  enr <- do.call(rbind, enr)
  enr_sd <- lapply(res_list, function(ires){
    ires$Enrichment_std_error
  })
  enr_sd <- do.call(rbind, enr_sd)
  
  enrstat <- lapply(res_list, function(ires){
    ires$enrstat
  })
  enrstat <- do.call(rbind, enrstat)
  enrstat_sd <- lapply(res_list, function(ires){
    ires$enrstat_sd
  })
  enrstat_sd <- do.call(rbind, enrstat_sd)
  
  tau <- lapply(res_list, function(ires){
    ires$tau
  })
  tau <- do.call(rbind, tau)
  tau_sd <- lapply(res_list, function(ires){
    ires$tau_sd
  })
  tau_sd <- do.call(rbind, tau_sd)
  
  Prop_SNPs <- lapply(res_list, function(ires){
    ires$Prop._SNPs
  })
  Prop_SNPs <- do.call(rbind, Prop_SNPs)
  Prop_SNPs <- Prop_SNPs[,1]
  
  enr_meta = NULL
  tau_meta = NULL
  for (i in 1:nrow(enr)){
    test1 = meta.summaries(enr[i,],enr_sd[i,],method="random")
    if (Prop_SNPs[i]==1) {
      enr_meta = rbind(enr_meta,c(test1$summary,test1$se.summary,NA)) # case of the base annotation
    } else {
      test2 = meta.summaries(enrstat[i,],enrstat_sd[i,],method="random")
      enr_meta = rbind(enr_meta,c(test1$summary,test1$se.summary,2*pnorm(-abs(test2$summary/test2$se.summary))))
    }
    test = meta.summaries(tau[i,],tau_sd[i,],method="random")
    tau_meta = rbind(tau_meta,c(test$summary,test$se.summary,2*pnorm(-abs(test$summary/test$se.summary))))
  }
  
  out = cbind(names(res_list),enr_meta,tau_meta)
  colnames(out) = c("Annotation", "Enrichment","Enrichment_std_error","Enrichment_pval","Coefficient","Coefficient_std_error","Coefficient_pval")
  out
}

make_plot_df <- function(net, lambda, trait_set){
  all_genes <- get_summ_res(net, lambda, "_all_genes", trait_set)
  baseline <- get_summ_res(net, lambda, "_baseline", trait_set)
  all_genes$Control <- "All genes"
  baseline$Control <- "Baseline"
  res <- all_genes
  res <- rbind(all_genes, baseline)
  
  degree_all_meta <- data.frame(metaAnalysis(res, "Degree", "All genes"))
  betweeness_all_meta <- data.frame(metaAnalysis(res, "Betweeness", "All genes"))
  closeness_all_meta <- data.frame(metaAnalysis(res, "Closeness", "All genes"))
  maxweight_all_meta <- data.frame(metaAnalysis(res, "Maximum weight", "All genes"))
  pagerank_all_meta <- data.frame(metaAnalysis(res, "Page rank", "All genes"))
  strength_all_meta <- data.frame(metaAnalysis(res, "Strength", "All genes"))
  
  degree_all_meta <- degree_all_meta[degree_all_meta$Annotation == "L2_1", ]
  betweeness_all_meta <- betweeness_all_meta[betweeness_all_meta$Annotation == "L2_1", ]
  closeness_all_meta <- closeness_all_meta[closeness_all_meta$Annotation == "L2_1", ]
  maxweight_all_meta <- maxweight_all_meta[maxweight_all_meta$Annotation == "L2_1", ]
  pagerank_all_meta <- pagerank_all_meta[pagerank_all_meta$Annotation == "L2_1", ]
  strength_all_meta <- strength_all_meta[strength_all_meta$Annotation == "L2_1", ]
  
  degree_all_meta$Centrality <- "Degree"
  betweeness_all_meta$Centrality <- "Betweeness"
  closeness_all_meta$Centrality <- "Closeness"
  maxweight_all_meta$Centrality <- "Maximum weight"
  pagerank_all_meta$Centrality <- "Page rank"
  strength_all_meta$Centrality <- "Strength"
  
  all_meta <- rbind(degree_all_meta, betweeness_all_meta, closeness_all_meta, 
                    maxweight_all_meta, pagerank_all_meta, strength_all_meta)
  all_meta$Control <- "All genes"
  
  degree_baseline_meta <- data.frame(metaAnalysis(res, "Degree", "Baseline"))
  betweeness_baseline_meta <- data.frame(metaAnalysis(res, "Betweeness", "Baseline"))
  closeness_baseline_meta <- data.frame(metaAnalysis(res, "Closeness", "Baseline"))
  maxweight_baseline_meta <- data.frame(metaAnalysis(res, "Maximum weight", "Baseline"))
  pagerank_baseline_meta <- data.frame(metaAnalysis(res, "Page rank", "Baseline"))
  strength_baseline_meta <- data.frame(metaAnalysis(res, "Strength", "Baseline"))
  
  degree_baseline_meta <- degree_baseline_meta[degree_baseline_meta$Annotation == "L2_2", ]
  betweeness_baseline_meta <- betweeness_baseline_meta[betweeness_baseline_meta$Annotation == "L2_2", ]
  closeness_baseline_meta <- closeness_baseline_meta[closeness_baseline_meta$Annotation == "L2_2", ]
  maxweight_baseline_meta <- maxweight_baseline_meta[maxweight_baseline_meta$Annotation == "L2_2", ]
  pagerank_baseline_meta <- pagerank_baseline_meta[pagerank_baseline_meta$Annotation == "L2_2", ]
  strength_baseline_meta <- strength_baseline_meta[strength_baseline_meta$Annotation == "L2_2", ]
  
  degree_baseline_meta$Centrality <- "Degree"
  betweeness_baseline_meta$Centrality <- "Betweeness"
  closeness_baseline_meta$Centrality <- "Closeness"
  maxweight_baseline_meta$Centrality <- "Maximum weight"
  pagerank_baseline_meta$Centrality <- "Page rank"
  strength_baseline_meta$Centrality <- "Strength"
  
  baseline_meta <- rbind(degree_baseline_meta, betweeness_baseline_meta, closeness_baseline_meta, 
                         maxweight_baseline_meta, pagerank_baseline_meta, strength_baseline_meta)
  baseline_meta$Control <- "Baseline"
  
  meta <- rbind(all_meta, baseline_meta)
  meta$net <- net
  meta$Lambda <- lambda
  
  res <- res[res$Category %in% c("L2_2", "L2_1"), ]
  res$net <- net
  res$Lambda <- lambda
  list(res, meta)
}


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


metaAnalysis_of_select_traits <- function(res, traits){
  res <- res[res$trait %in% traits, ]
  degree_all_meta <- data.frame(metaAnalysis(res, "Degree", "All genes"))
  betweeness_all_meta <- data.frame(metaAnalysis(res, "Betweeness", "All genes"))
  closeness_all_meta <- data.frame(metaAnalysis(res, "Closeness", "All genes"))
  maxweight_all_meta <- data.frame(metaAnalysis(res, "Maximum weight", "All genes"))
  pagerank_all_meta <- data.frame(metaAnalysis(res, "Page rank", "All genes"))
  strength_all_meta <- data.frame(metaAnalysis(res, "Strength", "All genes"))
  
  degree_all_meta <- degree_all_meta[degree_all_meta$Annotation == "L2_1", ]
  betweeness_all_meta <- betweeness_all_meta[betweeness_all_meta$Annotation == "L2_1", ]
  closeness_all_meta <- closeness_all_meta[closeness_all_meta$Annotation == "L2_1", ]
  maxweight_all_meta <- maxweight_all_meta[maxweight_all_meta$Annotation == "L2_1", ]
  pagerank_all_meta <- pagerank_all_meta[pagerank_all_meta$Annotation == "L2_1", ]
  strength_all_meta <- strength_all_meta[strength_all_meta$Annotation == "L2_1", ]
  
  degree_all_meta$Centrality <- "Degree"
  betweeness_all_meta$Centrality <- "Betweeness"
  closeness_all_meta$Centrality <- "Closeness"
  maxweight_all_meta$Centrality <- "Maximum weight"
  pagerank_all_meta$Centrality <- "Page rank"
  strength_all_meta$Centrality <- "Strength"
  
  all_meta <- rbind(degree_all_meta, betweeness_all_meta, closeness_all_meta, 
                    maxweight_all_meta, pagerank_all_meta, strength_all_meta)
  all_meta$Control <- "All genes"
  
  degree_baseline_meta <- data.frame(metaAnalysis(res, "Degree", "Baseline"))
  betweeness_baseline_meta <- data.frame(metaAnalysis(res, "Betweeness", "Baseline"))
  closeness_baseline_meta <- data.frame(metaAnalysis(res, "Closeness", "Baseline"))
  maxweight_baseline_meta <- data.frame(metaAnalysis(res, "Maximum weight", "Baseline"))
  pagerank_baseline_meta <- data.frame(metaAnalysis(res, "Page rank", "Baseline"))
  strength_baseline_meta <- data.frame(metaAnalysis(res, "Strength", "Baseline"))
  
  degree_baseline_meta <- degree_baseline_meta[degree_baseline_meta$Annotation == "L2_2", ]
  betweeness_baseline_meta <- betweeness_baseline_meta[betweeness_baseline_meta$Annotation == "L2_2", ]
  closeness_baseline_meta <- closeness_baseline_meta[closeness_baseline_meta$Annotation == "L2_2", ]
  maxweight_baseline_meta <- maxweight_baseline_meta[maxweight_baseline_meta$Annotation == "L2_2", ]
  pagerank_baseline_meta <- pagerank_baseline_meta[pagerank_baseline_meta$Annotation == "L2_2", ]
  strength_baseline_meta <- strength_baseline_meta[strength_baseline_meta$Annotation == "L2_2", ]
  
  degree_baseline_meta$Centrality <- "Degree"
  betweeness_baseline_meta$Centrality <- "Betweeness"
  closeness_baseline_meta$Centrality <- "Closeness"
  maxweight_baseline_meta$Centrality <- "Maximum weight"
  pagerank_baseline_meta$Centrality <- "Page rank"
  strength_baseline_meta$Centrality <- "Strength"
  
  baseline_meta <- rbind(degree_baseline_meta, betweeness_baseline_meta, closeness_baseline_meta, 
                         maxweight_baseline_meta, pagerank_baseline_meta, strength_baseline_meta)
  baseline_meta$Control <- "Baseline"
  
  meta <- rbind(all_meta, baseline_meta)
  meta$net <- res$net[1]
  meta$Lambda <- res$Lambda[1]
  meta
}

## Main code ##
all_consensus_0.14 <- make_plot_df("all_consensus", "0.14", "all_sumstats")
all_consensus_0.16 <- make_plot_df("all_consensus", "0.16", "all_sumstats")
all_consensus_0.18 <- make_plot_df("all_consensus", "0.18", "all_sumstats")
all_consensus_0.20 <- make_plot_df("all_consensus", "0.20", "all_sumstats")

blood_consensus_0.18 <- make_plot_df("blood_consensus", "0.18", "independent_sumstats")
blood_consensus_0.20 <- make_plot_df("blood_consensus", "0.20", "independent_sumstats")
blood_consensus_0.22 <- make_plot_df("blood_consensus", "0.22", "independent_sumstats")
blood_consensus_0.24 <- make_plot_df("blood_consensus", "0.24", "independent_sumstats")
blood_consensus_0.26 <- make_plot_df("blood_consensus", "0.26", "independent_sumstats")

blood_GTEx_0.24 <- make_plot_df("blood_GTEx", "0.24", "independent_sumstats")
blood_GTEx_0.26 <- make_plot_df("blood_GTEx", "0.26", "independent_sumstats")
blood_GTEx_0.28 <- make_plot_df("blood_GTEx", "0.28", "independent_sumstats")
blood_GTEx_0.30 <- make_plot_df("blood_GTEx", "0.30", "independent_sumstats")
blood_GTEx_0.32 <- make_plot_df("blood_GTEx", "0.32", "independent_sumstats")

cns_consensus_0.20 <- make_plot_df("CNS_consensus", "0.20", "all_sumstats")
cns_consensus_0.22 <- make_plot_df("CNS_consensus", "0.22", "all_sumstats")
cns_consensus_0.24 <- make_plot_df("CNS_consensus", "0.24", "all_sumstats")
cns_consensus_0.26 <- make_plot_df("CNS_consensus", "0.26", "all_sumstats")
cns_consensus_0.28 <- make_plot_df("CNS_consensus", "0.28", "all_sumstats")
cns_consensus_0.30 <- make_plot_df("CNS_consensus", "0.30", "all_sumstats")
cns_consensus_0.32 <- make_plot_df("CNS_consensus", "0.32", "all_sumstats")

cns_GTEx_0.22 <- make_plot_df("CNS_GTEx", "0.22", "all_sumstats")
cns_GTEx_0.24 <- make_plot_df("CNS_GTEx", "0.24", "all_sumstats")
cns_GTEx_0.26 <- make_plot_df("CNS_GTEx", "0.26", "all_sumstats")
cns_GTEx_0.28 <- make_plot_df("CNS_GTEx", "0.28", "all_sumstats")
cns_GTEx_0.30 <- make_plot_df("CNS_GTEx", "0.30", "all_sumstats")
cns_GTEx_0.32 <- make_plot_df("CNS_GTEx", "0.32", "all_sumstats")
cns_GTEx_0.34 <- make_plot_df("CNS_GTEx", "0.34", "all_sumstats")
cns_GTEx_0.36 <- make_plot_df("CNS_GTEx", "0.36", "all_sumstats")


blood_consensus_res <- rbind(blood_consensus_0.18[[2]], 
                             blood_consensus_0.20[[2]], 
                             blood_consensus_0.22[[2]], 
                             blood_consensus_0.24[[2]], 
                             blood_consensus_0.26[[2]])

blood_GTEx_res <- rbind(blood_GTEx_0.24[[2]], 
                        blood_GTEx_0.26[[2]], 
                        blood_GTEx_0.28[[2]], 
                        blood_GTEx_0.30[[2]], 
                        blood_GTEx_0.32[[2]])

blood_traits <- unique(blood_consensus_0.18[[1]]$trait)
all_consensus_blood_res <- rbind(metaAnalysis_of_select_traits(all_consensus_0.14[[1]], blood_traits), 
                                 metaAnalysis_of_select_traits(all_consensus_0.16[[1]], blood_traits), 
                                 metaAnalysis_of_select_traits(all_consensus_0.18[[1]], blood_traits), 
                                 metaAnalysis_of_select_traits(all_consensus_0.20[[1]], blood_traits))

blood_res <- rbind(all_consensus_blood_res, blood_consensus_res, blood_GTEx_res)
blood_res$Enrichment <- as.numeric(blood_res$Enrichment)
blood_res$Enrichment_std_error <- as.numeric(blood_res$Enrichment_std_error)
blood_res$Enrichment_pval <- as.numeric(blood_res$Enrichment_pval)
blood_res$Coefficient <- as.numeric(blood_res$Coefficient)
blood_res$Coefficient_std_error <- as.numeric(blood_res$Coefficient_std_error)
blood_res$Coefficient_pval<- as.numeric(blood_res$Coefficient_pval)
blood_res$net[blood_res$net == "all_consensus"] <- "Universal consensus"
blood_res$net[blood_res$net == "blood_consensus"] <- "Blood consensus"
blood_res$net[blood_res$net == "blood_GTEx"] <- "Blood (GTEx)"
blood_res$Control[blood_res$Control == "Baseline"] <- "All genes + Baseline"
blood_res_1 <- blood_res[blood_res$Control == "All genes", ]
blood_res_2 <- blood_res[blood_res$Control == "All genes + Baseline", ]


p1 <- ggplot() +
  geom_bar(data = blood_res_1[blood_res_1$net == "Blood consensus", ], mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = blood_res_2[blood_res_2$net == "Blood consensus", ], mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = blood_res_1[blood_res_1$net == "Blood consensus", ], mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = blood_res_2[blood_res_2$net == "Blood consensus", ], mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.3, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") 


p2 <- ggplot() +
  geom_bar(data = blood_res_1[blood_res_1$net == "Blood consensus", ], mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = blood_res_2[blood_res_2$net == "Blood consensus", ], mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = blood_res_1[blood_res_1$net == "Blood consensus", ], mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = blood_res_2[blood_res_2$net == "Blood consensus", ], mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.4, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")


title_blood_consensus <- ggdraw() + 
  draw_label(
    "Meta-Analysis of 9 Blood Traits (Blood consensus)",
    fontface = 'bold',
    x = 0,
    size = 16,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

legend_blood_consensus <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
plot_grid(title_blood_consensus, plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1)), legend_blood_consensus, nrow = 3, rel_heights =  c(0.1, 1, 0.1))

p3 <- ggplot() +
  geom_bar(data = blood_res_1[blood_res_1$net == "Blood (GTEx)", ], mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = blood_res_2[blood_res_2$net == "Blood (GTEx)", ], mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = blood_res_1[blood_res_1$net == "Blood (GTEx)", ], mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = blood_res_2[blood_res_2$net == "Blood (GTEx)", ], mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.3, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") 


p4 <- ggplot() +
  geom_bar(data = blood_res_1[blood_res_1$net == "Blood (GTEx)", ], mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = blood_res_2[blood_res_2$net == "Blood (GTEx)", ], mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = blood_res_1[blood_res_1$net == "Blood (GTEx)", ], mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = blood_res_2[blood_res_2$net == "Blood (GTEx)", ], mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.4, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

title_blood_GTEx <- ggdraw() + 
  draw_label(
    "Meta-Analysis of 9 Blood Traits (Blood GTEx)",
    fontface = 'bold',
    x = 0,
    size = 16,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

legend_blood_GTEx <- get_legend(p3)
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
plot_grid(title_blood_GTEx, plot_grid(p3, p4, nrow = 1, rel_widths = c(1, 1)), legend_blood_GTEx, nrow = 3, rel_heights =  c(0.1, 1, 0.1))


blood_res_mean_1 <- blood_res_1 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), mean)
blood_res_sd_1 <- blood_res_1 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), sd)
blood_res_mean_1$Enrichment_sd <- blood_res_sd_1$Enrichment
blood_res_mean_1$Coefficient_sd <- blood_res_sd_1$Coefficient

blood_res_mean_2 <- blood_res_2 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), mean)
blood_res_sd_2 <- blood_res_2 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), sd)
blood_res_mean_2$Enrichment_sd <- blood_res_sd_2$Enrichment
blood_res_mean_2$Coefficient_sd <- blood_res_sd_2$Coefficient

p5a <- ggplot() +
  geom_bar(data = blood_res_mean_1, mapping = aes(x = Centrality, y = Enrichment, fill = net), stat = 'identity', position = 'dodge', colour = "black", alpha = 0.3) +
  scale_fill_brewer(name = "Network", palette = "Set1") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = blood_res_mean_1, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_sd, ymax = Enrichment + Enrichment_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

p5b <- ggplot() +
  geom_bar(data = blood_res_mean_2, mapping = aes(x = Centrality, y = Enrichment, fill = net), stat = 'identity', position = 'dodge', colour = "black", alpha = 1) +
  scale_fill_brewer(name = "Network", palette = "Set1") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = blood_res_mean_2, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_sd, ymax = Enrichment + Enrichment_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

p6a <- ggplot() +
  geom_bar(data = blood_res_mean_1, mapping = aes(x = Centrality, y = Coefficient, fill = net), alpha = 0.3, stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(name = "Network", palette = "Set1") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = blood_res_mean_1, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_sd, ymax = Coefficient + Coefficient_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

p6b <- ggplot() +
  geom_bar(data = blood_res_mean_2, mapping = aes(x = Centrality, y = Coefficient, fill = net), alpha = 1, stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(name = "Network", palette = "Set1") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = blood_res_mean_2, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_sd, ymax = Coefficient + Coefficient_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")


title_blood_a <- ggdraw() + 
  draw_label(
    "Meta-Analysis of 9 Blood Related Traits (Conditioned on All genes)",
    fontface = 'bold',
    x = 0,
    size = 16,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

title_blood_b <- ggdraw() + 
  draw_label(
    "Meta-Analysis of 9 Blood Related Traits (Conditioned on All genes and Baseline)",
    fontface = 'bold',
    x = 0,
    size = 16,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

legend_blood_a <- get_legend(p5a)
p5a <- p5a + theme(legend.position = "none")
p6a <- p6a + theme(legend.position = "none")
plot_grid(title_blood_a, plot_grid(p5a, p6a, nrow = 1, rel_widths = c(1, 1)), legend_blood_a, nrow = 3, rel_heights =  c(0.1, 1, 0.1))

legend_blood_b <- get_legend(p5b)
p5b <- p5b + theme(legend.position = "none")
p6b <- p6b + theme(legend.position = "none")
plot_grid(title_blood_b, plot_grid(p5b, p6b, nrow = 1, rel_widths = c(1, 1)), legend_blood_b, nrow = 3, rel_heights =  c(0.1, 1, 0.1))


cns_indep_traits <- c("PASS_Alzheimers_Jansen2019","PASS_Epilepsy_Anney_2014","PASS_Parkinsons_23andMe_Corces2020",
                      "PASS_BIP_Stahl2019", "PASS_SmokingCessation_Liu2019","UKB_460K.body_WHRadjBMIz",
                      "PASS_BDSCZ_Ruderfer2018", "PASS_MDD_Wray2018", "PASS_DrinksPerWeek_Liu2019") #, 
#"PASS_ENIGMA2_MeanPutamen", "PASS_ReactionTime_Davies2018","PASS_Autism", "PASS_ADHD_Demontis2018")

cns_consensus_res <- rbind(metaAnalysis_of_select_traits(cns_consensus_0.20[[1]], cns_indep_traits),
                           metaAnalysis_of_select_traits(cns_consensus_0.22[[1]], cns_indep_traits),
                           metaAnalysis_of_select_traits(cns_consensus_0.24[[1]], cns_indep_traits),
                           metaAnalysis_of_select_traits(cns_consensus_0.26[[1]], cns_indep_traits), 
                           metaAnalysis_of_select_traits(cns_consensus_0.28[[1]], cns_indep_traits), 
                           metaAnalysis_of_select_traits(cns_consensus_0.30[[1]], cns_indep_traits), 
                           metaAnalysis_of_select_traits(cns_consensus_0.32[[1]], cns_indep_traits))
cns_consensus_res$Enrichment <- as.numeric(cns_consensus_res$Enrichment)
cns_consensus_res$Enrichment_std_error <- as.numeric(cns_consensus_res$Enrichment_std_error)
cns_consensus_res$Coefficient <- as.numeric(cns_consensus_res$Coefficient)
cns_consensus_res$Coefficient_std_error <- as.numeric(cns_consensus_res$Coefficient_std_error)

cns_GTEx_res <- rbind(metaAnalysis_of_select_traits(cns_GTEx_0.22[[1]], cns_indep_traits), 
                      metaAnalysis_of_select_traits(cns_GTEx_0.24[[1]], cns_indep_traits), 
                      metaAnalysis_of_select_traits(cns_GTEx_0.26[[1]], cns_indep_traits), 
                      metaAnalysis_of_select_traits(cns_GTEx_0.28[[1]], cns_indep_traits),
                      metaAnalysis_of_select_traits(cns_GTEx_0.30[[1]], cns_indep_traits), 
                      metaAnalysis_of_select_traits(cns_GTEx_0.32[[1]], cns_indep_traits), 
                      metaAnalysis_of_select_traits(cns_GTEx_0.34[[1]], cns_indep_traits)) #, 
#metaAnalysis_of_select_traits(cns_GTEx_0.36[[1]], cns_indep_traits))

all_consensus_cns_res <- rbind(metaAnalysis_of_select_traits(all_consensus_0.14[[1]], cns_indep_traits), 
                               metaAnalysis_of_select_traits(all_consensus_0.16[[1]], cns_indep_traits), 
                               metaAnalysis_of_select_traits(all_consensus_0.18[[1]], cns_indep_traits), 
                               metaAnalysis_of_select_traits(all_consensus_0.20[[1]], cns_indep_traits))

cns_res <- rbind(all_consensus_cns_res, cns_consensus_res ,cns_GTEx_res)
cns_res$Enrichment <- as.numeric(cns_res$Enrichment)
cns_res$Enrichment_std_error <- as.numeric(cns_res$Enrichment_std_error)
cns_res$Enrichment_pval <- as.numeric(cns_res$Enrichment_pval)
cns_res$Coefficient <- as.numeric(cns_res$Coefficient)
cns_res$Coefficient_std_error <- as.numeric(cns_res$Coefficient_std_error)
cns_res$Coefficient_pval<- as.numeric(cns_res$Coefficient_pval)
cns_res$net[cns_res$net == "all_consensus"] <- "Universal consensus"
cns_res$net[cns_res$net == "CNS_consensus"] <- "CNS consensus"
cns_res$net[cns_res$net == "CNS_GTEx"] <- "CNS (GTEx)"
cns_res$Control[cns_res$Control == "Baseline"] <- "All genes + Baseline"
cns_res_1 <- cns_res[cns_res$Control == "All genes", ]
cns_res_2 <- cns_res[cns_res$Control == "All genes + Baseline", ]

cns_res_mean_1 <- cns_res_1 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), mean)
cns_res_sd_1 <- cns_res_1 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), sd)
cns_res_mean_1$Enrichment_sd <- cns_res_sd_1$Enrichment
cns_res_mean_1$Coefficient_sd <- cns_res_sd_1$Coefficient

cns_res_mean_2 <- cns_res_2 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), mean)
cns_res_sd_2 <- cns_res_2 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), sd)
cns_res_mean_2$Enrichment_sd <- cns_res_sd_2$Enrichment
cns_res_mean_2$Coefficient_sd <- cns_res_sd_2$Coefficient


p7 <- ggplot() +
  geom_bar(data = cns_res_1[cns_res_1$net == "CNS (GTEx)", ], mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = cns_res_2[cns_res_2$net == "CNS (GTEx)", ], mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = cns_res_1[cns_res_1$net == "CNS (GTEx)", ], mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = cns_res_2[cns_res_2$net == "CNS (GTEx)", ], mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.3, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") 


p8 <- ggplot() +
  geom_bar(data = cns_res_1[cns_res_1$net == "CNS (GTEx)", ], mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = cns_res_2[cns_res_2$net == "CNS (GTEx)", ], mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = cns_res_1[cns_res_1$net == "CNS (GTEx)", ], mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = cns_res_2[cns_res_2$net == "CNS (GTEx)", ], mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.4, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

title_cns_GTEx <- ggdraw() + 
  draw_label(
    paste0("Meta-Analysis of ", length(cns_indep_traits),  " CNS Traits (CNS GTEx)"),
    fontface = 'bold',
    x = 0,
    size = 16,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

legend_CNS_GTEx <- get_legend(p7)
p7 <- p7 + theme(legend.position = "none")
p8 <- p8 + theme(legend.position = "none")
plot_grid(title_cns_GTEx, plot_grid(p7, p8, nrow = 1, rel_widths = c(1, 1)), legend_CNS_GTEx, nrow = 3, rel_heights =  c(0.1, 1, 0.2))

p9 <- ggplot() +
  geom_bar(data = cns_res_1[cns_res_1$net == "CNS consensus", ], mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = cns_res_2[cns_res_2$net == "CNS consensus", ], mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = cns_res_1[cns_res_1$net == "CNS consensus", ], mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = cns_res_2[cns_res_2$net == "CNS consensus", ], mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.3, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") 


p10 <- ggplot() +
  geom_bar(data = cns_res_1[cns_res_1$net == "CNS consensus", ], mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = cns_res_2[cns_res_2$net == "CNS consensus", ], mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = cns_res_1[cns_res_1$net == "CNS consensus", ], mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = cns_res_2[cns_res_2$net == "CNS consensus", ], mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.4, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

title_cns_consensus <- ggdraw() + 
  draw_label(
    paste0("Meta-Analysis of ", length(cns_indep_traits),  " CNS Traits (CNS consensus)"),
    fontface = 'bold',
    x = 0,
    size = 16,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

legend_CNS_consensus <- get_legend(p9)
p9 <- p9 + theme(legend.position = "none")
p10 <- p10 + theme(legend.position = "none")
plot_grid(title_cns_consensus , plot_grid(p9, p10, nrow = 1, rel_widths = c(1, 1)), legend_CNS_consensus, nrow = 3, rel_heights =  c(0.1, 1, 0.2))

p11a <- ggplot() +
  geom_bar(data = cns_res_mean_1, mapping = aes(x = Centrality, y = Enrichment, fill = net), stat = 'identity', position = 'dodge', colour = "black", alpha = 0.3) +
  scale_fill_brewer(name = "Network", palette = "Set1") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = cns_res_mean_1, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_sd, ymax = Enrichment + Enrichment_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

p11b <- ggplot() +
  geom_bar(data = cns_res_mean_2, mapping = aes(x = Centrality, y = Enrichment, fill = net), stat = 'identity', position = 'dodge', colour = "black", alpha = 1) +
  scale_fill_brewer(name = "Network", palette = "Set1") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = cns_res_mean_2, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_sd, ymax = Enrichment + Enrichment_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

p12a <- ggplot() +
  geom_bar(data = cns_res_mean_1, mapping = aes(x = Centrality, y = Coefficient, fill = net), alpha = 0.3, stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(name = "Network", palette = "Set1") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = cns_res_mean_1, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_sd, ymax = Coefficient + Coefficient_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

p12b <- ggplot() +
  geom_bar(data = cns_res_mean_2, mapping = aes(x = Centrality, y = Coefficient, fill = net), alpha = 1, stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(name = "Network", palette = "Set1") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = cns_res_mean_2, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_sd, ymax = Coefficient + Coefficient_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")


title_cns_a <- ggdraw() + 
  draw_label(
    "Meta-Analysis of 9 CNS Related Traits (Conditioned on All genes)",
    fontface = 'bold',
    x = 0,
    size = 16,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

title_cns_b <- ggdraw() + 
  draw_label(
    "Meta-Analysis of 9 CNS Related Traits (Conditioned on All genes and Baseline)",
    fontface = 'bold',
    x = 0,
    size = 16,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

legend_cns_a <- get_legend(p11a)
p11a <- p11a + theme(legend.position = "none")
p12a <- p12a + theme(legend.position = "none")
plot_grid(title_cns_a, plot_grid(p11a, p12a, nrow = 1, rel_widths = c(1, 1)), legend_cns_a, nrow = 3, rel_heights =  c(0.1, 1, 0.1))

legend_cns_b <- get_legend(p11b)
p11b <- p11b + theme(legend.position = "none")
p12b <- p12b + theme(legend.position = "none")
plot_grid(title_cns_b, plot_grid(p11b, p12b, nrow = 1, rel_widths = c(1, 1)), legend_cns_b, nrow = 3, rel_heights =  c(0.1, 1, 0.1))

saveRDS(blood_res, paste0(resDir, "blood_res.rds"))
saveRDS(cns_res, paste0(resDir, "cns_res.rds"))

