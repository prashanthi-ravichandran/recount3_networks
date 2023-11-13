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
  if(trait_set == "independent_sumstats"){
    res_files <- res_files[!res_files %in% c("PASS_BMI1", "PASS_Height1", "PASS_Type_2_Diabetes", 
                                             "PASS_Years_of_Education2", "PASS_Ever_Smoked")]
  }
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
                        "all_consensus" = 42, 
                        "normal_consensus" = 42, 
                        "blood_consensus" = 9, 
                        "CNS_consensus" = 8, 
                        "blood_GTEx" = 9,
                        "CNS_GTEx" = 8)
    }else{
      nTraits <- switch(study, 
                        "all_consensus" = 176,
                        "CNS_consensus" = 47)
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

## Main code ##
trait_set <- "UKBB"
all_consensus_0.14 <- make_plot_df("all_consensus", "0.14", trait_set)
all_consensus_0.16 <- make_plot_df("all_consensus", "0.16", trait_set)
all_consensus_0.18 <- make_plot_df("all_consensus", "0.18", trait_set)
all_consensus_0.20 <- make_plot_df("all_consensus", "0.20", trait_set)

normal_consensus_0.14 <- make_plot_df("normal_consensus", "0.14", trait_set)
normal_consensus_0.16 <- make_plot_df("normal_consensus", "0.16", trait_set)
normal_consensus_0.18 <- make_plot_df("normal_consensus", "0.18", trait_set)
normal_consensus_0.20 <- make_plot_df("normal_consensus", "0.20", trait_set)

normal_consensus_res <- rbind(normal_consensus_0.14[[2]], 
                              normal_consensus_0.16[[2]], 
                              normal_consensus_0.18[[2]], 
                              normal_consensus_0.20[[2]])

normal_consensus_res$Enrichment <- as.numeric(normal_consensus_res$Enrichment)
normal_consensus_res$Enrichment_std_error <- as.numeric(normal_consensus_res$Enrichment_std_error)
normal_consensus_res$Enrichment_pval <- as.numeric(normal_consensus_res$Enrichment_pval)
normal_consensus_res$Coefficient <- as.numeric(normal_consensus_res$Coefficient)
normal_consensus_res$Coefficient_pval <- as.numeric(normal_consensus_res$Coefficient_pval)
normal_consensus_res$Coefficient_std_error <- as.numeric(normal_consensus_res$Coefficient_std_error)
normal_consensus_res$Lambda <- factor(normal_consensus_res$Lambda, levels = c("0.14", "0.16", "0.18", "0.20"))
normal_consensus_res$Control[normal_consensus_res$Control == "Baseline"] <- "All genes + Baseline"
normal_consensus_res_1 <- normal_consensus_res[normal_consensus_res$Control == "All genes", ]
normal_consensus_res_2 <- normal_consensus_res[normal_consensus_res$Control == "All genes + Baseline", ]

saveRDS(normal_consensus_res, paste0(resDir, "normal_consensus_", trait_set, ".rds"))

p1 <- ggplot() +
  geom_bar(data = normal_consensus_res_1, mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = normal_consensus_res_2, mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = normal_consensus_res_1, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = normal_consensus_res_2, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.3, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") 


p2 <- ggplot() +
  geom_bar(data = normal_consensus_res_1, mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = normal_consensus_res_2, mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = normal_consensus_res_1, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = normal_consensus_res_2, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.4, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

if(trait_set == "UKBB"){
  title_normal <- ggdraw() + 
    draw_label(
      "Meta-Analysis of 219 UKBB Traits (Non-cancerous consensus)",
      fontface = 'bold',
      x = 0,
      size = 16,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )}

if(trait_set == "independent_sumstats"){
  title_normal <- ggdraw() + 
    draw_label(
      "Meta-Analysis of 42 Independent Traits (Non-cancerous consensus)",
      fontface = 'bold',
      x = 0,
      size = 16,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )}


legend_normal <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
plot_grid(title_normal, plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1)), legend_normal, nrow = 3, rel_heights =  c(0.1, 1, 0.1))


all_consensus_res <- rbind(all_consensus_0.14[[2]], 
                           all_consensus_0.16[[2]], 
                           all_consensus_0.18[[2]], 
                           all_consensus_0.20[[2]])

all_consensus_res$Enrichment <- as.numeric(all_consensus_res$Enrichment)
all_consensus_res$Enrichment_std_error <- as.numeric(all_consensus_res$Enrichment_std_error)
all_consensus_res$Enrichment_pval <- as.numeric(all_consensus_res$Enrichment_pval)
all_consensus_res$Coefficient <- as.numeric(all_consensus_res$Coefficient)
all_consensus_res$Coefficient_pval <- as.numeric(all_consensus_res$Coefficient_pval)
all_consensus_res$Coefficient_std_error <- as.numeric(all_consensus_res$Coefficient_std_error)
all_consensus_res$Lambda <- factor(all_consensus_res$Lambda, levels = c("0.14", "0.16", "0.18", "0.20"))
all_consensus_res$Control[all_consensus_res$Control == "Baseline"] <- "All genes + Baseline"
all_consensus_res_1 <- all_consensus_res[all_consensus_res$Control == "All genes", ]
all_consensus_res_2 <- all_consensus_res[all_consensus_res$Control == "All genes + Baseline", ]

p3 <- ggplot() +
  geom_bar(data = all_consensus_res_1, mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = all_consensus_res_2, mapping = aes(x = Centrality, y = Enrichment, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = all_consensus_res_1, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = all_consensus_res_2, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_std_error, ymax = Enrichment + Enrichment_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.3, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")


p4 <- ggplot() +
  geom_bar(data = all_consensus_res_1, mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = all_consensus_res_2, mapping = aes(x = Centrality, y = Coefficient, fill = Lambda, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(palette = "BrBG") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.4) +
  geom_errorbar(data = all_consensus_res_1, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = all_consensus_res_2, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_std_error, ymax = Coefficient + Coefficient_std_error, fill = Lambda), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.4, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

if(trait_set == "UKBB"){
  title_all <- ggdraw() + 
    draw_label(
      "Meta-Analysis of 219 UKBB Traits (Universal consensus)",
      fontface = 'bold',
      x = 0,
      size = 16,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )}

if(trait_set == "independent_sumstats"){
  title_all <- ggdraw() + 
    draw_label(
      "Meta-Analysis of 42 Independent Traits (Universal consensus)",
      fontface = 'bold',
      x = 0,
      size = 16,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )}


saveRDS(all_consensus_res, paste0(resDir, "all_consensus_", trait_set, ".rds"))

legend_all <- get_legend(p3)
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
plot_grid(title_all, plot_grid(p3, p4, nrow = 1, rel_widths = c(1, 1, 0.2)), legend_all,  nrow = 3, rel_heights =  c(0.1, 1, 0.1))

consensus_res <- rbind(all_consensus_res, normal_consensus_res)
consensus_res$net[consensus_res$net == "all_consensus"] <- "Universal consensus"
consensus_res$net[consensus_res$net == "normal_consensus"] <- "Non-cancerous consensus"
consensus_res$net <- factor(consensus_res$net, levels = c("Universal consensus", "Non-cancerous consensus"))
consensus_res_1 <- consensus_res[consensus_res$Control == "All genes", ]
consensus_res_2 <- consensus_res[consensus_res$Control == "All genes + Baseline", ]

consensus_res_mean_1 <- consensus_res_1 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), mean)
consensus_res_sd_1 <- consensus_res_1 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), sd)
consensus_res_mean_1$Enrichment_sd <- consensus_res_sd_1$Enrichment
consensus_res_mean_1$Coefficient_sd <- consensus_res_sd_1$Coefficient

consensus_res_mean_2 <- consensus_res_2 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), mean)
consensus_res_sd_2 <- consensus_res_2 %>% group_by(Centrality, Control, net) %>% summarise_at(vars("Enrichment", "Coefficient", ), sd)
consensus_res_mean_2$Enrichment_sd <- consensus_res_sd_2$Enrichment
consensus_res_mean_2$Coefficient_sd <- consensus_res_sd_2$Coefficient


p5 <- ggplot() +
  geom_bar(data = consensus_res_mean_1, mapping = aes(x = Centrality, y = Enrichment, fill = net, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = consensus_res_mean_2, mapping = aes(x = Centrality, y = Enrichment, fill = net, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(name = "Network", palette = "Set2") + theme_classic2() + xlab("") + ylab("Enrichment") +
  geom_hline(yintercept=1, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = consensus_res_mean_1, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_sd, ymax = Enrichment + Enrichment_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = consensus_res_mean_2, mapping = aes(x = Centrality, ymin = Enrichment - Enrichment_sd, ymax = Enrichment + Enrichment_sd, fill = net), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.3, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

p6 <- ggplot() +
  geom_bar(data = consensus_res_mean_1, mapping = aes(x = Centrality, y = Coefficient, fill = net, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  geom_bar(data = consensus_res_mean_2, mapping = aes(x = Centrality, y = Coefficient, fill = net, alpha = Control), stat = 'identity', position = 'dodge', colour = "black") +
  scale_fill_brewer(name = "Network", palette = "Set2") + theme_classic2() + xlab("") + ylab("Coefficient") +
  geom_hline(yintercept=0, colour = "black", linetype = 2, linewidth = 0.6) +
  geom_errorbar(data = consensus_res_mean_1, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_sd, ymax = Coefficient + Coefficient_sd, fill = net), width = 0.2, position=position_dodge(0.9)) +
  geom_errorbar(data = consensus_res_mean_2, mapping = aes(x = Centrality, ymin = Coefficient - Coefficient_sd, ymax = Coefficient + Coefficient_sd, fill = net), width = 0.2, position=position_dodge(0.9)) + 
  scale_alpha_manual(values = c("All genes" = 0.3, "All genes + Baseline" = 1)) + 
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom")

if(trait_set == "UKBB"){
  title_consensus <- ggdraw() + 
    draw_label(
      "Meta-Analysis of 219 UKBB Traits",
      fontface = 'bold',
      x = 0,
      size = 16,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )}

if(trait_set == "independent_sumstats"){
  title_consensus <- ggdraw() + 
    draw_label(
      "Meta-Analysis of 42 Independent Traits",
      fontface = 'bold',
      x = 0,
      size = 16,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )}


legend_consensus <- get_legend(p5)
p5 <- p5 + theme(legend.position = "none")
p6 <- p6 + theme(legend.position = "none")
plot_grid(title_consensus, plot_grid(p5, p6, nrow = 1, rel_widths = c(1, 1, 0.2)), legend_consensus,  nrow = 3, rel_heights =  c(0.1, 1, 0.1))

saveRDS(consensus_res, paste0(resDir, "consensus_", trait_set, ".rds"))





