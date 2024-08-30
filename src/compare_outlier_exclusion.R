rm(list = ls())
library(coop)
library(parallel)
library(stringr)
library(stringi)

all_samples_outlier_detection <- readRDS("/data/abattle4/prashanthi/recount3_paper/data/all_samples_outlier_detection.rds")
all_samples_outlier_detection$tissue[all_samples_outlier_detection$tissue == "Prostrate"] <- "Prostate"
unique_tissues <- unique(all_samples_outlier_detection$tissue)

all_samples_outlier_detection$source <- "TCGA"
all_samples_outlier_detection$source[all_samples_outlier_detection$study == "GTEx"] <- "GTEx"
all_samples_outlier_detection$source[grep("DRP", all_samples_outlier_detection$study)] <- "SRA"
all_samples_outlier_detection$source[grep("ERP", all_samples_outlier_detection$study)] <- "SRA"
all_samples_outlier_detection$source[grep("SRP", all_samples_outlier_detection$study)] <- "SRA"

all_samples_outlier_detection$z <- all_samples_outlier_detection$D1_z & all_samples_outlier_detection$D2_z
all_samples_outlier_detection$mad <- all_samples_outlier_detection$D1_mad & all_samples_outlier_detection$D2_mad
all_samples_outlier_detection$tukey <- all_samples_outlier_detection$D1_tukey & all_samples_outlier_detection$D2_tukey

all_samples_outlier_detection$score <- rowSums(all_samples_outlier_detection[ ,c("z", "mad", "tukey", "maha_z", "maha_mad", "maha_tukey")])

tissue_dat <- lapply(unique_tissues, function(tissue){
  tissue_outlier_detection <- all_samples_outlier_detection[all_samples_outlier_detection$tissue == tissue, ]
  df_sra <- tissue_outlier_detection[tissue_outlier_detection$source == "SRA", ]
  df_GTEx <- tissue_outlier_detection[tissue_outlier_detection$source == "GTEx", ]
  df_TCGA <- tissue_outlier_detection[tissue_outlier_detection$source == "TCGA", ]
  unique_sra <- unique(df_sra$study)
  data <- list()
  ind <- 1
  if(length(unique_sra) > 0){
    for(j in unique_sra){
      print(j)
      jstudy <- readRDS(paste0("/data/abattle4/prashanthi/recount3_paper/data/all_consensus/corrected_expression/", j, ".rds"))
      jstudy <- jstudy[ ,colnames(jstudy) %in% df_sra$sample_id]
      data[[ind]] <- as.matrix(jstudy)
      ind <- ind + 1
    }
  }
  unique_GTEx <- unique(df_GTEx$tissue_subtype)
  unique_GTEx <- str_to_title(unique_GTEx)
  unique_GTEx <- gsub(" ", "_", unique_GTEx)
  unique_GTEx[unique_GTEx == "Blood"] <- "Whole_Blood"
  unique_GTEx[unique_GTEx == "Brain_Spinal_Cord_Cervical_C-1"] <- "Brain_Spinal_cord_cervical_c-1" 
  unique_GTEx[unique_GTEx == "Brain_Frontal_Cortex_Ba9"] <- "Brain_Frontal_Cortex_BA9" 
  unique_GTEx[unique_GTEx == "Brain_Anterior_Cingulate_Cortex_Ba24"] <- "Brain_Anterior_cingulate_cortex_BA24"
  unique_GTEx[unique_GTEx == "Brain_Caudate_Basal_Ganglia"] <- "Brain_Caudate_basal_ganglia"
  unique_GTEx[unique_GTEx == "Cerebellum"] <- "Brain_Cerebellum"
  unique_GTEx[unique_GTEx == "Brain_Nucleus_Accumbens_Basal_Ganglia"] <- "Brain_Nucleus_accumbens_basal_ganglia"
  unique_GTEx[unique_GTEx == "Brain_Putamen_Basal_Ganglia"] <- "Brain_Putamen_basal_ganglia"
  unique_GTEx[unique_GTEx == "Brain_Substantia_Nigra"] <- "Brain_Substantia_nigra"
  unique_GTEx[unique_GTEx == "Cells_Cultured_Fibroblasts"] <- "Cells_Cultured_fibroblasts"
  unique_GTEx[unique_GTEx == "Skeletal_Muscle"] <- "Muscle_Skeletal"
  unique_GTEx[unique_GTEx == "Skin_Sun_Exposed_Lower_Leg"] <- "Skin_Sun_Exposed_Lower_leg"
  unique_GTEx <- unique_GTEx[!unique_GTEx == "Kidney_Medulla"]
  unique_GTEx[unique_GTEx == "Prostrate"] <- "Prostate"
  if(length(unique_GTEx) > 0){
    for(j in unique_GTEx){
      print(j)
      jstudy <- readRDS(paste0("/data/abattle4/prashanthi/recount3_paper/data/all_consensus/corrected_expression/", j, ".rds"))
      jstudy <- jstudy[ ,colnames(jstudy) %in% df_GTEx$sample_id]
      data[[ind]] <- as.matrix(jstudy)
      ind <- ind + 1
    }
  }
  unique_tcga <- unique(df_TCGA$study)
  if(length(unique_tcga) > 0){
    for(j in unique_tcga){
      print(j)
      jstudy <- readRDS(paste0("/data/abattle4/prashanthi/recount3_paper/data/all_consensus/corrected_expression/", j, ".rds"))
      jstudy <- jstudy[ ,colnames(jstudy) %in% df_TCGA$sample_id]
      data[[ind]] <- as.matrix(jstudy)
      ind <- ind + 1
    }
  }
  
  common_genes <- rownames(data[[1]])
  for(p in c(2:length(data))){
    common_genes <- intersect(common_genes, rownames(data[[p]]))
  }
  
  data <- lapply(data, function(idat){
    idat <- as.matrix(idat[rownames(idat) %in% common_genes, ])
    idat <- as.matrix(idat[match(common_genes, rownames(idat)), ])
    idat
  })
  data <- do.call(cbind, data)
  data <- limma::normalizeQuantiles(data)
  data
})

names(tissue_dat) <- unique_tissues
diff_none_excl_vs_any <- c()
diff_one_vs_any <- c()
diff_two_vs_any <- c()
diff_three_vs_any <- c()
for(t in c(1:length(unique_tissues))){
  cat(t, "\n")
  cov_list <- list()
  cind <- 1
  for(thresholds in c(0, 3, 4, 5, 6)){
    select_samples <- all_samples_outlier_detection[all_samples_outlier_detection$tissue == unique_tissues[t], ]
    select_samples <- select_samples[select_samples$score >= thresholds, ]
    dat <- tissue_dat[[t]]
    dat <- dat[ ,colnames(dat) %in% select_samples$sample_id]
    cov_list[[cind]] <- covar(scale(t(dat)))
    cind <- cind + 1
  }
  diff_none_excl_vs_any[t] <- norm((cov_list[[1]] - cov_list[[5]]), type = "F")
  diff_one_vs_any[t] <- norm((cov_list[[4]] - cov_list[[5]]), type = "F")
  diff_two_vs_any[t] <- norm((cov_list[[3]] - cov_list[[5]]), type = "F")
  diff_three_vs_any[t] <- norm((cov_list[[2]] - cov_list[[5]]), type = "F")
}

outlier_forbenius <- data.frame("tissue" = unique_tissues, 
                                "None_vs_any" = diff_none_excl_vs_any, 
                                "One_vs_any" = diff_one_vs_any, 
                                "Two_vs_any" = diff_two_vs_any, 
                                "Three_vs_any" = diff_three_vs_any)

random.seed <- 100
diff_other_tissue <- c()
alt_tissue <- c()
for(t in c(1:length(unique_tissues))){
  cat(t, "\n")
  dat <- tissue_dat[[t]]
  select_samples <- all_samples_outlier_detection[all_samples_outlier_detection$tissue == unique_tissues[t], ]
  select_samples <- select_samples[select_samples$score >= 6, ]
  dat <- dat[ ,colnames(dat) %in% select_samples$sample_id]
  
  tissue_ind <- c(1:length(unique_tissues))
  tissue_ind <- tissue_ind[!tissue_ind == t]
  alt_id <- sample(tissue_ind, 1) 
  alt_tissue[t] <- unique_tissues[alt_id]
  ext_dat <- tissue_dat[[alt_id]]
  select_samples_alt <- all_samples_outlier_detection[all_samples_outlier_detection$tissue == unique_tissues[alt_id], ]
  select_samples_alt <- select_samples_alt[select_samples_alt$score >= 6, ]
  ext_dat <- ext_dat[ ,colnames(ext_dat) %in% select_samples_alt$sample_id]
  
  common_genes <- intersect(rownames(dat), rownames(ext_dat))
  dat <- dat[rownames(dat) %in% common_genes, ]
  ext_dat <- ext_dat[rownames(ext_dat) %in% common_genes, ]
  dat <- dat[match(common_genes, rownames(dat)), ]
  ext_dat <- ext_dat[match(common_genes, rownames(ext_dat)), ]
  dat_cov <- covar(scale(t(dat)))
  ext_cov <- covar(scale(t(ext_dat)))
  diff_other_tissue[t] <-  norm((dat_cov - ext_cov), type = "F")
}


outlier_forbenius$Other_tissue_vs_any <- diff_other_tissue
outlier_forbenius$Other_tissue <- alt_tissue

saveRDS(outlier_forbenius, "/data/abattle4/prashanthi/recount3_paper/results/outlier_forbenius_new.rds")
outlier_forbenius_df <- reshape2::melt(outlier_forbenius[, c(1:6)], id.vars = "tissue")
outlier_forbenius_df$variable <- as.character(outlier_forbenius_df$variable)
outlier_forbenius_df$variable[outlier_forbenius_df$variable == "None_vs_any"] <- "No outlier exclusions"
outlier_forbenius_df$variable[outlier_forbenius_df$variable == "One_vs_any"] <- "Include outliers identified by <= 1 method"
outlier_forbenius_df$variable[outlier_forbenius_df$variable == "Two_vs_any"] <- "Include outliers identified by <= 2 methods"
outlier_forbenius_df$variable[outlier_forbenius_df$variable == "Three_vs_any"] <- "Include outliers identified by <= 3 methods"
outlier_forbenius_df$variable[outlier_forbenius_df$variable == "Other_tissue_vs_any"] <- "Non-outlier samples (unrelated tissue)"

outlier_forbenius_df$variable <- factor(outlier_forbenius_df$variable, levels = c("Include outliers identified by <= 1 method",
                                                                                  "Include outliers identified by <= 2 methods", 
                                                                                  "Include outliers identified by <= 3 methods", 
                                                                                  "No outlier exclusions", 
                                                                                  "Non-outlier samples (unrelated tissue)"))

ggplot(outlier_forbenius_df, aes(x = variable, y = value, color = variable, fill = variable)) + 
  geom_boxplot(alpha=0.3, outlier.size = 0) + geom_jitter() + theme_classic() + scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + xlab("") + ylab("Forbenius norm") + guides(fill = guide_legend(title = "Covariance compared"), color = guide_legend(title = "Covariance compared")) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                                                               legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
                                                                               axis.title = element_text(size = 12, face = "bold"), legend.position = "right") + ggtitle("Forbenius norm of difference in empirical covariance")

p1 <- ggplot(outlier_forbenius_df[outlier_forbenius_df$tissue == "Blood", ], aes(x = variable, y = value, fill = variable)) + geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Dark2") + theme_classic() +  xlab("") + ylab("Forbenius norm") + guides(fill = guide_legend(title = "Covariance compared"), color = guide_legend(title = "Covariance compared")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ggtitle("Blood")+ theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                                                                                   legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
                                                                                                   axis.title = element_text(size = 12, face = "bold"), legend.position = "right")
p2 <- ggplot(outlier_forbenius_df[outlier_forbenius_df$tissue == "Central nervous system", ], aes(x = variable, y = value, fill = variable)) + geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Dark2") + theme_classic() +  xlab("") + ylab("Forbenius norm") + guides(fill = guide_legend(title = "Covariance compared"), color = guide_legend(title = "Covariance compared")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ggtitle("Central nervous system")+ theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                                                                                                  legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
                                                                                                                  axis.title = element_text(size = 12, face = "bold"), legend.position = "right")

p3 <- ggplot(outlier_forbenius_df[outlier_forbenius_df$tissue == "Liver", ], aes(x = variable, y = value, fill = variable)) + geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Dark2") + theme_classic() +  xlab("") + ylab("Forbenius norm") + guides(fill = guide_legend(title = "Covariance compared"), color = guide_legend(title = "Covariance compared")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ggtitle("Liver")+ theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                                                                                 legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
                                                                                                 axis.title = element_text(size = 12, face = "bold"), legend.position = "right")

p4 <- ggplot(outlier_forbenius_df[outlier_forbenius_df$tissue == "Skin", ], aes(x = variable, y = value, fill = variable)) + geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Dark2") + theme_classic() +  xlab("") + ylab("Forbenius norm") + guides(fill = guide_legend(title = "Covariance compared"), color = guide_legend(title = "Covariance compared")) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ggtitle("Skin")+ theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                                                                                 legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
                                                                                                 axis.title = element_text(size = 12, face = "bold"), legend.position = "right")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend_sample <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
plot_grid(plot_grid(p1, p2, p3,p4, nrow = 2), legend_sample, ncol= 2, rel_widths = c(1, 0.8))
