# Description
# This script reads in the t-SNE dimensionality reduction for each tissue found in 
# data/tissue_xcell_tsne obtained by 09-xCell_deconvolution_and_dimRed.R
# And for each tissue, identify outliers using Z score, MAD, and the Mahalanobis distance
rm(list = ls())
#Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(ruler)
library(viridis)

#Inputs
# 1. data/tissue_xcell_tsne

#Outputs
# 1. data/all_samples_outlier_detection.rds
# 2. data/filter_mahalanobis_mad.rds
# 3. data/filter_mahalanobis_tukey.rds
# 4. data/filter_two_fail.rds

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
plot_single_tissues <- function(filename, df){
  plotDir <- paste0(homeDir, "plots/tissue_outlier/")
  pdf(paste0(plotDir, filename))
  par(mfrow = c(3, 3))
  plot(df$D1[df$tissue == "Adipose"], df$D2[df$tissue == "Adipose"],
       col = alpha(df$col_tissue[df$tissue == "Adipose"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Adipose", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Airway"], df$D2[df$tissue == "Airway"],
       col = alpha(df$col_tissue[df$tissue == "Airway"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Airway", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Blood"], df$D2[df$tissue == "Blood"],
       col = alpha(df$col_tissue[df$tissue == "Blood"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Blood", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Breast"], df$D2[df$tissue == "Breast"],
       col = alpha(df$col_tissue[df$tissue == "Breast"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Breast", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Cardiac"], df$D2[df$tissue == "Cardiac"],
       col = alpha(df$col_tissue[df$tissue == "Cardiac"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Cardiac", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Central nervous system"], df$D2[df$tissue == "Central nervous system"],
       col = alpha(df$col_tissue[df$tissue == "Central nervous system"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "CNS", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Colon"], df$D2[df$tissue == "Colon"],
       col = alpha(df$col_tissue[df$tissue == "Colon"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Colon", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Esophagus"], df$D2[df$tissue == "Esophagus"],
       col = alpha(df$col_tissue[df$tissue == "Esophagus"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Esophagus", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Eye"], df$D2[df$tissue == "Eye"],
       col = alpha(df$col_tissue[df$tissue == "Eye"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Eye", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Fibroblasts"], df$D2[df$tissue == "Fibroblasts"],
       col = alpha(df$col_tissue[df$tissue == "Fibroblasts"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Fibroblasts", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "hESCs"], df$D2[df$tissue == "hESCs"],
       col = alpha(df$col_tissue[df$tissue == "hESCs"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "hESCs", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Intestine"], df$D2[df$tissue == "Intestine"],
       col = alpha(df$col_tissue[df$tissue == "Intestine"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Intestine", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "iPSCs"], df$D2[df$tissue == "iPSCs"],
       col = alpha(df$col_tissue[df$tissue == "iPSCs"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "iPSCs", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Kidney"], df$D2[df$tissue == "Kidney"],
       col = alpha(df$col_tissue[df$tissue == "Kidney"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Kidney", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Liver"], df$D2[df$tissue == "Liver"],
       col = alpha(df$col_tissue[df$tissue == "Liver"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Liver", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Lung"], df$D2[df$tissue == "Lung"],
       col = alpha(df$col_tissue[df$tissue == "Lung"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Lung", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Multipotent cells"], df$D2[df$tissue == "Multipotent cells"],
       col = alpha(df$col_tissue[df$tissue == "Multipotent cells"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Multipotent cells", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue_subtype == "skeletal muscle"], df$D2[df$tissue_subtype == "skeletal muscle"],
       col = alpha(df$col_tissue[df$tissue_subtype == "skeletal muscle"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Skeletal muscle", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Nervous system"], df$D2[df$tissue == "Nervous system"],
       col = alpha(df$col_tissue[df$tissue == "Nervous system"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Nervous system", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Pancreas"], df$D2[df$tissue == "Pancreas"],
       col = alpha(df$col_tissue[df$tissue == "Pancreas"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Pancreas", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Prostrate"], df$D2[df$tissue == "Prostrate"],
       col = alpha(df$col_tissue[df$tissue == "Prostrate"], 0.9), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Prostrate", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Skin"], df$D2[df$tissue == "Skin"],
       col = alpha(df$col_tissue[df$tissue == "Skin"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Skin", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Stomach"], df$D2[df$tissue == "Stomach"],
       col = alpha(df$col_tissue[df$tissue == "Stomach"], 0.9), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Stomach", xlim = c(-50, 50), ylim = c(-50, 50))
  plot(df$D1[df$tissue == "Vascular"], df$D2[df$tissue == "Vascular"],
       col = alpha(df$col_tissue[df$tissue == "Vascular"], 0.5), pch = 20, xlab = "Dimension 1", ylab = "Dimension 2", main = "Vascular", xlim = c(-50, 50), ylim = c(-50, 50))
  dev.off()
}

isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}

isnt_out_mad <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - median(x, na.rm = na.rm)) <= thres * mad(x, na.rm = na.rm)
}

isnt_out_tukey <- function(x, k = 1.5, na.rm = TRUE) {
  quar <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- diff(quar)
  
  (quar[1] - k * iqr <= x) & (x <= quar[2] + k * iqr)
}

maha_dist <- . %>% select_if(is.numeric) %>%
  mahalanobis(center = colMeans(.), cov = cov(.))
isnt_out_maha <- function(tbl, isnt_out_f, ...) {
  tbl %>% maha_dist() %>% isnt_out_f(...)
}

isnt_out_funs <- list(
  z = isnt_out_z,
  mad = isnt_out_mad,
  tukey = isnt_out_tukey
)

find_outliers <- function(tissue){
  df <- readRDS(paste0(datDir, tissue))
  col_outlier <- df %>% transmute_if(is.numeric, isnt_out_funs)
  md_outlier <- df %>% 
    transmute(maha = maha_dist(.)) %>%
    transmute_at(vars(maha = maha), isnt_out_funs)
  plotDir <- paste0(homeDir, "plots/tissue_outlier/")
  pdf(paste0(plotDir, gsub(".rds", ".pdf", tissue)))
  par(mfrow = c(3, 2))
  par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
  plot(df$D1, df$D2, col = alpha(c("red", "blue")[(col_outlier$D1_z & col_outlier$D2_z) + 1], 0.5), main = "z-score", pch = 20, xlab = "D1", ylab = "D2", xlim = c(-50, 50), ylim = c(-50, 50))
  legend("topright", inset=c(-0.4,0), legend = c("outlier", "not outlier"), col = alpha(c("red", "blue"), 0.5), pch = 20, box.lty=0)
  par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
  plot(df$D1, df$D2, col = alpha(c("red", "blue")[(col_outlier$D1_mad & col_outlier$D2_mad) + 1], 0.5), main = "MAD", pch = 20, xlab = "D1", ylab = "D2", xlim = c(-50, 50), ylim = c(-50, 50))
  legend("topright", inset=c(-0.4,0), legend = c("outlier", "not outlier"), col = alpha(c("red", "blue"), 0.5), pch = 20, box.lty=0)
  par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
  plot(df$D1, df$D2, col = alpha(c("red", "blue")[(col_outlier$D1_tukey & col_outlier$D2_tukey) + 1], 0.5), main = "Tukey", pch = 20, xlab = "D1", ylab = "D2", xlim = c(-50, 50), ylim = c(-50, 50))
  legend("topright", inset=c(-0.4,0), legend = c("outlier", "not outlier"), col = alpha(c("red", "blue"), 0.5), pch = 20, box.lty=0)
  par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
  plot(df$D1, df$D2, col = alpha(c("red", "blue")[md_outlier$maha_z + 1], 0.5), main = "Mahalanobis z-score", pch = 20, xlab = "D1", ylab = "D2", xlim = c(-50, 50), ylim = c(-50, 50))
  legend("topright", inset=c(-0.4,0), legend = c("outlier", "not outlier"), col = alpha(c("red", "blue"), 0.5), pch = 20, box.lty=0)
  par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
  plot(df$D1, df$D2, col = alpha(c("red", "blue")[md_outlier$maha_mad + 1], 0.5), main = "Mahalanobis MAD", pch = 20, xlab = "D1", ylab = "D2", xlim = c(-50, 50), ylim = c(-50, 50))
  legend("topright", inset=c(-0.4,0), legend = c("outlier", "not outlier"), col = alpha(c("red", "blue"), 0.5), pch = 20, box.lty=0)
  par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
  plot(df$D1, df$D2, col = alpha(c("red", "blue")[md_outlier$maha_tukey + 1], 0.5), main = "Mahalanobis Tukey", pch = 20, xlab = "D1", ylab = "D2", xlim = c(-50, 50), ylim = c(-50, 50))
  legend("topright", inset=c(-0.4,0), legend = c("outlier", "not outlier"), col = alpha(c("red", "blue"), 0.5), pch = 20, box.lty=0)
  dev.off()
  df <- cbind(df, col_outlier, md_outlier)
  df
}

datDir <- paste0(homeDir, "data/tissue_xcell_tsne/")
tissues <- list.files(datDir)
tissues <- tissues[!tissues %in% c("musculoskeletal.rds", "immune.rds")]

outlier_annot <- lapply(tissues, function(tissue){
  find_outliers(tissue)
})

names(outlier_annot) <- gsub(".rds", "", tissues)

filter_mahalanobis_mad <- lapply(outlier_annot, function(df){
  df <- df[df$maha_mad, ]
})

filter_mahalanobis_tukey <- lapply(outlier_annot, function(df){
  df <- df[df$maha_tukey, ]
})

filter_two_fail <- lapply(outlier_annot, function(df){
  df$z <- df$D1_z & df$D2_z
  df$mad <- df$D1_mad & df$D2_mad
  df$tukey <- df$D1_tukey & df$D2_tukey
  df$score <- rowSums(df[ ,15:20])
  df$outlier <- df$score > 4
  df[df$outlier, ]
})

all_samples <- data.frame(do.call(rbind, outlier_annot))
rownames(all_samples) <- c(1:dim(all_samples)[1])

saveRDS(all_samples, paste0(homeDir, "data/all_samples_outlier_detection.rds"))

filter_mahalanobis_mad <- data.frame(do.call(rbind, filter_mahalanobis_mad))
rownames(filter_mahalanobis_mad) <- c(1:dim(filter_mahalanobis_mad)[1])

filter_mahalanobis_tukey <- data.frame(do.call(rbind, filter_mahalanobis_tukey))
rownames(filter_mahalanobis_tukey) <- c(1:dim(filter_mahalanobis_tukey)[1])

filter_two_fail <- data.frame(do.call(rbind, filter_two_fail))
rownames(filter_two_fail) <- c(1:dim(filter_two_fail)[1])

plot_single_tissues("all_tissues_maha_mad.pdf", filter_mahalanobis_mad)
plot_single_tissues("all_tissues_maha_tukey.pdf", filter_mahalanobis_tukey)
plot_single_tissues("all_tissues_two_fail.pdf", filter_two_fail)


summary_samples <- data.frame(x = as.numeric(table(all_samples$tissue)), 
                              Tissue = names(table(all_samples$tissue)), 
                              Method = "Total samples")
summary_samples <- rbind(summary_samples, 
                         data.frame(x = as.numeric(table(filter_mahalanobis_mad$tissue)), 
                                    Tissue = names(table(filter_mahalanobis_mad$tissue)), 
                                    Method = "w/o outliers Mahalanobis MAD"))
summary_samples <- rbind(summary_samples, 
                         data.frame(x = as.numeric(table(filter_mahalanobis_tukey$tissue)), 
                                    Tissue = names(table(filter_mahalanobis_tukey$tissue)), 
                                    Method = "w/o outliers Mahalanobis Tukey"))
summary_samples <- rbind(summary_samples, 
                         data.frame(x = as.numeric(table(filter_two_fail$tissue)), 
                                    Tissue = names(table(filter_two_fail$tissue)), 
                                    Method = "w/o outliers based on 2 metrics"))

plotDir <- paste0(homeDir, "plots/tissue_outlier/")
pdf(paste0(plotDir, "number_samples.pdf"))
ggplot(summary_samples, aes(x = Tissue, y = x, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") + theme_classic() + 
  coord_flip() + scale_fill_viridis(discrete = TRUE) + scale_x_discrete(limits = names(table(all_samples$tissue)[order(table(all_samples$tissue))]))
dev.off()

saveRDS(filter_mahalanobis_mad, paste0(homeDir, "data/filter_mahalanobis_mad.rds"))
saveRDS(filter_mahalanobis_tukey, paste0(homeDir, "data/filter_mahalanobis_tukey.rds"))
saveRDS(filter_two_fail, paste0(homeDir, "data/filter_two_fail.rds"))


