# Description
# This script reads in the t-SNE dimensionality reduction for samples belonging to the immune system
# found in data/tissue_xcell_tsne obtained by 09-xCell_deconvolution_and_dimRed.R
# And for each immune subcluster, identify outliers using Z score, MAD, and the Mahalanobis distance
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(ruler)

#Inputs
# 1. data/tissue_xcell_tsne/immune.rds

#Outputs
# 1. data/immune_samples_outlier_detectection.rds
# 2. data/filter_two_fail_immune.rds




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

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
datDir <- paste0(homeDir, "data/tissue_xcell_tsne/")
df <- readRDS(paste0(datDir, "immune.rds"))

plotDir <- paste0(homeDir, "plots/immune_preprocessing/")
pdf(paste0(plotDir, "all_samples.pdf"))
par(mfrow = c(2, 2))
for(i in unique(df$tissue_subtype)){
  plot(df$D1[!df$tissue_subtype == i], df$D2[!df$tissue_subtype == i], main = i, 
       xlab = "D1", ylab = "D2", pch = 20, col = alpha("grey", 0.1), xlim = c(-50, 50), ylim = c(-50, 50))
  points(df$D1[df$tissue_subtype == i], df$D2[df$tissue_subtype == i],
         xlab = "D1", ylab = "D2", pch = 20, col = alpha("#CC66FF", 0.6))
}

plt.colors <-c("blue", "red", "green")
set.seed(123)
km.res <- kmeans(df[ ,1:2], 3, nstart = 25)
plot(df$D1, df$D2, main = "Clustering", 
     xlab = "D1", ylab = "D2", pch = 20, col = alpha(plt.colors[km.res$cluster], 0.4), xlim = c(-50, 50), ylim = c(-50, 50))
legend("topright", legend = c("1", "2", "3"), fill = plt.colors)
dev.off()

cluster_1 <- df[km.res$cluster == 1, ]
cluster_2 <- df[km.res$cluster == 2, ]
cluster_3 <- df[km.res$cluster == 3, ]

find_outliers <- function(df, title){
  col_outlier <- df %>% transmute_if(is.numeric, isnt_out_funs)
  md_outlier <- df %>% 
    transmute(maha = maha_dist(.)) %>%
    transmute_at(vars(maha = maha), isnt_out_funs)
  plotDir <- paste0(homeDir, "plots/immune_preprocessing/")
  pdf(paste0(plotDir, title, ".pdf"))
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

cluster_1 <- find_outliers(cluster_1, "cluster_1")
cluster_2 <- find_outliers(cluster_2, "cluster_2")
cluster_3 <- find_outliers(cluster_3, "cluster_3")

outlier_annot <- list(cluster_1, cluster_2, cluster_3)
df <- rbind(cluster_1, cluster_2, cluster_3)
saveRDS(df, paste0(homeDir, "data/immune_samples_outlier_detectection.rds"))

filter_two_fail <- lapply(outlier_annot, function(df){
  df$z <- df$D1_z & df$D2_z
  df$mad <- df$D1_mad & df$D2_mad
  df$tukey <- df$D1_tukey & df$D2_tukey
  df$score <- rowSums(df[ ,15:20])
  df$outlier <- df$score > 4
  df[df$outlier, ]
})

cluster_id <- c(rep("1", dim(filter_two_fail[[1]])[1]), 
                rep("2", dim(filter_two_fail[[2]])[1]),
                rep("3", dim(filter_two_fail[[3]])[1]))
filter_two_fail <- data.frame(do.call(rbind, filter_two_fail))
filter_two_fail$cluster_id <- cluster_id

saveRDS(filter_two_fail, 
        paste0(homeDir, "data/filter_two_fail_immune.rds"))