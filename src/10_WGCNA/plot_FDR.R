rm(list = ls())
library(ggplot2)

net <- "GTEx/consensus"
#agg_levels <- c("net_1", "net_20", "net_40", "net_60", "net_65")
agg_levels <- c("net_1", "net_10", "net_20", "net_30", "net_100", "net_300", "net_566")
#agg_levels <- c("net_1", "net_3", "net_10", "net_30", "net_49")
res_list <- list()
for(i in c(1:length(agg_levels))){
  res_list[[i]] <- readRDS(paste0("/data/abattle4/prashanthi/recount3_paper/results/WGCNA/networks/", net ,"/", agg_levels[i], 
                                  "/fdr_results.rds"))
  res_list[[i]]$nStudies <- agg_levels[i]
}

res_df <- do.call(rbind, res_list)
res_df$nStudies <- factor(res_df$nStudies, levels = agg_levels)

pal.colors <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Blues"))(9)[-c(1:4)]

res_df$F1_score <- (2*res_df$precision*res_df$recall)/(res_df$precision + res_df$recall)

res_df$nStudies <- gsub("net_", "", res_df$nStudies)
res_df$nStudies <- as.numeric(res_df$nStudies)
#res_df$nStudies <- factor(res_df$nStudies, levels = c(1, 10, 20, 30, 100, 300, 566))
res_df$nStudies <- factor(res_df$nStudies, levels = c(1, 20, 40, 60, 65))
#res_df$nStudies <- factor(res_df$nStudies, levels = c(1, 3, 10, 30, 49))

p1 <- ggplot(res_df, aes(x = nEdges, y = precision, color = nStudies)) + geom_point(size = 1) + theme_classic() +
  scale_colour_manual(name = "Number of studies\n aggregated", values = pal.colors) + xlab("Number of edges") +geom_line()+
  ylab("Precision") + ggtitle("Blood") + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                       legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
                                       axis.title = element_text(size = 12, face = "bold"), legend.position = "right")  + xlim(c(5000, 60000)) + ylim(c(0, 0.22))


p2 <- ggplot(res_df, aes(x = nEdges, y = F1_score, color = nStudies)) + geom_point(size = 1) + theme_classic() +
  scale_colour_manual(name = "Number of studies\n aggregated", values = pal.colors) + xlab("Number of edges") +geom_line()+
  ylab("F1 Score") + ggtitle("Blood") + theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
                                                                legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
                                                                axis.title = element_text(size = 12, face = "bold"), legend.position = "right")  + xlim(c(5000, 60000))+ ylim(c(0, 0.01))


