library(ggplot2)
rm(list = ls())
resDir <- "/data/abattle4/prashanthi/recount3_paper/results/weighted_cov_networks/"
context <- "skin/all/net_20"

res_odds_perm <- readRDS(paste0(resDir, context, "/res_odds_perm_2.rds"))
lambda <- res_odds_perm$lambda
odds_actual <- res_odds_perm$odds_actual
odds_perm <- data.frame(res_odds_perm$odds_perm)
summary_df <- data.frame("Lambda" = lambda, "Odds_ratio" = odds_actual)
summary_df <- cbind(summary_df, odds_perm)
summary_df <- reshape2::melt(summary_df, id.vars = "Lambda")
colnames(summary_df) <- c("Lambda", "Type", "Odds")
summary_df$Type <- as.character(summary_df$Type)
summary_df$Type[grep("perm", summary_df$Type)] <- "Permuted"
summary_df$Type[summary_df$Type == "Odds_ratio"] <- "Actual"
summary_df$Type <- factor(summary_df$Type, levels = c("Actual", "Permuted"))

ggplot(summary_df, aes(x = Lambda, y = Odds, fill = Type, color = Type)) + geom_boxplot(outlier.shape = 21, outlier.size = 0.5) + 
  geom_hline(yintercept = 1, linetype = "dashed") + theme_classic() + ggtitle("Skin") +
  theme(plot.title = element_text(face="bold", size = 16), axis.text = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"), legend.position = "bottom") + ylab("Odds ratio")


