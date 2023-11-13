rm(list = ls())
library(data.table)
# Read frequency of all SNPs
homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
datDir <- paste0(homeDir, "data/s_LDSC/")
resDir <- paste0(homeDir, "results/s_LDSC/")

SNP_list <- lapply(c(1:22), function(i){
  df <- read.table(paste0(datDir, "1000G_Phase3_frq/1000G.EUR.QC.", i, ".frq"), header = T)
  df <- df[df$MAF >= 0.05, ]
  df$SNP
})
SNP_list <- unlist(SNP_list)

inputArgs <- commandArgs(TRUE)
network <- inputArgs[1]
centrality <- inputArgs[2]
control <- inputArgs[3]

if(control == "baseline"){
  annot1 <- lapply(c(1:22), function(i){
    df <- fread(paste0(datDir, "all-genes/SAMGENIC.", i, ".annot.gz"))
    df <- df[df$SNP %in% SNP_list, ]
    df[ ,5]
  })
  annot1 <- do.call(rbind, annot1)
  
  annot2 <- lapply(c(1:22), function(i){
    df <- fread(paste0(datDir, "1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.", i, ".annot.gz"))
    df <- df[df$SNP %in% SNP_list, ]
    df[ ,5:101]
  })
  annot2 <- do.call(rbind, annot2)
  annot <- cbind(annot1, annot2)
  sd_control <- apply(annot, 2, sd)
}else{
  annot <- lapply(c(1:22), function(i){
    df <- fread(paste0(datDir, "all-genes/SAMGENIC.", i, ".annot.gz"))
    df <- df[df$SNP %in% SNP_list, ]
    df[ ,5]
  })
  annot <- do.call(rbind, annot)
  sd_control <- apply(annot, 2, sd)
}

annot_net <- lapply(c(1:22), function(i){
  df <- fread(paste0(resDir, 
                     network, "/", centrality, "/ldscore/", network, ".", i, ".annot.gz"))
  df <- df[df$SNP %in% SNP_list, ]
  df[ ,5]
})

annot_net <- do.call(rbind, annot_net)
sd_net <- apply(annot_net, 2, sd)

sd <- c(sd_control, sd_net)

if(control == "baseline"){
  write.table(sd, paste0(resDir, network, "/", centrality, "/ldscore/", network, "_", control, ".sd"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}else{
  write.table(sd, paste0(resDir, network, "/", centrality, "/ldscore/", network, "_", control, ".sd"), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}



