rm(list = ls())
library(recount3)
# This script reads in all the count data use in the universal consensus networks 
# and prepares the corresponding metadata by combining the labels available in GTEx and TCGA
# with the manually curated labels for SRA

# Input
# 1. RSE objects from data/all_consensus/uncorrected_expression obtained from 05_process_expr_all_consensus.R
# 2. data/manual_annotations.tsv

homeDir <- "/data/abattle4/prashanthi/recount3_paper/"
uncorrectedDir <- paste0(homeDir, "data/all_consensus/uncorrected_expression/")

read_study <- function(study){
  projects <- list.files(paste0(uncorrectedDir, study, "/"))
  rse_list <- lapply(projects, function(iproj){
    readRDS(paste0(uncorrectedDir, study, "/", iproj))
  })
  names(rse_list) <- gsub(".rds", "", projects)
  rse_list
}


studies <- c("GTEx", "sra", "tcga")

sra <- read_study("sra")

sra[["SRP096986"]] <- NULL # This is a single cell study
sra[["SRP135684"]] <- NULL # this is a single cell study
sra[["SRP166966"]] <- NULL # this is a single cell study
sra[["SRP200058"]] <- NULL # this is a single cell study
sra[["SRP067502"]] <- NULL # this is a single cell study
nsra <- lapply(sra, function(rse){
  dim(rse)[2]
})

nsra <- unlist(nsra)
sra <- sra[nsra >= 20]

gtex <- read_study("GTEx")
tcga <- read_study("tcga")

sra_counts <- lapply(sra, function(rse){
  assays(rse)$RPKM
})
gtex_counts <- lapply(gtex, function(rse){
  assays(rse)$RPKM
})
tcga_counts <- lapply(tcga, function(rse){
  assays(rse)$RPKM
})

all_counts <- c(sra_counts, gtex_counts, tcga_counts)
rm(sra_counts)
rm(gtex_counts)
rm(tcga_counts)

common.genes <- rownames(all_counts[[1]])
for(i in c(2:length(all_counts))){
  common.genes <- intersect(common.genes, rownames(all_counts[[i]]))
}

all_counts <- lapply(all_counts, function(counts){
  counts <- counts[rownames(counts) %in% common.genes, ]
  counts <- counts[match(common.genes, rownames(counts)), ]
  counts
})

counts <- all_counts[[1]]
for(i in c(2:length(all_counts))){
  counts <- cbind(counts, all_counts[[i]])
}

samples <- colnames(counts)

gtex_metadata <- cbind(gtex[[1]]$external_id, gtex[[1]]$gtex.smtsd)
for(i in c(2:length(gtex))){
  df <- cbind(gtex[[i]]$external_id, gtex[[i]]$gtex.smtsd)
  gtex_metadata <- rbind(gtex_metadata, df)
}
gtex_metadata <- data.frame(gtex_metadata)
colnames(gtex_metadata) <- c("sample", "tissue.category")
gtex_metadata$study <- "GTEx"
gtex_metadata$cancer <- "none"

gtex_metadata <- gtex_metadata[ , c("study", "sample", "tissue.category", "cancer")]

tcga_metadata <- cbind(tcga[[1]]$study ,colnames(tcga[[1]]), tcga[[1]]$tcga.gdc_cases.project.primary_site, tcga[[1]]$tcga.cgc_sample_sample_type)
for(i in c(2:length(tcga))){
  df <- cbind(tcga[[i]]$study ,colnames(tcga[[i]]), tcga[[i]]$tcga.gdc_cases.project.primary_site, tcga[[i]]$tcga.cgc_sample_sample_type)
  tcga_metadata <- rbind(tcga_metadata, df)
}

tcga_metadata <- data.frame(tcga_metadata)
colnames(tcga_metadata) <- c("study", "sample", "tissue.category", "cancer")
tcga_metadata$cancer[tcga_metadata$cancer == "Solid Tissue Normal"] <- "none"
tcga_metadata$cancer[!tcga_metadata$cancer == "none"] <- "cancer"

tissue_df <- read.delim(homeDir, "/data/manual_annotations.tsv")
tissue_df <- tissue_df[ ,c("study", "sample", "tissue.category", "cancer")]
tissue_df$cancer[tissue_df$cancer == "normal"] <- "none"
tissue_df$cancer[tissue_df$study == "SRP093727"] <- "none"
tissue_df$cancer[tissue_df$study == "SRP126741"] <- "none"
tissue_df$cancer[tissue_df$study == "SRP097611"] <- "cancer"
tissue_df <- tissue_df[!tissue_df$study == "SRP067502", ]
all_metadata <- rbind(tissue_df, gtex_metadata, tcga_metadata)
source <- c(rep("sra", dim(tissue_df)[1]), rep("GTEx", dim(gtex_metadata)[1]), rep("tcga", dim(tcga_metadata)[1]))
all_metadata <- cbind(all_metadata, source)
all_metadata <- all_metadata[match(colnames(counts), all_metadata$sample), ]

saveRDS(all_metadata, paste0(homeDir, "data/all_metadata.rds"))
saveRDS(counts, paste0(homeDir, "data/all_counts.rds"))

