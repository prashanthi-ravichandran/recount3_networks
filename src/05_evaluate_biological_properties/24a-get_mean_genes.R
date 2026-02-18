rm(list = ls())

library(Matrix)
library(dplyr)
datDir <- "/data/abattle4/prashanthi/recount3_paper/data/"
resDir <- "/data/abattle4/prashanthi/recount3_paper/results/"
net_name <- "skin/GTEx/"

if(net_name == "all_consensus"){
  sources <- c("GTEx", "sra", "tcga")
  expr_all <- lapply(sources, function(isource){
    studies <- list.files(paste0(datDir, net_name, "/uncorrected_expression/", isource, "/"))
    expr <- lapply(studies, function(istudy){
      idat <- readRDS(paste0(datDir, net_name, "/uncorrected_expression/", isource,  "/", istudy))
      idat@assays@data$RPKM
    })
    studies <- gsub(".rds", "", studies)
    names(expr) <- studies 
    nSamples <- lapply(expr, function(iexpr){
      dim(iexpr)[2]
    })
    nSamples <- unlist(nSamples)
    expr <- expr[nSamples > 15]
    if(length(expr) > 1){
      common.genes <- rownames(expr[[1]])
      for(i in c(2:length(expr))){
        common.genes <- intersect(common.genes, rownames(expr[[i]]))
      }
      expr <- lapply(expr, function(iexpr){
        iexpr[rownames(iexpr) %in% common.genes, ]
      })
      expr_mat <- do.call(cbind, expr)
    }else{expr_mat <- expr[[1]]}
    
    expr_mat })
    common.genes <- rownames(expr_all[[1]])
    for(i in c(2:length(expr_all))){
      common.genes <- intersect(common.genes, rownames(expr_all[[i]]))
    }
    expr_all <- lapply(expr_all, function(iexpr){
      iexpr <- iexpr[rownames(iexpr) %in% common.genes, ]
      iexpr <- iexpr[match(common.genes, rownames(iexpr)), ]
      iexpr
    })
    expr_mat <- do.call(cbind, expr_all)
}else{
  if(net_name == "normal_consensus"){
    sources <- c("GTEx", "sra")
    expr_all <- lapply(sources, function(isource){
      studies <- list.files(paste0(datDir, net_name, "/uncorrected_expression/", isource, "/"))
      expr <- lapply(studies, function(istudy){
        idat <- readRDS(paste0(datDir, net_name, "/uncorrected_expression/", isource,  "/", istudy))
        idat@assays@data$RPKM
      })
      studies <- gsub(".rds", "", studies)
      names(expr) <- studies 
      nSamples <- lapply(expr, function(iexpr){
        dim(iexpr)[2]
      })
      nSamples <- unlist(nSamples)
      expr <- expr[nSamples > 15]
      
      if(length(expr) > 1){
        common.genes <- rownames(expr[[1]])
        for(i in c(2:length(expr))){
          common.genes <- intersect(common.genes, rownames(expr[[i]]))
        }
        expr <- lapply(expr, function(iexpr){
          iexpr[rownames(iexpr) %in% common.genes, ]
        })
        expr_mat <- do.call(cbind, expr)
      }else{
        expr_mat <- expr[[1]]}
      expr_mat
    })
    common.genes <- rownames(expr_all[[1]])
    for(i in c(2:length(expr_all))){
      common.genes <- intersect(common.genes, rownames(expr_all[[i]]))
    }
    expr_all <- lapply(expr_all, function(iexpr){
      iexpr <- iexpr[rownames(iexpr) %in% common.genes, ]
      iexpr <- iexpr[match(common.genes, rownames(iexpr)), ]
      iexpr
    })
    expr_mat <- do.call(cbind, expr_all)
}else{
    studies <- list.files(paste0(datDir, net_name, "uncorrected_expression/"))
    expr <- lapply(studies, function(istudy){
      idat <- readRDS(paste0(datDir, net_name, "uncorrected_expression/", istudy))
      idat@assays@data$RPKM
    })
    
    studies <- gsub(".rds", "", studies)
    names(expr) <- studies 
    nSamples <- lapply(expr, function(iexpr){
      dim(iexpr)[2]
    })
    nSamples <- unlist(nSamples)
    expr <- expr[nSamples > 15]
    
    if(length(expr) > 1){
      common.genes <- rownames(expr[[1]])
      for(i in c(2:length(expr))){
        common.genes <- intersect(common.genes, rownames(expr[[i]]))
      }
      expr <- lapply(expr, function(iexpr){
        iexpr <- iexpr[rownames(iexpr) %in% common.genes, ]
        iexpr[match(common.genes, rownames(iexpr)), ]
      })
      expr_mat <- do.call(cbind, expr)
    }else{
      expr_mat <- expr[[1]]
    }
  }
}

expr_mean <- rowMeans(expr_mat)

expr_mean <- data.frame(names(expr_mean) , expr_mean)
row.names(expr_mean) <- c(1:dim(expr_mat)[1])
colnames(expr_mean) <- c("gene_id", "expr")

geneData <-  readRDS(paste0(datDir, "geneData.rds"))
symbol <- c()
for(i in c(1:dim(expr_mean)[1])){
  symbol[i] <- geneData$gene_name[geneData$gene_id == expr_mean$gene_id[i]]
}

expr_mean$gene_name <- symbol
expr_mean$quantiles <- ntile(expr_mean$expr, 5)

saveRDS(expr_mean, paste0(resDir, "mean_expr/skin_GTEx.rds"))
