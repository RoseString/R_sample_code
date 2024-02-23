#!/usr/bin/env Rscript
# By Dan (Meredith) Sun
# Created: 01/27/2020
# Last updated: 11/28/2023
library(DESeq2)
library(biomaRt)
library(ggpubr)
library(openxlsx)
library(BiocParallel)
library(stringr)
library(apeglm)
#library(IHW)


### Methods
## A Method to fetch gene names and gene descriptions from biomaRT
biomaRT_convertID <- function(geneIDs, species) {
  print('Get gene names and descriptions for Ensembl IDs...')

  if (species == 'mouse') {
    dataset <- 'mmusculus_gene_ensembl'
  } else if (species == 'human') {
    dataset <- 'hsapiens_gene_ensembl'
  } else if (species == 'rat') {
    dataset <- 'rnorvegicus_gene_ensembl'
  } else if (species == 'yeast') {
    dataset <- 'scerevisiae_gene_ensembl'
  } else if (species == 'celegans') {
    dataset <- 'celegans_gene_ensembl'
  } else if (species == 'drosophila') {
    dataset <- 'dmelanogaster_gene_ensembl'
  } else if (species == 'frog') {
    dataset <- 'xtropicalis_gene_ensembl'
  } else if (species == 'zebrafish') {
    dataset <- 'drerio_gene_ensembl'
  }

  if (species == 'human') {
    converted.df <- tryCatch(
      {
        print("Attempting useast mirror...")
        mart <- useDataset(dataset, useMart('ENSEMBL_MART_ENSEMBL', host = 'https://useast.ensembl.org'))
        getBM(filters='ensembl_gene_id',
              attributes=c('ensembl_gene_id', 'hgnc_symbol','description','chromosome_name','start_position','end_position','gene_biotype'),
              values=geneIDs,
              mart=mart)
      },
      error = function(e) {
        tryCatch({
          print("Failed. Attempting asia mirror...")
          mart <- useDataset(dataset, useMart('ENSEMBL_MART_ENSEMBL', host = 'https://asia.ensembl.org'))
          getBM(filters='ensembl_gene_id',
                attributes=c('ensembl_gene_id', 'hgnc_symbol','description','chromosome_name','start_position','end_position','gene_biotype'),
                values=geneIDs,
                mart=mart)
        },
        error = function(e2) {
          print("Failed. Attempting global mirror...")
          mart <- useDataset(dataset, useMart('ENSEMBL_MART_ENSEMBL', host = 'https://www.ensembl.org'))
          getBM(filters='ensembl_gene_id',
                attributes=c('ensembl_gene_id', 'hgnc_symbol','description','chromosome_name','start_position','end_position','gene_biotype'),
                values=geneIDs,
                mart=mart)
        })
      }
    )
  } else {
    converted.df <- tryCatch(
      {
        print("Attempting useast mirror...")
        mart <- useDataset(dataset, useMart('ENSEMBL_MART_ENSEMBL', host = 'https://useast.ensembl.org'))
        getBM(filters='ensembl_gene_id',
              attributes=c('ensembl_gene_id', 'external_gene_name','description','chromosome_name','start_position','end_position','gene_biotype'),
              values=geneIDs,
              mart=mart)
      },
      error = function(e) {
        tryCatch({
          print("Failed. Attempting uswest mirror...")
          mart <- useDataset(dataset, useMart('ENSEMBL_MART_ENSEMBL', host = 'https://uswest.ensembl.org'))
          getBM(filters='ensembl_gene_id',
                attributes=c('ensembl_gene_id', 'external_gene_name','description','chromosome_name','start_position','end_position','gene_biotype'),
                values=geneIDs,
                mart=mart)
        },
        error = function(e2) {
          print("Failed. Attempting asia mirror...")
          mart <- useDataset(dataset, useMart('ENSEMBL_MART_ENSEMBL', host = 'https://asia.ensembl.org'))
          getBM(filters='ensembl_gene_id',
                attributes=c('ensembl_gene_id', 'external_gene_name','description','chromosome_name','start_position','end_position','gene_biotype'),
                values=geneIDs,
                mart=mart)
        })
      }
    )
  }
  return(converted.df)
}

## A method to do differential expression analysis when there is only one condition
deseq_1way <- function(samples, counts, gp1, gp2, cond, species, biomaRT_IDs = NULL, adjustVars = NULL, sizeFactors = NULL) {
  print(paste(gp2, 'vs', gp1))
  coldata <- samples[samples[[cond]] %in% c(gp1, gp2), , drop=F]
  coldata <- coldata[rownames(coldata) %in% colnames(counts), ,drop=F] # if not every sample in the condition is in the count data
  cts <- counts[, rownames(coldata)]
  coldata[[cond]] <- droplevels(factor(coldata[[cond]]))
  for (adjustVar in adjustVars) {
    coldata[[adjustVar]] <- droplevels(factor(coldata[[adjustVar]]))
    if (length(levels(coldata[[adjustVar]])) < 2) {
      adjustVars = adjustVars[adjustVars != adjustVar]
    }
  }
  print(coldata)
  if (is.null(adjustVars) | length(adjustVars) == 0) {
    fml <- paste('~', cond)
  } else {
    fml <- paste('~', paste(adjustVars,collapse=" + "), '+', cond)
  }
  print(paste('Model design:', fml))
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = formula(fml))
  dds[[cond]] <- factor(dds[[cond]], levels=c(gp1, gp2))
  if (is.null(sizeFactors)) {
    dds <- estimateSizeFactors(dds)
  } else {
    SFs <- sizeFactors[names(sizeFactors) %in% rownames(coldata)]
    print(SFs)
    sizeFactors(dds) <- SFs
  }
  dds <- DESeq(dds)
  res <- lfcShrink(dds, coef=paste0(cond,'_',gp2,'_vs_',gp1), type="apeglm") # for more accurate estimation of log2FC
  #resOri <- results(dds, filterFun=ihw) # for more accurate adjusted p values (based on mean)
  #res$padj <- resOri$padj
  baseMeanConds <- sapply(levels(dds[[cond]]), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds[[cond]] == lvl]))

  # MA plot
  png(paste0(gp2,'_vs_',gp1,'_MAplot.png'), width=600, height=500)
  plotMA(res, ylim=c(-4,4))
  abline(h=-1, lty='dashed',col='yellow', lwd=2)
  abline(h=1, lty='dashed',col='yellow', lwd=2)
  dev.off()

  res <- as.data.frame(res)
  print(head(res))
  res <- cbind(res,baseMeanConds)
  res$geneID <- rownames(res)
  if (is.null(biomaRT_IDs)) {
    biomaRT_IDs.df <- biomaRT_convertID(res$geneID, species)
  } else {
    biomaRT_IDs.df <- biomaRT_IDs
  }
  res.final <- merge(res, biomaRT_IDs.df, by.x = 'geneID', by.y = 'ensembl_gene_id', all.x = T)
  #res.final <- merge(res, biomaRT_IDs.df, by.x = 'geneID', by.y = 'hgnc_symbol', all.x = T)
  res.final <- res.final[,c(1,9,10,11,12,13,14,2,7,8,3,4,5,6)]
  colnames(res.final) <- c('geneID','geneName','geneDescription','chrom','start','end', 'biotype', 'baseMean', gp1, gp2,
                           'log2FoldChange','lfcSE','pvalue','padj')
  res.final$padj <- p.adjust(res.final$pvalue, method = 'BH')
  res.final <- res.final[order(res.final$pvalue),]
  print(head(res.final))
  # # Volcano plot
  # pdf(file = paste0(gp1,'_vs_',gp2,'_volcano.pdf'), width = 6, height = 6)
  # with(subset(res.final, padj < 0.05), plot(log2FoldChange, -log10(pvalue), pch=20, col ='red', main=paste0(gp1,'_vs_',gp2), xlim=c(-10, 10), ylim = c(0, max(-log10(pvalue)))))
  # with(subset(res.final, padj < 0.05 & log2FoldChange < 0), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  # with(subset(res.final, padj >= 0.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))
  # abline(v=-1, lty='dashed',col='yellow', lwd=2)
  # abline(v=1, lty='dashed',col='yellow', lwd=2)
  # dev.off()

  if (is.null(biomaRT_IDs)) {
    return(list(res.final, biomaRT_IDs.df))
  } else {
    return(res.final)
  }
}

## A method to do differential expression analysis when there are two conditions
deseq_2way <- function(samples, counts, gpA, gpB1, gpB2, condA, condB, species, biomaRT_IDs = NULL, adjustVars = NULL, sizeFactors = NULL) {
  print(paste(gpA, 'samples: ', gpB2, 'vs', gpB1))
  coldata <- samples[samples[[condA]] == gpA & samples[[condB]] %in% c(gpB1, gpB2), ]
  cts <- counts[, rownames(coldata)]
  coldata[[condA]] <- droplevels(factor(coldata[[condA]]))
  coldata[[condB]] <- droplevels(factor(coldata[[condB]]))

  if (length(levels(coldata[[condB]])) == 0 | length(levels(coldata[[condB]])) == 1) {
    print("Not enough samples for testing. Skipping..")
    return(NULL)
  }

  for (adjustVar in adjustVars) {
    coldata[[adjustVar]] <- droplevels(factor(coldata[[adjustVar]]))
    if (length(levels(coldata[[adjustVar]])) < 2) {
      adjustVars = adjustVars[adjustVars != adjustVar]
    }
  }
  print(coldata)
  if (is.null(adjustVars) | length(adjustVars) == 0) {
    fml <- paste('~', condB)
  } else {
    fml <- paste('~', paste(adjustVars,collapse=" + "), '+', condB)
  }
  print(paste('Model design:', fml))
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = formula(fml))
  dds[[condB]] <- factor(dds[[condB]], levels=c(gpB1, gpB2)) # to check
  
  if (is.null(sizeFactors)) {
    dds <- estimateSizeFactors(dds)
  } else {
    SFs <- sizeFactors[names(sizeFactors) %in% rownames(coldata)]
    print(SFs)
    sizeFactors(dds) <- SFs
  }
  dds <- DESeq(dds)
  res <- lfcShrink(dds, coef=paste0(condB,'_',gpB2,'_vs_',gpB1), type="apeglm") # for more accurate estimation of log2FC
  #resOri <- results(dds, filterFun=ihw) # for more accurate adjusted p values (based on mean)
  #res$padj <- resOri$padj
  baseMeanConds <- sapply(levels(dds[[condB]]), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds[[condB]] == lvl]))
  # MA plot
  png(paste0(gpA,'__',gpB2,'_vs_',gpB1,'_MAplot.png'), width=600, height=500)
  plotMA(res, ylim=c(-4,4))
  abline(h=-1, lty='dashed',col='yellow', lwd=2)
  abline(h=1, lty='dashed',col='yellow', lwd=2)
  dev.off()

  res <- as.data.frame(res)
  res <- cbind(res,baseMeanConds)
  res$geneID <- rownames(res)
  if (is.null(biomaRT_IDs)) {
    biomaRT_IDs.df <- biomaRT_convertID(res$geneID, species)
  } else {
    biomaRT_IDs.df <- biomaRT_IDs
  }
  res.final <- merge(res, biomaRT_IDs.df, by.x = 'geneID', by.y = 'ensembl_gene_id', all.x = T)
  #res.final <- merge(res, biomaRT_IDs.df, by.x = 'geneID', by.y = 'hgnc_symbol', all.x = T)
  res.final <- res.final[,c(1,9,10,11,12,13,14,2,7,8,3,4,5,6)]
  colnames(res.final) <- c('geneID','geneName','geneDescription','chrom','start','end','biotype', 'baseMean', gpB1, gpB2,
                              'log2FoldChange','lfcSE','pvalue','padj')
  res.final$padj <- p.adjust(res.final$pvalue, method = 'BH')
  res.final <- res.final[order(res.final$pvalue),]
  print(head(res.final))
  # # Volcano plot
  # pdf(file = paste0(gpA,'_',gpB1,'_vs_',gpB2,'_volcano.pdf'), width = 7, height = 7)
  # with(subset(res.final, padj < 0.05), plot(log2FoldChange, -log10(pvalue), pch=20, col ='red', main=paste0(gpA,': ',gpB1,'_vs_',gpB2), xlim=c(-10, 10), ylim = c(0, max(-log10(pvalue)))))
  # with(subset(res.final, padj < 0.05 & log2FoldChange < 0), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  # with(subset(res.final, padj >= 0.05), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))
  # abline(v=-1, lty='dashed',col='yellow', lwd=2)
  # abline(v=1, lty='dashed',col='yellow', lwd=2)
  # dev.off()

  if (is.null(biomaRT_IDs)) {
    return(list(res.final, biomaRT_IDs.df))
  } else {
    return(res.final)
  }
}

# A method to draw PCA when there is only one condition
drawPCA_1way <- function(pcaData, percentVar, cond) {
  p <- ggscatter(pcaData, x = 'PC1', y = 'PC2', color = cond,
                 size = 4, label = 'name', repel = TRUE, font.label = c(10, "plain"),
                 xlab = paste0("PC1: ", percentVar[1], "% variance"),
                 ylab = paste0("PC2: ", percentVar[2], "% variance"),
                 palette = 'npg') +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14))
  return(p)
}

# A method to draw PCA when there are two conditions
drawPCA_2way <- function(pcaData, percentVar, cond1, cond2) {
  p <- ggscatter(pcaData, x = 'PC1', y = 'PC2', color = cond1, shape = cond2,
                 size = 4, label = 'name', repel = TRUE, font.label = c(10, "plain"),
                 xlab = paste0("PC1: ", percentVar[1], "% variance"),
                 ylab = paste0("PC2: ", percentVar[2], "% variance"),
                 palette = 'npg') +
    theme_linedraw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14))
  return(p)
}


### Main
## 1. Read counts and sample data
args <- commandArgs(TRUE)
samples <- read.csv(args[1], check.names=F, header=T, row.names = 1)
counts <- read.csv(args[2], check.names=F, header=T, row.names = 1)
testVars <- unlist(strsplit(args[3], ' '))
adjustVars <- unlist(strsplit(args[4], ' '))
if (length(adjustVars) == 1) {
  if (adjustVars == "none") {
    adjustVars = NULL
  }
}
species <- args[5]
cores <- args[6]
levels1 <- unlist(strsplit(args[7], ' ')) # ordered levels for testVars or testVars1
levels2 <- unlist(strsplit(args[8], ' ')) # ordered levels for testVars2

if (length(args) == 7 | length(args) == 8) {
  sizeFactors = NULL
} else {
  print('Custom size factors')
  sizeFactorCounts <- read.csv(args[9], check.names=F, header=T, row.names = 1) # count data used to estimate pre-defined size factors
  sizeFactorCounts <- sizeFactorCounts[,colnames(sizeFactorCounts) %in% rownames(samples)] # if not every sample in the counts is in the sample list
  print(head(sizeFactorCounts))
  if (length(testVars) == 1) {
    dds <- DESeqDataSetFromMatrix(countData = sizeFactorCounts, colData = samples, design = formula(paste('~', testVars)))
  } else if (length(testVars) == 2) {
    dds <- DESeqDataSetFromMatrix(countData = sizeFactorCounts, colData = samples, design = formula(paste('~', testVars[1], '+', testVars[2])))
  }
  print(paste('Obtain size factors using the count file', args[9]))
  dds <- estimateSizeFactors(dds)
  sizeFactors <- sizeFactors(dds)
  print(sizeFactors)
}

print(paste('Running DESeq2 using', cores, 'cores'))
register(MulticoreParam(cores))

## 2. Differential gene expression analysis
# Only one condition (one-way)
if (length(testVars) == 1) {
  cond <- testVars
  #cond_gps <- levels(factor(samples[[cond]]))
  cond_gps <- levels1
  numSigGenes <- data.frame(comb=NA, total=NA, AgtB=NA, BgtA=NA)

  biomaRT_IDs <- NULL
  gp_res <- list()

  cond_gps_combn <- combn(cond_gps, 2, simplify = F)
  for (i in 1:length(cond_gps_combn)) {
    res_name <- str_trunc(paste0(cond_gps_combn[[i]][2], '_vs_', cond_gps_combn[[i]][1]), 31, ellipsis='')
    if (is.null(biomaRT_IDs)) {
      temp_list <- deseq_1way(samples, counts, cond_gps_combn[[i]][1], cond_gps_combn[[i]][2], cond, species, adjustVars = adjustVars, sizeFactors = sizeFactors)
      gp_res[[res_name]] <- temp_list[[1]]
      biomaRT_IDs <- temp_list[[2]]
    } else {
      gp_res[[res_name]] <- deseq_1way(samples, counts, cond_gps_combn[[i]][1], cond_gps_combn[[i]][2], cond, species, biomaRT_IDs, adjustVars = adjustVars, sizeFactors = sizeFactors)
    }
    n1 <- nrow(subset(gp_res[[res_name]], padj < 0.05))
    n2 <- nrow(subset(gp_res[[res_name]], log2FoldChange > 0 & padj < 0.05))
    n3 <- nrow(subset(gp_res[[res_name]], log2FoldChange < 0 & padj < 0.05))
    numSigGenes <- rbind(numSigGenes, data.frame(comb=res_name, total=n1, AgtB=n2, BgtA=n3))
  }
  write.xlsx(gp_res, paste0(cond,'_comparison_DESeq2_res.xlsx'), colNames = T, keepNA = T, na.string = 'NA', colWidths = 'auto', overwrite = T)

  numSigGenes <- numSigGenes[-c(1), ]
  write.table(numSigGenes, file = 'numSigGenes.txt', quote = F, row.names = F, sep = '\t')
}

# Two conditions (two-way)
if (length(testVars) == 2) {
  condA <- testVars[1]
  condB <- testVars[2]
  #condA_gps <- levels(factor(samples[[condA]]))
  #condB_gps <- levels(factor(samples[[condB]]))
  condA_gps <- levels1
  condB_gps <- levels2
  numSigGenes <- data.frame(comb=NA, total=NA, AgtB=NA, BgtA=NA)

  biomaRT_IDs <- NULL
  gpA_res <- list()
  for (gpA in condA_gps) {
    condB_gps_combn <- combn(condB_gps, 2, simplify = F)
    for (i in 1:length(condB_gps_combn)) {
      res_name <- str_trunc(paste0(gpA, '__', condB_gps_combn[[i]][2], '_vs_', condB_gps_combn[[i]][1]), 31, ellipsis='') # excel only allows a max of 31 characters as tab names
      if (is.null(biomaRT_IDs)) {
        temp_list <- deseq_2way(samples, counts, gpA, condB_gps_combn[[i]][1], condB_gps_combn[[i]][2], condA, condB, species, adjustVars = adjustVars, sizeFactors = sizeFactors)
        if (! is.null(temp_list)) {
          gpA_res[[res_name]] <- temp_list[[1]]
          biomaRT_IDs <- temp_list[[2]]
        } else {
          next
        }
      } else {
        temp_list <- deseq_2way(samples, counts, gpA, condB_gps_combn[[i]][1], condB_gps_combn[[i]][2], condA, condB, species, biomaRT_IDs, adjustVars = adjustVars, sizeFactors = sizeFactors)
        if (! is.null(temp_list)) {
          gpA_res[[res_name]] <- temp_list
        } else {
          next
        }
      }
      n1 <- nrow(subset(gpA_res[[res_name]], padj < 0.05))
      n2 <- nrow(subset(gpA_res[[res_name]], log2FoldChange > 0 & padj < 0.05))
      n3 <- nrow(subset(gpA_res[[res_name]], log2FoldChange < 0 & padj < 0.05))
      numSigGenes <- rbind(numSigGenes, data.frame(comb=res_name, total=n1, AgtB=n2, BgtA=n3))
    }
  }
  write.xlsx(gpA_res, paste0(condB,'_comparison_DESeq2_res.xlsx'), colNames = T, keepNA = T, na.string = 'NA', colWidths = 'auto', overwrite = T)

  gpB_res <- list()
  for (gpB in condB_gps) {
    condA_gps_combn <- combn(condA_gps, 2, simplify = F)
    for (i in 1:length(condA_gps_combn)) {
      res_name <- str_trunc(paste0(gpB, '__', condA_gps_combn[[i]][2], '_vs_', condA_gps_combn[[i]][1]), 31, ellipsis='') # excel only allows a max of 31 characters as tab names
      if (is.null(biomaRT_IDs)) {
        temp_list <- deseq_2way(samples, counts, gpB, condA_gps_combn[[i]][1], condA_gps_combn[[i]][2], condB, condA, species, adjustVars = adjustVars, sizeFactors = sizeFactors)
        if (! is.null(temp_list)) {
          gpB_res[[res_name]] <- temp_list[[1]]
          biomaRT_IDs <- temp_list[[2]]
        } else {
          next
        }
      } else {
        temp_list <- deseq_2way(samples, counts, gpB, condA_gps_combn[[i]][1], condA_gps_combn[[i]][2], condB, condA, species, biomaRT_IDs, adjustVars = adjustVars, sizeFactors = sizeFactors)
        if (! is.null(temp_list)) {
          gpB_res[[res_name]] <- temp_list
        } else {
          next
        }
      }
      n1 <- nrow(subset(gpB_res[[res_name]], padj < 0.05))
      n2 <- nrow(subset(gpB_res[[res_name]], log2FoldChange > 0 & padj < 0.05))
      n3 <- nrow(subset(gpB_res[[res_name]], log2FoldChange < 0 & padj < 0.05))
      numSigGenes <- rbind(numSigGenes, data.frame(comb=res_name, total=n1, AgtB=n2, BgtA=n3))
    }
  }
  write.xlsx(gpB_res, paste0(condA,'_comparison_DESeq2_res.xlsx'), colNames = T, keepNA = T, na.string = 'NA', colWidths = 'auto', overwrite = T)

  numSigGenes <- numSigGenes[-c(1), ]
  write.table(numSigGenes, file = 'numSigGenes.txt', quote = F, row.names = F, sep = '\t')
}


## 3. PCA analysis
print("Obtain DESeq2 run summary files..")
samples <- samples[rownames(samples) %in% colnames(counts), ,drop=F] # if not every sample in the condition is in the count data
counts <- counts[, rownames(samples)]
if (length(testVars) == 1) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = formula(paste('~', testVars)))
} else if (length(testVars) == 2) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples, design = formula(paste('~', testVars[1], '+', testVars[2])))
}
if (length(args) == 7 | length(args) == 8) {
  dds <- estimateSizeFactors(dds)
} else {
  sizeFactors(dds) <- sizeFactors
}
## 4. Save R objects
saveRDS(dds, 'dds.rds')
write.csv(counts(dds, normalized = F), 'raw_counts.csv', quote = F)
write.csv(counts(dds, normalized = T), 'normalized_counts.csv', quote = F)
write.table(sizeFactors(dds), 'size_factors.csv', sep = ',', quote = F, col.names = F)

dds <- DESeq(dds)


if (nrow(dds) >= 1000){
  vsd <- vst(dds, blind=T)
} else {
  vsd <- varianceStabilizingTransformation(dds)
}
rv <- rowVars(assay(vsd)) # calculate variance (or whatever measure you prefer) per gene
rv_top <- order(rv, decreasing = TRUE)[1:500] # sort, so that the most variable genes will be on top of the object
allgenes <- rownames(dds@assays@data$counts)
topVarGenes <- allgenes[rv_top]
topVarGenes_biomaRT.df <- biomaRT_convertID(topVarGenes, species)

pcaData <- plotPCA(vsd, intgroup=testVars, returnData = T)
percentVar <- round(100 * attr(pcaData, "percentVar"))


png('samples_PCA.png', width=600, height=500)

if (length(testVars) == 1) {
  cond <- testVars
  p <- drawPCA_1way(pcaData, percentVar, cond)
} else if (length(testVars) == 2) {
  condA <- testVars[1]
  condB <- testVars[2]
  condA_gps <- levels(samples[[condA]])
  condB_gps <- levels(samples[[condB]])
  if (length(condA_gps) >= length(condB_gps)) {
    p <- drawPCA_2way(pcaData, percentVar, condA, condB)
  } else {
    p <- drawPCA_2way(pcaData, percentVar, condB, condA)
  }
}
print(p)
dev.off()

## 4. Save R objects (cont.)
write.xlsx(topVarGenes_biomaRT.df, 'topVarGenes.xlsx', colNames = T, keepNA = T, na.string = 'NA', colWidths = 'auto', overwrite = T)

