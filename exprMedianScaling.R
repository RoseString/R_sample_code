setwd("DIR/TO/EXPRESSION_DATA")

chickenCondDir <- "chicken_cond/"
titCondDir <- "blue_tit_cond/"
ostrichCondDir <- "ostrich_cond/"

TPM_scale <- function(tissue){
  df <- read.table("chicken_tit_ostrich.TPM", header = T) # gene expression (measured by TPM) matrix (row -> gene; column -> sample)
  rownames(df) <- df$Gene
  df$Gene <- NULL
  expr.df <- log2(df+1)
  
  ### Extract corresponding samples
  chickenCond <- read.table(paste(chickenCondDir, tissue,"_adult.cond", sep = ''), col.names = c("Sample", "Sex", "Type")) # each "XXX.cond" contains samples of a particular tissue/developmental stage;
  titCond <- read.table(paste(titCondDir, tissue,"_adult.cond", sep = ''), col.names = c("Sample", "Sex", "Type"))   # (row -> sample; 1st column -> sample name; other columns are not necessary in this program)
  ostrichCond <- read.table(paste(ostrichCondDir, tissue,"_adult.cond", sep = ''), col.names = c("Sample", "Sex", "Type"))
  
  samples <- c(as.character(chickenCond$Sample), as.character(titCond$Sample), as.character(ostrichCond$Sample))
  expr.df <- expr.df[,samples, drop=F] # raw expression data for selected samples
  
  boxplot(expr.df, main=paste(tissue, "\nBefore Scaling"))
  
  ### For each sample, calculate IQR of TPM values
  lowers <- list()
  uppers <- list()
  for (sample in samples) {
    expr.expressed.df <- expr.df[expr.df[,sample] >= 0,] # filter out unexpressed genes for this step
    lowers[[sample]] = quantile(expr.expressed.df[,sample])[[2]]
    uppers[[sample]] = quantile(expr.expressed.df[,sample])[[4]]
  }
  
  ### For each sample, get genes that have TPM values within the IQR
  iqrGenes <- list()
  for (sample in samples) {
    expr <- expr.df[,sample]
    goodRows <- which(expr > lowers[[sample]] & expr < uppers[[sample]])
    iqrGenes[[sample]] <- rownames(expr.df[goodRows,])
  }
  
  ### Across samples, rank genes by their frequencies within IQRs
  iqrTimes <- list()
  genes <- row.names(expr.df)
  for (gene in genes) {
    iqrTimes[[gene]] = 0
  }
  
  for (sample in samples) {
    print(sample)
    for (iqrGene in iqrGenes[[sample]]) {
      iqrTimes[[iqrGene]] = iqrTimes[[iqrGene]] + 1
    }
  }
  
  iqrTimes.df <- as.data.frame(unlist(iqrTimes))
  colnames(iqrTimes.df) <- "Times"
  iqrTimes.ranked.df <- iqrTimes.df[order(iqrTimes.df$Times, decreasing = T), , drop = F]
  
  
  ### Only for genes with the conserved ranks (present in all samples), get the median value of expression in each sample
  ## The code can be altered to select top conserved genes 
  consGenes <- rownames(subset(iqrTimes.ranked.df, Times == length(expr.df)))
  print(paste(length(consGenes), "genes used for median normalization in", tissue))
  consMedian <- list()
  for (sample in samples) {
    expr <- expr.df[,sample, drop=F]
    consMedian[[sample]] <- median(expr[row.names(expr) %in% consGenes,])
  }
  
  ### Calculate scaling factors for each sample (will scale median to mean of medians for genes with the most conserved ranks)
  scalingFactors <- list()
  for (sample in samples) {
    scalingFactors[[sample]] <- consMedian[[sample]] / mean(unlist(consMedian))
  }
  
  ### For each sample, recalculate log2(TPM) by dividing expression by the scaling factors
  expr.scaled.df <- expr.df
  for (sample in samples) {
    expr.scaled.df[, sample] <- round(expr.scaled.df[, sample] / scalingFactors[[sample]], 6)
  }
  
  ### Output scaled results
  expr.scaled.df <- data.frame(Gene = rownames(expr.scaled.df), expr.scaled.df)
  write.table(expr.scaled.df, paste("chicken_tit_ostrich_", tissue, ".scaledLog2.TPM", sep = ""), row.names = F, quote = F, sep = "\t") # output scaled gene expression data
  
  boxplot(expr.scaled.df[,-1], main=paste(tissue, "\nAfter Scaling"))
}

par(mfrow=c(2,1)) # To visualize gene expression across samples before and after scaling
TPM_scale("brain")
