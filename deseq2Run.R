library(DESeq2)

setwd("DIR/TO/COUNT_DATA")
conds <- c("blastoderm", "brain_E18", "brain_adult", "liver_E18", "liver_E19", "liver_adult", "kidney_E18", "kidney_adult",
           "bursa_E18", "heart_E18", "heart_E19", "heart_adult", "spleen_E18", "spleen_E19", "lung_E18","muscle_E18")

deseqRun <- function(cond) {
  cts <- as.matrix(read.csv("count_matrix.csv", row.names="gene_id")) # the csv file contains the raw count matrix
  coldata <- read.csv(paste0("cond/",cond,".cond"), sep="\t", row.names=1, header=F) # each "XXX.cond" contains samples of a particular tissue/developmental stage; # (1st column -> sample name; 2nd column -> sex; 3rd column -> read type)
  colnames(coldata) <- c("sex","type")
  
  rownames(coldata) <- sub("fb", "", rownames(coldata))
  all(rownames(coldata) %in% colnames(cts))
  cts <- cts[, rownames(coldata)]
  all(rownames(coldata) == colnames(cts))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ sex)
  dds <- DESeq(dds, fitType='local')
  
  res <- results(dds)
  baseMeanMF <- sapply(levels(dds$sex), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$condition == lvl]))
  res <- cbind(as.data.frame(res), baseMeanMF)
  res$Gene <- rownames(res)
  res <- res[,c(9,1,8,7,2,3,4,5,6)]
  
  resOrdered <- res[order(res$padj), ]
  write.table(resOrdered, file=paste0("OUTDIR/",cond,".DESeq2"), quote=F, sep="\t", row.names=F)
}

for (cond in conds) {
  deseqRun(cond)
}
