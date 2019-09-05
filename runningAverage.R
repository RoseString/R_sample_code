library(gtools)

setwd("DIR/TO/CHICKEN_EXPRESSION_DATA")
conds <- c("blastoderm", "brain_E18", "brain_adult", "liver_E18", "liver_E19", "liver_adult", "kidney_E18", "kidney_adult",
           "bursa_E18", "heart_E18", "heart_E19", "heart_adult", "spleen_E18", "spleen_E19", "lung_E18","muscle_E18")

expr.df <- read.table("chicken.scaledLog2.wCoord.TPM", header = T) # First three columns -> gene coordinates (Chrom, Start, End); other columns -> samples; rows -> genes 
rownames(expr.df) <- expr.df$Gene
expr.df$Gene <- NULL
zGenes <- rownames(subset(expr.df, Chrom == "NC_006127.4"))

## A function to calculate running average of numbers and adjust the values to the center coordinate of genes within a window
runningAvg <- function(values, winSize) {
  runAvg <- as.numeric(running(values, width=winSize, allow.fewer = F, pad = T))
  if (winSize %% 2 == 0) {
    runAvg <- c(runAvg[(winSize/2):length(runAvg)], rep(NA, winSize/2-1))
  } else {
    runAvg <- c(runAvg[((winSize+1)/2):length(runAvg)], rep(NA, (winSize+1)/2-1))
  }
  return(runAvg)
}

w <- 10 # window size

par(mfrow=c(4,4))
for (cond in conds){
  # Format male and female gene expression data (sorted by Z coordinates)
  cond.df <- read.table(paste("chicken_cond/", cond, ".cond", sep=""), col.names = c("Sample", "Sex", "Type")) # each "XXX.cond" contains samples of a particular tissue/developmental stage; # (row -> sample; 1st column -> sample name; 2nd column -> sex; 3rd column is not necessary in this program)
  samples <- as.character(cond.df$Sample)
  maleSamples <- as.character(cond.df[cond.df$Sex == "Male",]$Sample)
  femaleSamples <- as.character(cond.df[cond.df$Sex == "Female",]$Sample)
    
  maleExpr.df <- expr.df[zGenes, maleSamples]
  femaleExpr.df <- expr.df[zGenes, femaleSamples]
  chickenMaleExprAvg <- apply(maleExpr.df, 1, mean)
  chickenFemaleExprAvg <- apply(femaleExpr.df, 1, mean)
  coord <- subset(expr.df, Chrom == "NC_006127.4")[,"Start", drop=F]
    
  expr.Z.df <- cbind(coord, chickenMaleExprAvg, chickenFemaleExprAvg)
  expr.Z.ordered.df <- subset(expr.Z.df[with(expr.Z.df, order(Start)),], chickenMaleExprAvg > 1 | chickenFemaleExprAvg > 1)
  
  # Calculate running average of sex differences in gene expression
  run_avgs <- runningAvg(expr.Z.ordered.df$chickenMaleExprAvg - expr.Z.ordered.df$chickenFemaleExprAvg, w)
  
  # Visualize the running average of sex differences in gene expression across the Z chromosome
  plot(expr.Z.ordered.df$Start/1000000, run_avgs,
       type = "p", pch = 19, cex = 0.3, col = "aquamarine4", main=cond, cex.lab=1.5, cex.main=1.6, xlab="", ylab = "")
}

