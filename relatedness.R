setwd("/nv/hp10/dsun33/sparrow/analysis/19.WGBS_sparrow/0.SNP_calling")
library(SNPRelate)
library(SeqArray)

### Convert VCF to the gds format (LONG runtime!!!)
vcf.fn <- "WGS15.no2Z.snp.vcf.gz" # should exclude sex chromosomes and chr2 before clustering
seqVCF2GDS(vcf.fn, "tmp.gds", parallel = 8) # If using 4 cores, takes ~4hrs for a large WTSP VCF with 14 samples (will save to a gds file by merging 3 tmp_tmpXXX)!

### New session
(genofile <- seqOpen("tmp.gds"))

### Filter variants in LD (Takes ~3hrs for the LD filtering) NOT PERFORMED!!!
#set.seed(12345)
#snpset <- snpgdsLDpruning(genofile, autosome.only = FALSE, ld.threshold=0.2, num.thread = 8)

### Kinship calculation using KING
ibd <- snpgdsIBDKING(genofile, autosome.only = F, type = "KING-robust", num.thread = 4)
rownames(ibd$kinship) <- ibd$sample.id
colnames(ibd$kinship) <- ibd$sample.id
saveRDS(ibd$kinship, file = "kinship_King.RDS")

### IBS calculation
ibs <- snpgdsIBS(genofile, num.thread=4, autosome.only=F)
rownames(ibs$ibs) <- ibs$sample.id
colnames(ibs$ibs) <- ibs$sample.id
saveRDS(ibs$ibs, file = "ibs.RDS")

### PCA
snpset.id <- seqGetData(genofile, "variant.id")
pca <- snpgdsPCA(genofile, snp.id=snpset.id, autosome.only=F, num.thread = 4)
saveRDS(pca, file = "pca.RDS")

seqClose(genofile)


### 1. LOCAL computer (PCA analysis)
setwd("~/mountpoint/sparrow/analysis/19.WGBS_sparrow/0.SNP_calling/")
pca <- readRDS("pca.RDS")

pc.percent <- pca$varprop*100
pc.percent <- round(pc.percent, 4)
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  stringsAsFactors = FALSE)

par(mfrow = c(1,4))
barplot(pc.percent, names.arg=1:length(pc.percent), xlab = "Principal Componenet", ylab = "Percent Variance Explained", col = "burlywood2")
plot(tab$EV1, tab$EV2, xlab="PC1 (8.82% Variance Explained)", ylab="PC2 (8.41% Variance Explained)", col = "cyan4")
text(tab$EV1, tab$EV2, labels = tab$sample.id, pos = 4, cex = 0.65, col = "cyan4")

plot(tab$EV1, tab$EV3, xlab="PC1 (8.82% Variance Explained)", ylab="PC3 (7.92% Variance Explained)", col = "cyan4")
text(tab$EV1, tab$EV3, labels = tab$sample.id, pos = 4, cex = 0.65, col = "cyan4")

plot(tab$EV2, tab$EV3, xlab="PC2 (8.41% Variance Explained)", ylab="PC3 (7.92% Variance Explained)", col = "cyan4")
text(tab$EV2, tab$EV3, labels = tab$sample.id, pos = 4, cex = 0.65, col = "cyan4")


### 2. LOCAL computer (Kinship analysis with KING)
library(ggplot2)
library(ggpubr)
library(reshape2)

kinship <- readRDS("kinship_King.RDS")
kinship <- subset(kinship, select = -c(`N2-C10-12-TSF-1`,`N3-C10-01-TSF-4`,`13-02-TSF`))
kinship <- subset(kinship, ! rownames(kinship) %in% c("N2-C10-12-TSF-1","N3-C10-01-TSF-4","13-02-TSF"))

kinship.dat <- melt(kinship)
kinship.dat[kinship.dat$value < 0,]$value <-  0

ggplot(data = kinship.dat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,0.5), space = "Lab", 
                       name="Kinship Coefficient") +
  rremove("xlab") + rremove("ylab") +
  theme(axis.text.x = element_text(angle = 90, size = 10, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.position="top")


### 3. LOCAL computer (Kinship analysis with KING)
ibs <- readRDS("ibs.RDS")

loc <- cmdscale(1 - ibs, k = 3)
x <- loc[, 1]
y <- loc[, 2]
z <- loc[, 3]

par(mfrow=c(1,3))
plot(x, y, xlab="Dimension 1", ylab="Dimension 2", col = "deeppink4")
text(x, y, labels = colnames(ibs), pos = 4, cex = 0.65, col = "deeppink4")

plot(x, z, xlab="Dimension 1", ylab="Dimension 3", main = "Multidimensional Scaling Analysis (IBS)", col = "deeppink4")
text(x, z, labels = colnames(ibs), pos = 4, cex = 0.65, col = "deeppink4")

plot(y, z, xlab="Dimension 2", ylab="Dimension 3", col = "deeppink4")
text(y, z, labels = colnames(ibs), pos = 4, cex = 0.65, col = "deeppink4")

