library(data.table)
library(IRanges)
library(GenomicRanges)
library(AnnotationHub)
library(rtracklayer)
library(ggpubr)
library(reshape2)
setwd("~/mountpoint/XCI/1.Exploratory_Analysis/1.XaXi_meth/")

### Method 1: Normalization based on autosomes
## Calculate Xi methylation levels
M.df <- fread("1527_NeuN.X.cov")
F.df <- fread("AN15240_NeuN.X.cov")

df <- merge(M.df, F.df, by = "V2") # In total, 1,092,960 CpGs overlapping between the male and female samples

Xi.M <- df$V5.y - df$V5.x
Xi.U <- df$V6.y - df$V6.x
df$Xi.M <- Xi.M
df$Xi.U <- Xi.U
df$Xa.M <- df$V5.x
df$Xa.U <- df$V6.x
df <- subset(df, Xi.M > 0 & Xi.U > 0 & Xi.M+Xi.U >= 5) # kept 740,038 CpGs (passed Xi > 5 reads cutoff; Xi M and U reads not negative after calculation)
XaXi.df <- data.frame(chrom="chrX", start=df$V2, end=df$V3.x, 
                      Xa=df$Xa.M*100/(df$Xa.M+df$Xa.U), Xa_M=df$Xa.M, Xa_U=df$Xa.U,  
                      Xi=df$Xi.M*100/(df$Xi.M+df$Xi.U), Xi_M=df$Xi.M, Xi_U=df$Xi.U)

## Convert Xa and Xi information along with Coordinates into a GRanges object
meth.df <- makeGRangesFromDataFrame(XaXi.df, keep.extra.columns = T)

## Overlap human hg19 promoter and genebody ranges with DNA methylation data
aHub <- AnnotationHub()
aHub <- subset(aHub, species == "Homo sapiens")
genes <- query(aHub, "RefSeq")[[1]]
Xpromoter <- keepSeqlevels(promoters(genes, upstream = 1000, downstream = 500, use.names = T), "chrX", pruning.mode="coarse")
Xgenebody <- keepSeqlevels(genes, "chrX", pruning.mode="coarse")

meth.Xpromoter <- subsetByOverlaps(meth.df, Xpromoter)
meth.Xgenebody <- subsetByOverlaps(meth.df, Xgenebody)

## Plot promoter and genebody DNA methylation data
# Promoter
meth.Xpromoter.df <- melt(subset(data.frame(meth.Xpromoter), select=c("seqnames","start","end","Xa","Xi")), 
                          id.vars = c("seqnames","start","end"), 
                          variable.name = "chrom", value.name = "ml")
p1 <- ggboxplot(meth.Xpromoter.df, "chrom", "ml", fill = "chrom", palette = "jco", main = "Promoter", width = 0.5) +
  rremove("legend") + rremove("xlab") + rremove("ylab") + rremove("x.ticks") +
  theme(legend.spacing.x = unit(0.4, 'cm'), 
        legend.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=13),
        plot.title = element_text(size=15, hjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "lines"))

# Genebody
meth.Xgenebody.df <- melt(subset(data.frame(meth.Xgenebody), select=c("seqnames","start","end","Xa","Xi")), 
                          id.vars = c("seqnames","start","end"), 
                          variable.name = "chrom", value.name = "ml")
p2 <- ggboxplot(meth.Xgenebody.df, "chrom", "ml", fill = "chrom", palette = "jco", main = "genebody", width = 0.5) +
  rremove("legend") + rremove("xlab") + rremove("ylab") + rremove("x.ticks") +
  theme(legend.spacing.x = unit(0.4, 'cm'), 
        legend.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=13),
        plot.title = element_text(size=15, hjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "lines"))

figure <- ggarrange(p1, p2, ncol = 2, nrow = 1)
annotate_figure(figure, left = text_grob("5mC (%)", rot = 90, size = 14))


## Statistical test
wilcox.test(meth.Xpromoter$Xa, meth.Xpromoter$Xi, paired = T) #  p-value < 2.2e-16
wilcox.test(meth.Xgenebody$Xa, meth.Xgenebody$Xi, paired = T) #  p-value < 2.2e-16








### Method 2: Calculation (2*Female_Meth - Male_Meth)
## Calculate Xi methylation levels
M.df <- fread("1527_NeuN.X.cov")
F.df <- fread("AN15240_NeuN.X.cov")

df <- merge(M.df, F.df, by = "V2") # In total, 1,092,960 CpGs overlapping between the male and female samples

XaXi.df <- data.frame(chrom="chrX", start=df$V2, end=df$V3.x, 
                      Xa=df$V4.x,
                      Xi=2*df$V4.y-df$V4.x)
XaXi.df[XaXi.df$Xi>100,]$Xi <- 100
XaXi.df[XaXi.df$Xi<0,]$Xi <- 0

## Convert Xa and Xi information along with Coordinates into a GRanges object
meth.df <- makeGRangesFromDataFrame(XaXi.df, keep.extra.columns = T)

## Overlap human hg19 promoter and genebody ranges with DNA methylation data
aHub <- AnnotationHub()
aHub <- subset(aHub, species == "Homo sapiens")
genes <- query(aHub, "RefSeq")[[1]]
Xpromoter <- keepSeqlevels(promoters(genes, upstream = 1000, downstream = 500, use.names = T), "chrX", pruning.mode="coarse")
Xgenebody <- keepSeqlevels(genes, "chrX", pruning.mode="coarse")

meth.Xpromoter <- subsetByOverlaps(meth.df, Xpromoter)
meth.Xgenebody <- subsetByOverlaps(meth.df, Xgenebody)

## Plot promoter and genebody DNA methylation data
# Promoter
meth.Xpromoter.df <- melt(subset(data.frame(meth.Xpromoter), select=c("seqnames","start","end","Xa","Xi")), 
                          id.vars = c("seqnames","start","end"), 
                          variable.name = "chrom", value.name = "ml")
p1 <- ggboxplot(meth.Xpromoter.df, "chrom", "ml", fill = "chrom", palette = "jco", main = "Promoter", width = 0.5) +
  rremove("legend") + rremove("xlab") + rremove("ylab") + rremove("x.ticks") +
  theme(legend.spacing.x = unit(0.4, 'cm'), 
        legend.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=13),
        plot.title = element_text(size=15, hjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "lines"))

# Genebody
meth.Xgenebody.df <- melt(subset(data.frame(meth.Xgenebody), select=c("seqnames","start","end","Xa","Xi")), 
                          id.vars = c("seqnames","start","end"), 
                          variable.name = "chrom", value.name = "ml")
p2 <- ggboxplot(meth.Xgenebody.df, "chrom", "ml", fill = "chrom", palette = "jco", main = "genebody", width = 0.5) +
  rremove("legend") + rremove("xlab") + rremove("ylab") + rremove("x.ticks") +
  theme(legend.spacing.x = unit(0.4, 'cm'), 
        legend.text=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=13),
        plot.title = element_text(size=15, hjust = 0.5),
        plot.margin = unit(c(1,1,1,1), "lines"))

figure <- ggarrange(p1, p2, ncol = 2, nrow = 1)
annotate_figure(figure, left = text_grob("5mC (%)", rot = 90, size = 14))

## Statistical test
wilcox.test(meth.Xpromoter$Xa, meth.Xpromoter$Xi, paired = T) #  p-value < 2.2e-16
wilcox.test(meth.Xgenebody$Xa, meth.Xgenebody$Xi, paired = T) #  p-value < 2.2e-16


