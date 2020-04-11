# feature selection
# scRNA-seq data is often used in exploratory analyses to characterize heterogenity
# across cells.
# clustering and dim reduction compares cels based on their gene expression profiles
# which involves aggregating per gene differences into a single (dis)similarity metric
# between a pair of cells. 
# => we want to select genes that contain useful information about the biology of 
# the system while removeing the genes that contain random noise.
# => preserve interesting biological structure without the variance that obscures 
# that structure
# also reduces size of the data set and thus improves computational efficiency in 
# later steps
# => simplest approach is to select the most variable genes based on their 
# expression across the population. This assumes that genuine biological differences
# will manifest as increased variation in the affected genes, compared to other 
# genes that are only affected by technical noise or a baseline level of uninteresting
# biological variation. 
# => Highly variable genes (HVGs)

# data import: pbmc dataset 
#--- loading ---#
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

#--- gene-annotation ---#
library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)

# BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")

#--- cell-detection ---#
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]

#--- quality-control ---#
stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]

#--- normalization ---#
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

# data import: 416b dataset
#--- loading ---#
library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b") 
sce.416b$block <- factor(sce.416b$block)

#--- gene-annotation ---#
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
rowData(sce.416b)$ENSEMBL <- rownames(sce.416b)
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                                   keytype="GENEID", column="SYMBOL")
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                                    keytype="GENEID", column="SEQNAME")

library(scater)
rownames(sce.416b) <- uniquifyFeatureNames(rowData(sce.416b)$ENSEMBL, 
                                           rowData(sce.416b)$SYMBOL)

#--- quality-control ---#
mito <- which(rowData(sce.416b)$SEQNAME=="MT")
stats <- perCellQCMetrics(sce.416b, subsets=list(Mt=mito))
qc <- quickPerCellQC(stats, percent_subsets=c("subsets_Mt_percent",
                                              "altexps_ERCC_percent"), batch=sce.416b$block)
sce.416b <- sce.416b[,!qc$discard]

#--- normalization ---#
library(scran)
sce.416b <- computeSumFactors(sce.416b)
sce.416b <- logNormCounts(sce.416b)

# variance of the log-counts 
# simplest approach to quantifying per gene variation
# compute variance of log normalized expression values
# this has the advantage of basing the feature selection on the same log values 
# that are used for later steps. In particular, genes with the largest variances
# in log values will contribute the most to the euclidian distances between cells. 

# simply calculating the per-gene variance is not enough, feature selection 
# requires modeling of the mean-variance relationship. 
# The log transformation does not achieve perfecet variance stabilization, which 
# means that the variance of a gene is driven more by its abundance than its 
# underlying biological heterogenity. To account fo this effect we use modelGeneVar()
# to fit a trend to the variance with respect to abundance across all genes.

library(scran)
dec.pbmc <- modelGeneVar(sce.pbmc)

# Visualizing the fit:
fit.pbmc <- metadata(dec.pbmc)
plot(fit.pbmc$mean, fit.pbmc$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# we assume that the expression proviles of most genes are dominated by random 
# technical noise. Under this assumption our trend represents an estimate of the 
# technical noise as a function of abundance. 
# We then break down the total variance of each gene into the technical component
# => the fitted value of the trend at that genes abundance 
# and the biological component 
# => difference between total variance and technical component. 
# => this biological component represents the interesting variation for each gene 
# and is used as the metric for high variance gene selection. 

# Ordering by most interesting genes for inspection.
dec.pbmc[order(dec.pbmc$bio, decreasing=TRUE),] 

# alternative approach to quantification uses squared coeff of variation 
# of normalized expression values prior to log transform

dec.cv2.pbmc <- modelGeneCV2(sce.pbmc)

# this allows us the model the mean variance relationship when considering the 
# relevance of each gene
# => core assumption is that most genes contain random noise and that the trend 
# captures mostly technical variation. 
# => large cv^2 values that deviate from the trend are likely to represent genes 
# affected by biological structure. 

fit.cv2.pbmc <- metadata(dec.cv2.pbmc)
plot(fit.cv2.pbmc$mean, fit.cv2.pbmc$cv2, log="xy")
curve(fit.cv2.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# for each gene we quantify the deviation from the trend in terms of the ratio 
# of its cv^2 to the fitted value of trend at its abundance 
# this is more appropriate than directly subtracting the trend from the cv^2, as 
# the magnitude of the ratio is not affected by the mean. 

dec.cv2.pbmc[order(dec.cv2.pbmc$ratio, decreasing=TRUE),]

# both cv^2 and the variance of log counts are effective metrics for quantifying 
# variation in gene expression. 
# cv^2 tends to give higher rank to low abundance hvgs driven by upregulation in 
# rare subpopulations for which increase in variance on the raw scale is stronger 
# than that on the log scale
# however cv2 variance is less directly related to downstream procedures operating 
# on the log counts and the reliance on the ratio can assign high rank to 
# uninteresting genes with low absolute variance. 

# quantifying technical noise: 
# strictly speaking the noise we are not interested in is a compound of technical 
# noise and "uninteresting" biological noise caused by events like 
# transcriptional bursting
# this revised assumption is generally reasonable, but may be problematic in some 
# scenarios where many genes at a particular abundance are affected by a 
# biological process. 
# => strong upregulation may result in an enrichment of high variance genes at
# high abundances 
# => skews the fitted trend in abundance interval and compromises detection of 
# relevant genes. 
# avoid this problem by fitting a mean dependent trend to the variance of the 
# spike in transcripts if they are available 
# spike ins should not be affected by biological variation 

dec.spike.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC")
dec.spike.416b[order(dec.spike.416b$bio, decreasing=TRUE),]

plot(dec.spike.416b$mean, dec.spike.416b$total, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
fit.spike.416b <- metadata(dec.spike.416b)
points(fit.spike.416b$mean, fit.spike.416b$var, col="red", pch=16)
curve(fit.spike.416b$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# if spike in data is not available one can make a few distributional assumptions about
# the data
# for example, unique molecular identifiers usually exhibit poisson variation

set.seed(0010101)
dec.pois.pbmc <- modelGeneVarByPoisson(sce.pbmc)
dec.pois.pbmc <- dec.pois.pbmc[order(dec.pois.pbmc$bio, decreasing=TRUE),]
head(dec.pois.pbmc)

plot(dec.pois.pbmc$mean, dec.pois.pbmc$total, pch=16, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(metadata(dec.pois.pbmc)$trend(x), col="dodgerblue", add=TRUE)

# trends based purely on technical noise tend to yield large biological components 
# for highly expressed genes. 
# this often includes "house-keeping" genes for essential cellular components 
# such as ribosomal proteins which are usually considered to be uninteresting for 
# characterizing cellular heterogenity. 
# These observations suggest that a more accurate noise model does not necessarily
# yield a better ranking of high variance genes. 

### controlling for batch effects
# use modelgenevarwithspikes function to treat every plate (block) in the 416b 
# data as a different batch 
dec.block.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC", block=sce.416b$block)
head(dec.block.416b[order(dec.block.416b$bio, decreasing=TRUE),1:6])

# take advantage of the information provided in every batch
par(mfrow=c(1,2))
blocked.stats <- dec.block.416b$per.block
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(current$mean, current$total, main=i, pch=16, cex=0.5,
       xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit <- metadata(current)
  points(curfit$mean, curfit$var, col="red", pch=16)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2) 
}

# one may also use a design matrix to block specific genes like so
design <- model.matrix(~factor(block) + phenotype, colData(sce.416b))
dec.design.416b <- modelGeneVarWithSpikes(sce.416b, "ERCC", design=design)
dec.design.416b[order(dec.design.416b$bio, decreasing=TRUE),]

# in the end we always select by a given criterium 
# and proceed with care