# normalization
# aims to remove or dampen the systematic differences between libraries that may arise from 
# technical differences in cDNA capture or PCR amplification efficiency across cells

#--- loading ---#
library(scRNAseq)
sce.zeisel <- ZeiselBrainData()

library(scater)
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel, 
                                      id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))

#--- gene-annotation ---#
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
rowData(sce.zeisel) <- mapIds(org.Mm.eg.db, keys=rownames(sce.zeisel),
                              keytype="SYMBOL", column="ENSEMBL")

#--- quality-control ---#
stats <- perCellQCMetrics(sce.zeisel, subsets=list(
  Mt=rowData(sce.zeisel)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", 
                                              "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[,!qc$discard]

# ilbrary size normalization
library(scater)
lib.sf.zeisel <- librarySizeFactors(sce.zeisel)
summary(lib.sf.zeisel)

# 10 fold size differences 
# lets give 'er a log transform
hist(log10(lib.sf.zeisel), xlab="Log10[Size factor]", col='grey80')
# very nice

# deconvolution, idk what that means
library(scran)
set.seed(100)
# preclustering step
clust.zeisel <- quickCluster(sce.zeisel) 
table(clust.zeisel)
deconv.sf.zeisel <- calculateSumFactors(sce.zeisel, cluster=clust.zeisel)
summary(deconv.sf.zeisel)

plot(lib.sf.zeisel, deconv.sf.zeisel, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(sce.zeisel$level1class)))
abline(a=0, b=1, col="red")
