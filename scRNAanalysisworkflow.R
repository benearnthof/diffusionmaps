#scRNA-seq analysis workflow
# Before starting the analysis itself, some comments on experimental design may 
# be helpful. The most obvious question is the choice of technology, which can 
# be roughly divided into:
#   
# Droplet-based: 10X Genomics, inDrop, Drop-seq
# Plate-based with unique molecular identifiers (UMIs): CEL-seq, MARS-seq
# Plate-based with reads: Smart-seq2
# Other: sci-RNA-seq, Seq-Well

# sequencing data from scRNA-seq experiments must be converted into a matrix of 
# expression values that can be used for statistical analysis
# discrete data => Unique molecular identifiers or reads mapped to each gene in 
# each cell. 

# after quantification the count matrix is imported into R and converted to a 
# SingleCellExperiment object. 

# Data processing and downstream analysis: 
# 1. Compute quality control metrics to remove low quality cells
# => Total counts per cell, proportion of spike-in or mitochondrial reads and 
# number of detected features.
# 2. Convert counts into normalized expression values to eliminate cell specific
# biases like in capture efficiency. This allows for explicit comparisons across 
# cells in downstream steps like clustering.
# Also apply a transformation (typically log transformation) to adjust for the 
# mean-variance relationship. 
# 3. Perform feature selection to pick a subset of interesting features for 
# downstream analysis. 
# => Model Variance across cells for each gene and retain genes that are highly 
# variable. The aim is to reduce computational overhead and noise from 
# uninteresting genes. 
# 4. Apply dimensionality reduction to compact the data and further reduce noise.
# PCA is usually used initially followed by more aggressive methods like tSNE 
# for visualization. 
# 5. Cluster cells into groups according to similarities in their normalized 
# expression profiles. This aims to obtain groupings that serve as empirical 
# proxies for distinct biological states. 

# Example: Dataset from Macosko et al 2015
BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
BiocManager::install("scRNAseq")
library(scRNAseq)

# this will download a bunch of data, create a new directory if possible
sce <- MacoskoRetinaData()

# Quality control.
BiocManager::install("scater")
library(scater)
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
sce <- sce[, !filtered$discard]

# Normalization.
sce <- logNormCounts(sce)

# Feature selection.
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

# Dimensionality reduction.
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)
sce <- runUMAP(sce, dimred = 'PCA', external_neighbors=TRUE)

plotUMAP(sce, colour_by = "cluster")

# Clustering.
g <- buildSNNGraph(sce, use.dimred = 'PCA')
colData(sce)$clust <- factor(igraph::cluster_louvain(g)$membership)
plotUMAP(sce, colour_by = "clust")
plotPCA(sce, colour_by = "clust")
# colnames(sce) <- factor(igraph::cluster_louvain(g)$membership)
# does not work as of yet

# Visualization.
plotUMAP(sce, colour_by="cell.id")
plotPCA(sce, colour_by = "label")

# installing a mirror that exports the function colLabels
devtools::install_github("drisso/SingleCellExperiment")
# fails 
# need to update to R version 3.6.3
# BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
library(installr)
installr::updateR()
devtools::install_github("drisso/SingleCellExperiment")
library(SingleCellExperiment)

sce

# lets try the coloring example from plotUMAP
example_sce <- mockSCE()
example_sce <- logNormCounts(example_sce)
example_sce <- runPCA(example_sce)

plotPCA(example_sce)
plotPCA(example_sce, colour_by = "Mutation_Status")
plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment",
        size_by = "Mutation_Status")

plotPCA(example_sce, ncomponents = 4, colour_by = "Treatment",
        shape_by = "Mutation_Status")

example_sce <- runTSNE(example_sce)
plotTSNE(example_sce, run_args=list(perplexity = 10))

## Same for DiffusionMaps:
example_sce <- runDiffusionMap(example_sce)
plotDiffusionMap(example_sce)

## Same for MDS plots:
example_sce <- runMDS(example_sce)
plotMDS(example_sce)
