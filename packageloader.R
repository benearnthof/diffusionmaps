# packageloader
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
library(BiocManager)

## try http:// if https:// URLs are not supported

BiocManager::install("scater")
library(scater)

BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

BiocManager::install(c('scater', 'scran', 'uwot'))

counts_matrix <- data.frame(cell_1 = rpois(10, 10), 
                            cell_2 = rpois(10, 10), 
                            cell_3 = rpois(10, 30))
rownames(counts_matrix) <- paste0("gene_", 1:10)
counts_matrix <- as.matrix(counts_matrix) # must be a matrix object!

sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
sce <- scater::logNormCounts(sce)

counts_100 <- counts(sce) + 100
assay(sce, "counts_100") <- counts_100 # assign a new entry to assays slot
assays(sce) # new assay has now been added.

# adding metadata to the sce
cell_metadata <- data.frame(batch = c(1, 1, 2))
rownames(cell_metadata) <- paste0("cell_", 1:3)

# can either append metadata to sce or start from scratch 
# starting from scratch here, using the single cell experiment constructor
sce <- SingleCellExperiment(assays = list(counts = counts_matrix),
                            colData = cell_metadata)
colData(sce)
sce$batch

# some functions automatically add column metadata by returning a 
# singlecellexperiment with extra fields in colData slot
# for example: scater::addPerCellQC()

sce <- scater::addPerCellQC(sce)
colData(sce)[, 1:5]
# very useful

# alternatively one might want to add more data manually
sce$more_stuff <- runif(ncol(sce))
colnames(colData(sce))
sce$more_stuff

# colData is very useful vor subsetting the sce object
sce[, sce$batch == 1]
# only two rows returned, those from batch number one

# row data contains DataFrame where each row corresponds to a gene and 
# contains annotations like the transcript length or gene symbol.
# additionally there is a rowRanges slot that holds genomic coordinates 
# in the form of a GRanges data frame or GRanges list
# important to handle chromosome data like the start and end coordinates 
# of genes in the chromosomes. This is easy to query via the GenomicRanges
# framework.

rowRanges(sce)
rowData(sce)
# adding quality control data 
sce <- scater::addPerFeatureQC(sce)
rowData(sce)

# one may also consider adding annotation data manually
# BiocManager::install("AnnotationHub")
# BiocManager::install("ensembldb") # download depends on ensembldb
library(ensembldb)
library(AnnotationHub)
# will download ensembl annotation data to spice up our sce
edb <- AnnotationHub()[["AH73881"]] # Human, Ensembl v97.
genes(edb)[,2]

# subsetting can of course also be performed on the gene level
sce[c("gene_1", "gene_4"), ]
# or equivalently
sce[c(1,4), ]

# other study metadata that may not fit into these specific slots 
# can be stored nicely in the metadata slot; a named list
my_genes <- c("gene_1", "gene_5")
metadata(sce) <- list(favorite_genes = my_genes)
metadata(sce)
metadata(sce)$favorite_genes

your_genes <- c("gene_4", "gene_8")
metadata(sce)$your_genes <- your_genes
metadata(sce)

sce[metadata(sce)$your_genes, ] # nice subsetting is similar to sf objects in geostats

# recap: 
# main data => assays
# cell metadata => colData
# feature metadata => rowData
# all these slots are inherited from the summarizedExperiment class

# adding PCA data
sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
reducedDim(sce, "PCA")

# performing tSNE
sce <- scater::runTSNE(sce, perplexity = 0.1)
reducedDim(sce, "TSNE")

plot(reducedDim(sce, "PCA"))
plot(reducedDim(sce, "TSNE"))

# lets manually add uniform manifold approximation
u <- uwot::umap(t(logcounts(sce)), n_neighbors = 2)
reducedDim(sce, "UMAP_uwot") <- u
reducedDims(sce) # Now stored in the object.
reducedDim(sce, "UMAP_uwot")
plot(reducedDim(sce, "UMAP_uwot"))

# we can also add alternative experiments to the same sce object like so
spike_counts <- cbind(cell_1 = rpois(5, 10), 
                      cell_2 = rpois(5, 10), 
                      cell_3 = rpois(5, 30))
rownames(spike_counts) <- paste0("spike_", 1:5)
spike_se <- SummarizedExperiment(list(counts=spike_counts))
spike_se
# store this Summarized Experiment in the sce
altExp(sce, "spike") <- spike_se
altExps(sce)

# this ensures all relevant data to one experiment can be stored in a single object
# this also uses the same subsetting operations as mentioned above
sub <- sce[,1:2] # retain only two samples.
altExp(sub, "spike")
# dim: 5 2

# the sizeFactors() function allows to get or set a numeric vector of 
# per-cell scaling factors used for normalization. 
sce <- scran::computeSumFactors(sce)
sizeFactors(sce)
# alternatively we can manually add the size factors as shown below for library
# size derived factors

sizeFactors(sce) <- scater::librarySizeFactors(sce)
sizeFactors(sce)

# colLabels() allows to get or set a vector or factor of per cell labels
# typically corresponding to groupings assigned by unsupervised clustering
colLabels(sce) <- LETTERS[1:3]
colLabels(sce)
# does not work here