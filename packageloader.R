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
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = counts_matrix))

BiocManager::install("destiny")
library(destiny)
plot()