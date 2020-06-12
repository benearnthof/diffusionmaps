# run create_sce.R first then set wd to this folder
### DATA
d <- read.table("data.txt", header = T)
d <- d[!duplicated(d$ID),]
rownames(d) <- d$ID
d <- d[,2:ncol(d)]

### ANNOTATIONS
ann <- read.table("biase_cell_types.txt", stringsAsFactors = F)

### SINGLECELLEXPERIMENT
source("../create_sce.R")
sceset <- create_sce_from_normcounts(d, ann)
# convert ensembl ids into gene names
# gene symbols will be stored in the feature_symbol column of fData
sceset <- getBMFeatureAnnos(
  sceset, filters="ensembl_gene_id",
  biomart="ensembl", dataset="mmusculus_gene_ensembl")
# remove features with duplicated names
sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]

saveRDS(sceset, "biase.rds")

biase <- readRDS("biase.rds")

library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
chr.loc <- mapIds(ens.mm.v97, keys=rownames(biase),
                  keytype="GENEID", column="SEQNAME")
is.mito <- which(chr.loc=="MT")

assays(biase)$counts <- assays(biase)$normcounts
df <- perCellQCMetrics(biase, subsets=list(Mito=is.mito))
df

biase_qc <- addPerCellQC(biase, subsets=list(Mito=is.mito))
colnames(colData(biase_qc))

# identifying low quality gene reads
reasons <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent"))
colSums(as.matrix(reasons))

# four cells will be discarded
# First attempt with batch-specific thresholds.
colnames(colData(biase_qc))
discard.mito <- isOutlier(biase_qc$subsets_Mito_percent, type="higher")
# only two would be discarded
plotColData(biase_qc, x="cell_type1", y="subsets_Mito_percent",
                             colour_by=I(discard.mito))

####### to understand needed structure for plotColData
# library(scRNAseq)
# sce.grun <- GrunPancreasData()
# sce.grun <- addPerCellQC(sce.grun)
# # First attempt with batch-specific thresholds.
# # is.outlier takes colnames of coldata
# # plotColData too
# discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
#                           type="higher", batch=sce.grun$donor)
# with.blocking <- plotColData(sce.grun, x="donor", y="altexps_ERCC_percent",
#                              colour_by=I(discard.ercc))

# checking diagnostic plots
biase_3 <- biase_qc

colData(biase_3) <- cbind(colData(biase_3), df)
colnames(colData(biase_3))
biase_3$Type <- factor(biase_3$cell_type1)
biase_3$Discard <- reasons$discard

plotColData(biase_3, x="Type", y="sum", colour_by="Discard") + 
  scale_y_log10() + 
  ggtitle("Total Count of Reads per Cell: Biase") +
  ylab("Total Count") +
  scale_fill_manual(values = c("#0000FF", "#ff0000")) +
  labs(fill = "Discard")

plotColData(biase_3, x="Type", y="detected", colour_by="Discard") + 
  scale_y_log10() + 
  ggtitle("Detected Features per Cell: Biase") +
  ylab("Detected Features") +
  scale_fill_manual(values = c("#0000FF", "#ff0000")) +
  labs(fill = "Discard")

plotColData(biase_3, x="Type", y="subsets_Mito_percent", 
            colour_by="Discard") +
  ggtitle("Percentage of Mitochondrial RNA: Biase") +
  ylab("Mitochondrial RNA in %") +
  scale_fill_manual(values = c("#0000FF", "#ff0000")) +
  labs(fill = "Discard")

# trying out droplet based diagnostics
library(DropletUtils)
library(Matrix)

bcrank <- barcodeRanks(counts(biase_3))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

biase_3_filtered <- biase_3[,!reasons$discard]
biase_3 #56
biase_3_filtered #52

# normalization 
# library size normalization
lib_sf_biase <- librarySizeFactors(biase_3)
lib_sf_biase_filtered <- librarySizeFactors(biase_3_filtered)
summary(lib_sf_biase)
summary(lib_sf_biase_filtered)
# looks promising so far
hist(lib_sf_biase, xlab="Size factor", col='grey80', breaks = 20)
hist(lib_sf_biase_filtered, xlab="Size factor", col='grey80', breaks = 20)
hist(log10(lib_sf_biase), xlab="Log10[Size factor]", col='grey80', breaks = 20)

pca <- runPCA(biase_3)
colnames(colData(pca))
plt <- scater::plotPCA(pca, colour_by = "cell_type1") +
  ggtitle("PCA without Feature Selection: Biase") +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00")) +
  labs(fill = "Type")
  
plt

# Feature Selection
# quantifying per gene variation 
library(scran)
tmp <- modelGeneVar(biase_3)
fit <- metadata(tmp)
plot(fit$mean, fit$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression", main = "Quantifying per Gene Variation: Biase")
curve(fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# coefficient of variation 
cv2 <- modelGeneCV2(biase_3)
fit_cv2 <- metadata(cv2)
plot(fit_cv2$mean, fit_cv2$cv2, log = "xy")
curve(fit_cv2$trend(x), col="dodgerblue", add=TRUE, lwd=2)
# produces nonsensical plot

# modeling gene var with poisson model
set.seed(0010101)
pois <- modelGeneVarByPoisson(biase_3)

plot(pois$mean, pois$total, pch=16, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(metadata(pois)$trend(x), col="dodgerblue", add=TRUE)

# getting the top 1000 active genes (high variance genes)
hvg <- getTopHVGs(tmp, n=1000)
str(hvg)

# getting all above the line
hvg_2 <- getTopHVGs(tmp, var.threshold=0)
length(hvg_2)
# 10000 seems like a good tradeoff; less than half of original genes, yet 
# still more than enough to use to give secondary structure a chance to manifest 
# itself

biase_3_chosen <- biase_3[hvg_2,]
tmp_filtered <- modelGeneVar(biase_3_filtered)
hvg_2_filtered <- getTopHVGs(tmp_filtered, var.threshold = 0)
biase_3_filtered_chosen <- biase_3_filtered[hvg_2_filtered,]

dim(biase_3_chosen)
dim(biase_3_filtered_chosen)
# seems good so far

saveRDS(biase_3_chosen, file = "biase_chosen.RDS")
saveRDS(biase_3_filtered_chosen, file = "biase_filtered_chosen.RDS")

library(scater)
library(SingleCellExperiment)

biase_3_chosen <- readRDS("biase_chosen.RDS")
biase_3_filtered_chosen <- readRDS("biase_filtered_chosen.RDS")

# performing pca on chosen set
pca2 <- runPCA(biase_3_chosen)

plt2 <- scater::plotPCA(pca2, colour_by = "cell_type1") +
  ggtitle("PCA after Feature Selection: Biase") +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00")) +
  labs(fill = "Type")

plt2

pca3 <- runPCA(biase_3_filtered_chosen)

plt3 <- scater::plotPCA(pca3, colour_by = "cell_type1") +
  ggtitle("PCA after Feature Selection: Biase") +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00")) +
  labs(fill = "Type")

plt3

# Percentage of variance explained is tucked away in the attributes.
percent.var <- attr(reducedDim(pca3), "percentVar")
# BiocManager::install("PCAtools")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow

plot(percent.var[1:25], xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")

# applying dimensionality reduction
reducedDim(pca3, "PCA") <- reducedDim(pca3, "PCA")[,1:10]
ncol(reducedDim(pca3, "PCA"))

# diffusionmap
set.seed(1100101001)
diffmap1 <- runDiffusionMap(biase_3_filtered_chosen, ncomponents = 3)
plotReducedDim(diffmap1, dimred="DiffusionMap", colour_by="cell_type1", ncomponents = 3) +
  ggtitle("Diffusion Map after Feature Selection: Biase") +
  scale_fill_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00")) +
  labs(fill = "Type")

# for isomap => separate file
# BiocManager::install("RDRToolbox")
# library(RDRToolbox)
# casuses my laptop to crash
# test <- Isomap(assays(biase_3_filtered_chosen)$normcounts)

library(destiny)
?destiny

dm <- destiny::DiffusionMap(biase_3_filtered_chosen)
plot(dm, c(1,2))

dm@d

all.equal(rownames(dm@eigenvectors),colnames(biase_3_filtered_chosen))

eigs <- dm@eigenvectors
Diffmap_eucl <- dm

D_eucl <- dist(eigs)

# ja D_eucl ist die euklidische Distanzmatrix. 
# Diese muss jedoch als data.frame übergeben werden.
# res_diff_eucl müsste ein data.frame sein. 
# (Die Anzahl Spalten sollten Anzahl single-cells 
# des jeweiligen Datensatzes sein)

n_DC <- floor(0.04 * dim(D_eucl)[2]):ceiling(0.07 * dim(D_eucl)[2]) 
n_DC  <- sort(n_DC)

res_diff_eucl <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = dim(D_eucl)[1]))
res_diff_pear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = dim(D_eucl)[1]))
res_diff_spear <- as.data.frame(matrix(NA, nrow= n_DC[length(n_DC)], ncol = dim(D_eucl)[1]))


Diffmap_eucl <- destiny::DiffusionMap(biase_3_filtered_chosen)
Diffmap_pear <- destiny::DiffusionMap(biase_3_filtered_chosen, distance = "cosine")
Diffmap_spear <- destiny::DiffusionMap(biase_3_filtered_chosen, distance = "rankcor")
# biase: ktrue = 3
k_true <- 3
truth <- as.factor(colData(biase_3_filtered_chosen)$cell_type1)

for(j in n_DC[1]:n_DC[length(n_DC)]){ 
  
  tmp_eucl <- kmeans(Diffmap_eucl@eigenvectors[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
  res_diff_eucl[n_DC[length(n_DC)]-j+1, ] <- tmp_eucl$cluster
  tmp_pear <- kmeans(Diffmap_pear@eigenvectors[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
  res_diff_pear[n_DC[length(n_DC)]-j+1, ] <- tmp_pear$cluster
  tmp_spear <- kmeans(Diffmap_pear@eigenvectors[,1:j], centers = k_true, iter.max = 10^9, nstart = 1000)
  res_diff_spear[n_DC[length(n_DC)]-j+1, ] <- tmp_pear$cluster
}
#                                   Default value                 cosine distance         rank correlation distance
# Die Diffusionmaps jeweils auf der Euklidischen Distanz (D_eucl),1- Pearsonkorrelation , und 1-Spearmankorrelation anwenden

# Dann: 
neu_diff <- rbind( res_diff_eucl , res_diff_pear , res_diff_spear)
sc3_diff <- rbind(neu_diff)
cl_kmeans_dmap <- na.omit(sc3_diff)

consens_mat <- as.data.frame(matrix(NA, nrow = dim(D_eucl)[2], ncol = dim(D_eucl)[2]))
for(i in 1:dim(D_eucl)[2]){
  for(j in 1:dim(D_eucl)[2]){
    consens_mat[i,j]<- sum( cl_kmeans_dmap [,i]== cl_kmeans_dmap [,j])
    consens_mat[j,i]<- sum( cl_kmeans_dmap [,i]== cl_kmeans_dmap [,j])
  }
}

hc <- hclust(dist(consens_mat))

k_true  <- length(levels(truth))

res_2 <- cutree(hc, k=k_true)
mclust::adjustedRandIndex(truth, res_2)

# truth und k_true sind die vorliegenden Zelltypen sowie die jeweilige vorliegende Anzahl Kategorien von Zelltypen des jeweiligen Datensatzes

# BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
# input data format is a matrix where columns are samples and rows are features
# in this case the columns are the cells, the rows are the reduced dimensions
dm <- destiny::DiffusionMap(biase_3_filtered_chosen, n_eigs = 10)
plot(dm, c(1,2,3))
dim(dm@eigenvectors)

eigs <- dm@eigenvectors
rownames(eigs) <- colData(biase_3_filtered_chosen)$cell_type1

# columns are cells, rows are diffusion components
dta <- t(eigs)
# modified head function to display only first n rows AND n columns
head2d = function(x, n = 6L, ...) {
  indvecs = lapply(seq_along(dim(x)), function(i) {
    if (length(n) >= i) {
      ni = n[i]
    } else {
      ni =  dim(x)[i]
    }
    if (ni < 0L)
      ni = max(nrow(x) + ni, 0L)
    else
      ni = min(ni, dim(x)[i])
    seq_len(ni)
  })
  lstargs = c(list(x),indvecs, drop = FALSE)
  do.call("[", lstargs)
}

head2d(dta, n = c(5,5))
# seems good

# running consensus clustering
title = tempdir()

# install.packages("amap")
# library(amap)
# dta_pearson <- amap::Dist(dta, method = "pearson")
res <- ConsensusClusterPlus(dta, maxK = 4, reps = 100, pItem = 0.8, pFeature = 0.8,
                            title = title, clusterAlg = "kmdist", 
                            distance = "pearson", seed = 1,plot = "png")

# 5 entries with euclidean seems to work relatively well
# 10 entries with pearson distance seems to work really well
# 10 entries with spearman distance seems to work reasonably well

# generate cluster consensus and item consensus like so
icl = calcICL(res, title = title, plot="png")

tmp <- res[[3]]$consensusClass

mclust::adjustedRandIndex(truth, tmp)

# steps to repeat on all datasets: 
# do preprocessing
# run all 3 dim reduction algorithms
# select 5 10 15 20 25 dims from reduced dims
# run clustering on that with 3 different distance metrics
# repeat consensus clustering 100 times each 
# calculate ARI for everything
# plot results in pretty boxplot

dataset <- biase_3_filtered_chosen

# do methods of dim reduction and save 50 dimensions
# should all have the same structure
library(RDRToolbox) # for isomap

test <- "cell_type1"
labs <- subset(colData(dataset), select = test)[,1]

difres <- destiny::DiffusionMap(dataset, n_eigs = 25)@eigenvectors
rownames(difres) <- labs
difres <- t(difres)

pcares <- reducedDim(scater::runPCA(dataset, ncomponents = 25))
rownames(pcares) <- labs
pcares <- t(pcares)

isodata <- assays(dataset)$logcounts
isodata <- t(isodata)

isores <- RDRToolbox::Isomap(data = isodata, dims = 25, k = 10)
isores <- as.matrix(isores$dim25)
rownames(isores) <- labs
isores <- t(isores)

lapply(list(difres, pcares, isores), head2d, n = c(5,5))
# 3 matrices with the needed structure

# lets wrap that into a function
dimred_processing <- function(dataset, labs, n = 25) {
  labs <- subset(colData(dataset), select = labs)[,1]
  # diffusion mapping
  difres <- destiny::DiffusionMap(dataset, n_eigs = n)@eigenvectors
  rownames(difres) <- labs
  difres <- t(difres)
  # pca
  pcares <- reducedDim(scater::runPCA(dataset, ncomponents = n))
  rownames(pcares) <- labs
  pcares <- t(pcares)
  # isomap
  isodata <- assays(dataset)$logcounts
  isodata <- t(isodata)
  isores <- RDRToolbox::Isomap(data = isodata, dims = n, k = 10)
  isores <- as.matrix(isores[[1]])
  rownames(isores) <- labs
  isores <- t(isores)
  ret <- list(difres = difres, pcares = pcares, isores = isores)
  return(ret)
}

dta <- readRDS("biase_filtered_chosen.RDS")

rdims <- dimred_processing(dta, "cell_type1", n = 5)

library(purrr)
map(rdims, head2d, n = c(5,5))

# now we can use purrs listprocessing to avoid unwrapping lists
set.seed(1)

truth <- as.factor(colData(dta)$cell_type1)

calc_scores <- function(rdims, nclust, dist = c("euclidean", "pearson", "spearman"),
                        truth, iters = 25) {
  scores <- data.frame(difres = rep(0, iters),
                       pcares = rep(0, iters), 
                       isores = rep(0, iters))
  for (i in 1:iters) {
    # setting pFeatures = 1 so we dont exclude reduced dimensions
    # setting pItem to 0.8 allows resampling from cells
    ccp <- map(rdims, ConsensusClusterPlus, maxK = nclust, reps = 5, pItem = 0.8, 
               pFeature = 0.9, clusterAlg = "kmdist", distance = dist)
    # getting results for cluster number of interest
    classes <- map(ccp, `[[`, nclust)
    classes <- map(classes, `$`, "consensusClass")
    # calculating scores and storing them
    score <- map(classes, mclust::adjustedRandIndex, truth)
    score <- t(as.matrix(score))
    scores[i,] <- score
  }
  return(scores)
}

biase_eucl <- calc_scores(rdims, 3, dist = "euclidean", truth = truth)
biase_eucl2 <- calc_scores(rdims, 3, dist = "euclidean", truth = truth)
biase_eucl3 <- calc_scores(rdims, 3, dist = "euclidean", truth = truth)
biase_eucl4 <- calc_scores(rdims, 3, dist = "euclidean", truth = truth)

plt_scores <- function(scores) {
  mlt <- reshape2::melt(scores)
  plt <- ggplot(mlt, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    xlab("Method") +
    ylab("Adjusted Rand Index") +
    ggtitle("Adjusted Rand Index: Biase", subtitle = "Euclidean Distance") +
    scale_x_discrete(labels = c("difres" = "Diffusion Map", "pcares" = "PCA",
                                "isores" = "Isomap")) +
    theme_bw() +
    theme(legend.position = "none")
  return(plt)
}

plt_scores(biase_eucl4)

biase_spear <- calc_scores(rdims, 3, dist = "spearman", truth = truth)
plt_scores(biase_spear)

biase_pears <- calc_scores(rdims, 3, dist = "pearson", truth = truth)
plt_scores(biase_pears)

biase_scores_5dim <- list(eucl = biase_eucl4, spear = biase_spear, pears = biase_pears)

saveRDS(biase_scores_5dim, file = "biase_scores_5dim.RDS")

# looks okay so far



