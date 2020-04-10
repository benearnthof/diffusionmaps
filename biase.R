# biase
### DATA
d <- read.table("biase.txt", header = T)
d <- d[!duplicated(d$ID),]
rownames(d) <- d$ID
d <- d[,2:ncol(d)]

### ANNOTATIONS
ann <- read.table("biase_cell_types.txt", stringsAsFactors = F)

### SINGLECELLEXPERIMENT
# source("../utils/create_sce.R")
sceset <- create_sce_from_normcounts(d, ann)
# convert ensembl ids into gene names
# gene symbols will be stored in the feature_symbol column of fData
sceset <- getBMFeatureAnnos(
  sceset, filters="ensembl_gene_id",
  biomart="ensembl", dataset="mmusculus_gene_ensembl")
# remove features with duplicated names
sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
saveRDS(sceset, "biase.rds")
sceset <- readRDS("biase.rds")

require(destiny)
dm <- DiffusionMap(sceset, n_pcs = 50)
plot(dm)

dm2 <- DiffusionMap(sceset)
plot(dm2)

dm_sigma <- DiffusionMap(sceset, sigma = "global")
plot(dm_sigma)
plot(DiffusionMap(sceset, sigma = 120, n_pcs = 50))
plot(DiffusionMap(sceset, sigma = "local", n_pcs = 50, n_eigs = 2))
plot(DiffusionMap(sceset, sigma = "local", n_pcs = 50, n_eigs = 2, 
                  density_norm = FALSE))
plot(DiffusionMap(sceset, sigma = "local", n_pcs = 50, n_eigs = 2,
                  distance = "cosine"))

plot(DiffusionMap(sceset, sigma = "local", n_pcs = NULL, n_eigs = 2))
plot(DiffusionMap(sceset, sigma = "local", n_pcs = NULL, n_eigs = 2, n_local = 55))
plot(DiffusionMap(sceset, sigma = "local", n_pcs = NULL, n_eigs = 2, rotate = FALSE))

plot(DiffusionMap(sceset, sigma = "local", n_pcs = NULL, n_eigs = 2), col_by = )

# # lets try other package
# df <- sceset@assays@data@listData$normcounts
# D <- dist(df)
# dm <- diffuse(D)
# # uses 50 GB of RAM lmao

data(guo)
plot(DiffusionMap(guo, n_eigs = 2))
plot(DiffusionMap(guo, 13, censor_val = 15, censor_range = c(15, 40), verbose = TRUE))

plotColData

# diffusion map plot uses louvain clustering by default
g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
g <- add_edges(g, c(1,6, 1,11, 6, 11))
plot(g)
cluster_louvain(g)

install.packages("biocLite")
biocLite("scater")
