if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("destiny")

library(destiny)
# https://www.helmholtz-muenchen.de/fileadmin/ICB/software/destiny/zunder.tab.gz
zunder <- read.ExpressionSet("zunder.tab")
dm_zunder <- DiffusionMap(zunder, k = 2, verbose = TRUE)
plot(dm, c(2,5,-1), col.by = "Days", pal = blue2green2red(20))
# takes forever on my shitty laptop

data(guo)
dm <- DiffusionMap(guo, verbose = TRUE, censor_val = 15,
                   censor_range = c(15, 40))
dm
plot(dm)

palette(cube_helix(6))
plot(dm, col_by = "num_cells", legend_main = "Cell stage")
# works

covars <- data.frame(covar1 = letters[1:100])
dists <- dist(matrix(rnorm(100*10), 100))
dm_letters <- DiffusionMap(covars, distance = dists)
plot(dm_letters)

# visualization of the distribution of eigenvalues for a random covariance matrix
# to get an understanding of spectral decay
mat <- matrix(nrow = 1000, ncol = 1000)
mat[] <- rnorm(1e6)
cov <- cov(mat)
eigen(cov)
eigen <- eigen(cov)
plot(eigen$values)
plot(ecdf(eigen$values))
plot(density(eigen$values))
hist(eigen$values, breaks = 1000)

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.1, 0, 0.1), c("darkblue","white", "red"), transparency = 0.5)
Heatmap(cov, name = "Covariance", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE)
