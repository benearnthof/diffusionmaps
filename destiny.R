if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("destiny")

library(destiny)
# zunder <- read.ExpressionSet("zunder.tab")
# dm <- DiffusionMap(zunder, k = 1000)
# plot(dm, c(2,5,-1), col.by = "Days", pal = blue2green2red(20))
# takes forever on my shitty laptop

data(guo)
dm <- DiffusionMap(guo, verbose = TRUE, censor_val = 15, 
                   censor_range = c(15, 40))
dm
plot(dm)

palette(cube_helix(6))
plot(dm, col_by = "num_cells", legend_main = "Cell stage")
# works