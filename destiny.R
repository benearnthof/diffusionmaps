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
mat[] <- rnorm(1000*1000)
cov <- cov(mat)
eigen(cov)
eigen <- eigen(cov)
plot(eigen$values)
plot(ecdf(eigen$values))
plot(density(eigen$values))
hist(eigen$values, breaks = 100)
lines(y~x)

devtools::install_github("kassambara/ggcorrplot")
library(ggcorrplot)
ggcorrplot(cov)

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.2, 0, 0.2), c("darkblue","green", "red"), transparency = 0.5)
Heatmap(cov, name = "Covariance", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE)

x <- seq(from = 0.15, to = 4, length.out = 1000)
y <- (sqrt(4*x-x^2))/(2*pi*x)

# karsten data from imagereconstruction.R file
df <- data.frame(val = eigen$values, y = y, x = x)
ggplot(df, aes(x = val)) +
  geom_density(color = "blue", size = 1) +
  # geom_density(aes(x = karsten), color = "green", size = 1) + dominates entire plot
  geom_line(aes(x = x, y = y), color = "red", size = 1) +
  theme_bw() +
  ggtitle("Approx. Dichte der Eigenwerte einer zufälligen Covarianzmatrix (Blau) \n Dichte der Marchenko-Pastur Verteilung (Rot)") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

theta <- seq(from = 0, to = 4*pi, by = pi/29)
x <- 2 * cos(theta)
y <- 2 * sin(theta)
z <- theta/4
plot(x ~ y)
library(scatterplot3d)
df <- data.frame(x = x, y = y, z = z)
scatterplot3d(df)

# install.packages("plotly")
library(plotly)
p <- plot_ly(df, x = ~x, y = ~y, z = ~z)
p

type <- 1:nrow(df)
df$type = type
p <- plot_ly(df, x = ~x, y = ~y, z = ~z, marker = list(color = ~type, colorscale = c('#FFE1A1', '#683531'))) %>%
  add_markers()
p

dm <- DiffusionMap(df[,1:3], n_eigs = 10)
plot(dm)
plot(dm$DC2 ~ dm$DC1)
dc <- data.frame(DC2 = dm$DC2, DC1 = dm$DC1, DC3 = dm$DC3, type = type)
plot_ly(dc, x = ~DC1, y = ~DC2, z = ~DC3, marker = list(color = ~type, colorscale = c('#FFE1A1', '#683531'))) %>%
  add_markers()

# lets see other options
install.packages("diffusionMap")
library(diffusionMap)
n <- 2000
t <- runif(n) ^ 0.7 * 10
al <- 0.15
bet <- 0.5
x1 <- bet * exp(al * t) * cos(t) + rnorm(n, 0, 0.1)
y1 <- bet * exp(al * t) * sin(t) + rnorm(n, 0, 0.1)
plot(x1, y1, pch = 20, main = "Noisy spiral")
# compute diffusion map
d <- dist(cbind(x1, y1))
dmap <- diffuse(d, neigen = 10)
plot(dmap)
plot(t, dmap$X[,1], pch = 20, axes = FALSE, xlab = "spiral parameter", ylab = "1st diffusion coefficient")
box()
plot(1:10, dmap$eigenmult, typ='h', xlab="diffusion map dimension", ylab= "eigen-multipliers")

type = 1:2000
df <- data.frame(x = dmap$X[,1], y = dmap$X[,2], z = dmap$X[,3], type = type)
plot_ly(df, x = ~x, y = ~y, z = ~z, marker = list(color = ~type, colorscale = c('#FFE1A1', '#683531'))) %>%
  add_markers()

data(annulus)
plot(annulus,main="Annulus Data",pch=20,cex=.7)
D = dist(annulus) # use Euclidean 
distancedmap = diffuse(D,eps.val=.1) 
# compute diffusion map & plot 
print(distancedmap)
plot(distancedmap)

par(mfrow = c(2, 1))
plot(annulus, main = "Annulus Data", pch = 20, cex = 0.7)
dmap = diffuse(D, eps.val = 0.05)
k = 2 # number of clusters
dkmeans <- diffusionKmeans(dmap, k)
plot(annulus, main = "Colored by diffusion K-means clustering", pch = 20, cex = 0.7, col = dkmeans$part)

# lets try with other example data set
data("Chainlink")
lab.col = c(rep("red", 500), rep("blue", 500))
n = 1000
scatterplot3d(Chainlink$C1, Chainlink$C2, Chainlink$C3, color = lab.col, main = "Chainlink Data")
D <- dist(Chainlink)
dmap = diffuse(D, neigen = 3, eps.val = 0.1)
plot(dmap)
dkmeans <- diffusionKmeans(dmap, K = 2)
col.dkmeans <- ifelse(dkmeans$part == 1, "red", "blue")
scatterplot3d(Chainlink, color = col.dkmeans, main = "Chainlink Data")
table(dkmeans$part, lab.col)
