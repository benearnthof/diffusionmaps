# mnist digit recognition data set
# download data from http://yann.lecun.com/exdb/mnist/
download.file("http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz",
              "train-images-idx3-ubyte.gz")
download.file("http://yann.lecun.com/exdb/mnist/train-labels-idx1-ubyte.gz",
              "train-labels-idx1-ubyte.gz")
download.file("http://yann.lecun.com/exdb/mnist/t10k-images-idx3-ubyte.gz",
              "t10k-images-idx3-ubyte.gz")
download.file("http://yann.lecun.com/exdb/mnist/t10k-labels-idx1-ubyte.gz",
              "t10k-labels-idx1-ubyte.gz")

# gunzip the files
R.utils::gunzip("train-images-idx3-ubyte.gz")
R.utils::gunzip("train-labels-idx1-ubyte.gz")
R.utils::gunzip("t10k-images-idx3-ubyte.gz")
R.utils::gunzip("t10k-labels-idx1-ubyte.gz")

# helper function for visualization
show_digit = function(arr784, col = gray(12:1 / 12), ...) {
  image(matrix(as.matrix(arr784[-785]), nrow = 28)[, 28:1], col = col, ...)
}

# load image files
load_image_file = function(filename) {
  ret = list()
  f = file(filename, 'rb')
  readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  n    = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  nrow = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  ncol = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  x = readBin(f, 'integer', n = n * nrow * ncol, size = 1, signed = FALSE)
  close(f)
  data.frame(matrix(x, ncol = nrow * ncol, byrow = TRUE))
}

# load label files
load_label_file = function(filename) {
  f = file(filename, 'rb')
  readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  n = readBin(f, 'integer', n = 1, size = 4, endian = 'big')
  y = readBin(f, 'integer', n = n, size = 1, signed = FALSE)
  close(f)
  y
}

# load images
train = load_image_file("train-images-idx3-ubyte")
test  = load_image_file("t10k-images-idx3-ubyte")

# load labels
train$y = as.factor(load_label_file("train-labels-idx1-ubyte"))
test$y  = as.factor(load_label_file("t10k-labels-idx1-ubyte"))

# view test image
show_digit(train[123, ])

# testing classification on subset of training data
fit = randomForest::randomForest(y ~ ., data = train[1:1000, ])
fit$confusion
test_pred = predict(fit, test)
mean(test_pred == test$y)
table(predicted = test_pred, actual = test$y)

dta <- test[,1:784]
dta <- test
dta <- dta[1:500,]

library(destiny)
library(diffusionMap)

D <- stats::dist(dta)
dm <- diffuse(D, t = 10)
plot(dm)

type = test$y[1:500]
library(magrittr)
library(plotly)
df <- data.frame(x = dm$X[,1], y = dm$X[,2], z = dm$X[,3], type = type)
plot_ly(df, x = ~x, y = ~y, z = ~z, marker = list(color = ~type, colorscale = "Viridis", name = ~type)) %>%
  add_markers()



dm <- DiffusionMap(dta[1:10,], dist = D, sigma = "local")
data(guo)
guo
require(Biobase)
object <- new("ExpressionSet", exprs = as.matrix(dta))

dm <- DiffusionMap(object, k = 49)
require(colorRamps)
plot(dm, col.by = "y", pal = blue2green2red(10))

dc <- data.frame(DC2 = dm$DC2, DC1 = dm$DC1, DC3 = dm$DC3, type = type)
plot_ly(dc, x = ~DC1, y = ~DC2, z = ~DC3, marker = list(color = ~type, colorscale = c('#FFE1A1', '#683531'))) %>%
  add_markers()

library(readr)
library(Rtsne)


# The competition datafiles are in the directory ../input
# Read competition data files:

train$y <- as.factor(train$y)

# shrinking the size for the time limit
numTrain <- 5000
set.seed(1)
rows <- sample(1:nrow(train), numTrain)
train2 <- train[rows,]
# using tsne
set.seed(1) # for reproducibility
tsne <- Rtsne(train2[,-785], dims = 3, perplexity=30, verbose=TRUE, max_iter = 500)
# visualizing
colors = rainbow(length(unique(train2$y)))
names(colors) = unique(train2$y)
      plot(tsne$Y, t='n', main="tsne")
      text(tsne$Y, labels=train2$y, col=colors[train2$y])

  # compare with pca
pca = princomp(train2[,-785])$scores[,1:2]
plot(pca, t='n', main="pca")
text(pca, labels=train2$y,col=colors[train2$y])

library(scatterplot3d)

scatterplot3d(x=tsne$Y[,1],y=tsne$Y[,2],z=tsne$Y[,3],
              color = colors[train2$y])

library(rgl)
library(magick)

open3d()
par3d(windowRect = c(20, 30, 500, 500))
plot3d(x=tsne$Y[,1],y=tsne$Y[,2],z=tsne$Y[,3], 
       col=colors[train2$y], type="s",radius=0.5, xlab = "", ylab = "", zlab = "")
movie3d(spin3d(axis = c(1, 1, 1), rpm = 2), duration = 30, dir = getwd(), )


# trying uniform manifold approximation 
library(uwot)

# See function man page for help
?umap

# Non-numeric columns are ignored, so in a lot of cases you can pass a data
# frame directly to umap
iris_umap <- umap(iris, n_neighbors = 50, learning_rate = 0.5, init = "random")

# Load mnist from somewhere, e.g.
devtools::install_github("jlmelville/snedata")
mnist <- snedata::download_mnist()
mnist_umap <- umap(mnist, n_neighbors = 15, min_dist = 0.001, verbose = TRUE)
devtools::install_github("jlmelville/vizier")
embed_img <- function(X, Y, k = 15, ...) {
  args <- list(...)
  args$coords <- Y
  args$x <- X
  
  do.call(vizier::embed_plot, args)
}
embed_img(iris, iris_umap, pc_axes = TRUE, equal_axes = TRUE, alpha_scale = 0.5, title = "iris UMAP", cex = 1)
embed_img(mnist, mnist_umap, pc_axes = TRUE, equal_axes = TRUE, alpha_scale = 0.5, title = "MNIST UMAP", cex = 1)
# For high dimensional datasets (> 100-1000 columns) using PCA to reduce
# dimensionality is highly recommended to avoid the nearest neighbor search
# taking a long time. Keeping only 50 dimensions can speed up calculations
# without affecting the visualization much
mnist_umap_pca <- umap(mnist, pca = 50)
embed_img(mnist, mnist_umap_pca, pc_axes = TRUE, equal_axes = TRUE, alpha_scale = 0.5, title = "MNIST UMAP PCA", cex = 1)

# Use a specific number of threads
# mnist_umap <- umap(mnist, n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8)

# Use a different metric
mnist_umap_cosine <- umap(mnist, n_neighbors = 15, metric = "cosine", min_dist = 0.001, verbose = TRUE, n_threads = 8)

# If you are only interested in visualization, `fast_sgd = TRUE` gives a much faster optimization
mnist_umap_fast_sgd <- umap(mnist, n_neighbors = 15, metric = "cosine", min_dist = 0.001, verbose = TRUE, fast_sgd = TRUE)

# Supervised dimension reduction
mnist_umap_s <- umap(mnist, n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8,
                     y = mnist$Label, target_weight = 0.5)

# Add new points to an existing embedding
mnist_train <- head(mnist, 60000)
mnist_test <- tail(mnist, 10000)

# You must set ret_model = TRUE to return extra data we need
# coordinates are in mnist_train_umap$embedding
mnist_train_umap <- umap(mnist_train, verbose = TRUE, ret_model = TRUE)
mnist_test_umap <- umap_transform(mnist_test, mnist_train_umap, verbose = TRUE)

# Save the nearest neighbor data
mnist_nn <- umap(mnist, ret_nn = TRUE)
# coordinates are now in mnist_nn$embedding

# Re-use the nearest neighor data and save a lot of time
mnist_nn_spca <- umap(mnist, nn_method = mnist_nn$nn, init = spca)

# No problem to have ret_nn = TRUE and ret_model = TRUE at the same time
# Or just use the ret_extra parameter:
mnist_nn_and_model <- umap(mnist, ret_extra = c("model", "nn"))

# You can also get to the input fuzzy graph as a sparse matrix via "fgraph"
mnist_with_fgraph <- umap(mnist, ret_extra = c("fgraph"))
# equivalent for lvish is to use "P" (input probability matrix):
mnist_with_P <- lvish(mnist, ret_extra = c("P"))

# Calculate Petal and Sepal neighbors separately (uses intersection of the resulting sets):
iris_umap <- umap(iris, metric = list("euclidean" = c("Sepal.Length", "Sepal.Width"),
                                      "euclidean" = c("Petal.Length", "Petal.Width")))
# Can also use individual factor columns
iris_umap <- umap(iris, metric = list("euclidean" = c("Sepal.Length", "Sepal.Width"),
                                      "euclidean" = c("Petal.Length", "Petal.Width"),
                                      "categorical" = "Species"))