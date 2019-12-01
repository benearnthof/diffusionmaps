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
