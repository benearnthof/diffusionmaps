# isomap vignette
library(RDRToolbox)
sim = generateData(samples=20, genes=1000, diffgenes=100, diffsamples=10)
simData = sim[[1]]
simLabels = sim[[2]]

res <- Isomap(data=simData, dims=2, k=10)
res <- as.data.frame.matrix(res$dim2)
res$class <- as.factor(simLabels)

ggplot(res, aes(x = V1, y = V2, col = class)) +
  geom_point() +
  scale_fill_manual(c("red", "blue"))

simData <- assays(biase_3_filtered_chosen)$logcounts
simData <- t(simData)
simLabels <- colData(biase_3_filtered_chosen)$cell_type1

res2 <- Isomap(data = simData, dims = 25, k = 10)
res2 <- as.data.frame.matrix(res2$dim2)
res2$class <- as.factor(simLabels)

ggplot(res2, aes(x = V1, y = V2, col = class)) +
  geom_point() +
  ggtitle("Isomap after Feature Selection: Biase") +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00")) +
  labs(color = "Type") +
  theme_bw() +
  ylab("Dim 2") +
  xlab("Dim 1")
