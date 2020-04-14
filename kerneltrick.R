# illustrating the kernel trick
x <- matrix(rnorm(600),nc=2)
# y contains circle
y <- x/sqrt(rowSums(x^2))

plot(y[,1] ~ y[,2])

dta <- data.frame(x = y[,1], y = y[,2])
dta$x <- dta$x + rnorm(300, mean = 0, sd = 0.05)
dta$y <- dta$y + rnorm(300, mean = 0, sd = 0.05)
dta$class = rep(1, times = 300)

dot <- data.frame(class = rep(2, times = 200))
dot$x = rnorm(200, mean = 0, sd = 0.1)
dot$y = rnorm(200, mean = 0, sd = 0.1)

plot(dot$x~dot$y)

dta_full <- rbind(dta, dot)
plot(dta_full$y~dta_full$x)

dta_full$class <- as.factor(dta_full$class)
library(ggplot2)
ggplot(dta_full, aes(x = x, y = y, col = class)) +
  geom_point() +
  scale_color_manual(values = c('#FF0000', '#0000FF')) + 
  theme_bw() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("")

dta_full$z <- dta_full$x^2 + dta_full$y^2

scatterplot3d::scatterplot3d(x = dta_full$x, y = dta_full$y, z = dta_full$z)
library(plotly)
plot_ly(x = dta_full$x, y = dta_full$y, z = dta_full$z, mode = "markers",
        type = "scatter3d", marker = list(color = dta_full$class, 
                                          colorscale = c('#FF0000', '#0000FF'),
                                          showscale = TRUE))

plot_ly(dta_full, x = ~x, y = ~y, z = ~z, color = ~class, colors = c('#FF0000', '#0000FF'))
