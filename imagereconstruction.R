require(imager)
require(here)
require(dplyr)
require(tidyverse)
require(broom)
require(ggplot2)

img <- load.image(here("karsten2.jpeg"))
str(img)

img_df_long <- as.data.frame(img)

head(img_df_long)

img_df <- tidyr::pivot_wider(img_df_long, 
                             names_from = y, 
                             values_from = value)

# tmp <- img_df[1:800, 1:800]
# cov_karsten <- cov(tmp)
# eigen(cov_karsten)
# eigen_karsten <- eigen(cov_karsten)

dim(img_df)
img_df[1:5, 1:5]
img_pca <- img_df %>%
  dplyr::select(-x,) %>%
  prcomp(scale = TRUE, center = TRUE)

summary(img_pca)

pca_tidy <- tidy(img_pca, matrix = "pcs")

pca_tidy %>%
  ggplot(aes(x = PC, y = percent)) +
  geom_line() +
  labs(x = "Principal Component", y = "Anteil erklärter Varianz") +
  theme_bw() +
  ggtitle("PCA: Anteile erklärter Varianz")

reverse_pca <- function(n_comp = 20, pca_object = img_pca){
  # multiply matrix of rotated data by transpose of matrix of eigenvalues
  # to return original data values
  recon <- pca_object$x[, 1:n_comp] %*% t(pca_object$rotation[, 1:n_comp])
  # reverse scaling and centering
  if(all(pca_object$scale != FALSE)){
    # rescale by the reciprocal of the scaling factor
    recon <- scale(recon, center = FALSE, scale = 1/pca_object$scale)
  }
  if(all(pca_object$center != FALSE)){
    # add subtracted mean to reverse centering
    recon <- scale(recon, scale = FALSE, center = -1 * pca_object$center)
  }
  # create pivotable data frame for imager library
  recon_df <- data.frame(cbind(1:nrow(recon), recon))
  colnames(recon_df) <- c("x", 1:(ncol(recon_df)-1))
  # convert the data to long form 
  recon_df_long <- recon_df %>%
    tidyr::pivot_longer(cols = -x, 
                        names_to = "y", 
                        values_to = "value") %>%
    mutate(y = as.numeric(y)) %>%
    arrange(y) %>%
    as.data.frame()
  
  recon_df_long
}

n_pcs <- c(2, 10, 50, 100)
names(n_pcs) <- paste("Erste", n_pcs, "Komponenten", sep = "_")

# map reverse_pca() 
recovered_imgs <- map_dfr(n_pcs, 
                          reverse_pca, 
                          .id = "pcs") %>%
  mutate(pcs = stringr::str_replace_all(pcs, "_", " "), 
         pcs = factor(pcs, levels = unique(pcs), ordered = TRUE))

# reconstructing original image 
p <- ggplot(data = recovered_imgs, 
            mapping = aes(x = x, y = y, fill = value))
p_out <- p + geom_raster() + 
  scale_y_reverse() + 
  scale_fill_gradient(low = "black", high = "white") +
  facet_wrap(~ pcs, ncol = 2) + 
  guides(fill = FALSE) + 
  labs(title = "Rekonstruktion von Karsten aus den ersten \nN Hauptkomponenten") + 
  theme(strip.text = element_text(face = "bold", size = rel(1)),
        plot.title = element_text(size = rel(1)),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab("") +
  ylab("")


p_out

ggsave("karstenreconstructed.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 10, height = 10, units = "cm",
       dpi = 300, limitsize = TRUE)

