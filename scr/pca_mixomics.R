# AIM ---------------------------------------------------------------------
# PCA on Clinical Dataset using mixOmics

# lirbaries ---------------------------------------------------------------
library(tidyverse)
library(mixOmics)

# 1. Load data ------------------------------------------------------------
data <- read.csv("data/clinical_data.csv", row.names = "patient_id")

# 2. Run PCA --------------------------------------------------------------
# ncomp: number of principal components to compute
pca_res <- pca(data, ncomp = 6, center = TRUE, scale = TRUE)
print(pca_res)

# 3. Variance explained ---------------------------------------------------
pca_res$prop_expl_var$X * 100
pca_res$cum.var * 100


# 4. Loadings -------------------------------------------------------------
pca_res$loadings$X[, 1:2]


# 5. Scree plot -----------------------------------------------------------
pdf("out/scree_plot_mixomics.pdf", width = 3, height = 3)
plot(pca_res, main = "Scree Plot (mixOmics)")
dev.off()

# 6. Score plot (individuals)
pdf("out/score_plot_mixomics.png", width = 3, height = 3)
plotIndiv(
  pca_res,
  comp = c(1, 2),
  ind.names = TRUE,
  title = "PCA Score Plot — PC1 vs PC2",
  col = "steelblue",
  pch = 16,
  size.title = 1)
dev.off()


# 7. Loading plot (variables) ---------------------------------------------
pdf("out/loading_plot_mixomics.png", width = 3, height = 3)
plotLoadings(
  pca_res,
  comp = 1,
  title  = "Variable Loadings on PC1",
  size.title = 1)
dev.off()


# 8. Biplot (scores + loadings together) ----------------------------------
pdf("out/biplot_mixomics.png", width = 3, height = 3)
biplot(
  pca_res,
  comp = c(1, 2),
  ind.names = TRUE,
  var.names = TRUE,
  col = c("steelblue", "darkred"),
  cex = c(0.6, 0.9),
  main = "PCA Biplot (mixOmics)")
dev.off()


# 9. Correlation circle plot ----------------------------------------------
pdf("out/cor_circle_mixomics.png", width = 3, height = 3)
plotVar(
  pca_res,
  comp = c(1, 2),
  title = "Correlation Circle — PC1 vs PC2",
  col = "darkred",
  cex = 0.9)
dev.off()
