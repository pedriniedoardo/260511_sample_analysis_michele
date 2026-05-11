# AIM -----------------------------------------------------------------------
# PCA + k-means clustering on neurological biomarkers (GFAP, NfL, CHIT1, OSP)
# WITHOUT age correction — PCA runs directly on raw marker expression.
# Age is never included as a PCA variable (it is a confounder, not a marker).
# Compare results with pca_kmeans_biomarkers.R to assess age influence on clusters.

# Libraries -----------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)  # BiocManager::install("ComplexHeatmap")
library(circlize)        # colorRamp2()
library(cluster)         # silhouette()

markers <- c("GFAP", "NfL", "CHIT1", "OSP")

# 1. Load data --------------------------------------------------------------
data <- read.csv("data/biomarkers.csv", row.names = "patient_id")
cat("Samples:", nrow(data), "| Variables:", ncol(data), "\n")

# 2. PCA on raw markers (no age correction) ---------------------------------
pca_res   <- prcomp(data[, markers], center = TRUE, scale. = TRUE)
var_exp   <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc_labels <- paste0("PC", 1:length(var_exp), " (", round(var_exp * 100, 1), "%)")

cat("\nVariance explained per PC:\n")
print(data.frame(PC = pc_labels, Cumulative = round(cumsum(var_exp) * 100, 1)))

# Scree plot
scree_df <- data.frame(
  PC       = factor(pc_labels, levels = pc_labels),
  Variance = var_exp * 100
)
p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_line(aes(group = 1), color = "darkred", linewidth = 0.8) +
  geom_point(color = "darkred", size = 2.5) +
  labs(title = "Scree Plot (raw markers, no age correction)",
       x = NULL, y = "Variance Explained (%)") +
  theme_bw()
ggsave("out/plot/scree_biomarkers_noCorr.pdf", p_scree, width = 5, height = 4)

scores <- as.data.frame(pca_res$x[, 1:2])

# 3. Optimal k: elbow + silhouette ------------------------------------------
k_range <- 2:7

wss <- sapply(k_range, function(k)
  kmeans(scores, centers = k, nstart = 50, iter.max = 100)$tot.withinss)

sil <- sapply(k_range, function(k) {
  km <- kmeans(scores, centers = k, nstart = 50, iter.max = 100)
  mean(silhouette(km$cluster, dist(scores))[, 3])
})

opt_df <- data.frame(k = k_range, WSS = wss, Silhouette = sil)

p_elbow <- ggplot(opt_df, aes(x = k, y = WSS)) +
  geom_line(color = "steelblue", linewidth = 0.9) +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = 3, linetype = "dashed", color = "darkred") +
  scale_x_continuous(breaks = k_range) +
  labs(title = "Elbow Plot", x = "Number of clusters (k)", y = "Total within-cluster SS") +
  theme_bw()

p_sil <- ggplot(opt_df, aes(x = k, y = Silhouette)) +
  geom_line(color = "steelblue", linewidth = 0.9) +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(xintercept = 3, linetype = "dashed", color = "darkred") +
  scale_x_continuous(breaks = k_range) +
  labs(title = "Silhouette Score", x = "Number of clusters (k)", y = "Mean silhouette width") +
  theme_bw()

ggsave("out/plot/elbow_plot_noCorr.pdf",      p_elbow, width = 5, height = 4)
ggsave("out/plot/silhouette_plot_noCorr.pdf", p_sil,   width = 5, height = 4)

# 4. k-means (k = 3) --------------------------------------------------------
set.seed(42)
km <- kmeans(scores, centers = 3, nstart = 100, iter.max = 200)

data$cluster <- factor(km$cluster)

# Re-label by increasing mean GFAP: Group 1 = lowest, Group 3 = highest
cluster_gfap_rank <- data %>%
  group_by(cluster) %>%
  summarise(mean_GFAP = mean(GFAP)) %>%
  arrange(mean_GFAP) %>%
  mutate(label = factor(paste("Group", 1:n())))

data <- data %>%
  left_join(cluster_gfap_rank[, c("cluster", "label")], by = "cluster") %>%
  mutate(cluster = label) %>%
  select(-label)

cat("\nCluster sizes:\n")
print(table(data$cluster))

cat("\nCluster means (original expression):\n")
print(data %>% group_by(cluster) %>% summarise(across(all_of(markers), mean)))

cat("\nMean age per cluster (check for age-driven separation):\n")
print(data %>% group_by(cluster) %>% summarise(mean_age = mean(age), sd_age = sd(age)))

# 5. PCA score plot — color by cluster, size by age ------------------------
# Age encoded as point size reveals if cluster separation tracks age
scores$cluster    <- data$cluster
scores$patient_id <- rownames(data)
scores$age        <- data$age

p_pca <- ggplot(scores, aes(x = PC1, y = PC2, color = cluster, size = age)) +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = patient_id), size = 1.9, vjust = -0.8,
            color = "grey40", show.legend = FALSE) +
  stat_ellipse(aes(group = cluster), linewidth = 0.7, linetype = "dashed",
               show.legend = FALSE) +
  scale_color_manual(values = c("Group 1" = "#2196F3",
                                "Group 2" = "#FF9800",
                                "Group 3" = "#F44336")) +
  scale_size_continuous(range = c(1.5, 5), breaks = c(40, 55, 70)) +
  labs(
    title    = "PCA Score Plot — k-means clusters (no age correction)",
    subtitle = "Point size = age: inspect whether clusters align with age gradient",
    x        = pc_labels[1],
    y        = pc_labels[2],
    color    = "Cluster",
    size     = "Age (years)"
  ) +
  theme_bw()

ggsave("out/plot/pca_kmeans_scores_noCorr.pdf", p_pca, width = 7, height = 6)

# 6. Heatmap (ComplexHeatmap) -----------------------------------------------
mat   <- data %>% select(all_of(markers)) %>% as.matrix()
mat_z <- t(scale(mat))
colnames(mat_z) <- rownames(data)

col_order <- order(data$cluster)
mat_z     <- mat_z[, col_order]
data_ord  <- data[col_order, ]

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

age_col_fun <- colorRamp2(
  c(min(data_ord$age), max(data_ord$age)),
  c("#FFF9C4", "#E65100")
)

top_anno <- HeatmapAnnotation(
  Cluster = data_ord$cluster,
  Age     = data_ord$age,
  col     = list(
    Cluster = c("Group 1" = "#2196F3", "Group 2" = "#FF9800", "Group 3" = "#F44336"),
    Age     = age_col_fun
  ),
  annotation_label  = c("Cluster", "Age (years)"),
  annotation_height = unit(c(4, 3), "mm"),
  show_legend       = TRUE
)

sil_scores <- silhouette(as.integer(data_ord$cluster), dist(scores[col_order, 1:2]))[, 3]

right_anno <- rowAnnotation(
  Silhouette = anno_barplot(
    sil_scores,
    gp         = gpar(fill = ifelse(sil_scores > 0, "steelblue", "tomato")),
    axis_param = list(gp = gpar(fontsize = 7)),
    width      = unit(2, "cm")
  )
)

ht <- Heatmap(
  mat_z,
  name             = "z-score",
  col              = col_fun,
  top_annotation   = top_anno,
  right_annotation = right_anno,
  cluster_rows     = TRUE,
  cluster_columns  = FALSE,
  show_column_names = TRUE,
  column_names_gp  = gpar(fontsize = 6),
  row_names_gp     = gpar(fontsize = 10, fontface = "bold"),
  row_dend_side    = "left",
  column_title     = "Biomarker Expression by Cluster (z-scored)\nRaw PCA + k-means (k = 3) — no age correction",
  column_title_gp  = gpar(fontsize = 11),
  column_split     = data_ord$cluster,
  column_gap       = unit(2, "mm"),
  border           = TRUE,
  heatmap_legend_param = list(title = "z-score", direction = "horizontal")
)

pdf("out/plot/heatmap_clusters_noCorr.pdf", width = 9, height = 5)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
cat("Heatmap saved to out/plot/heatmap_clusters_noCorr.pdf\n")
