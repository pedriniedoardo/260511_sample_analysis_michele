# AIM -----------------------------------------------------------------------
# PCA + k-means clustering on neurological biomarkers (GFAP, NfL, CHIT1, OSP)
# Age is a known confounder: we regress it out before PCA so that clusters reflect disease biology, not age-driven expression differences.

# Libraries -----------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(viridis)
library(scales)
library(GGally)
library(ggrepel)

# specify the markers
markers <- c("GFAP", "NfL", "CHIT", "OPN")

# 1. Load data --------------------------------------------------------------
# load the new dataset fixed by Michele.
# the 0 have been converted to LOD/2
# addition of specific covariates in the datastet

data <- read.csv("data/data_new_full.csv",row.names = "sample_id") %>%
  # this is already only on MS patient
  # filter(pt == 1) %>%
  dplyr::select(age, NfL, GFAP, CHIT,OPN) %>%
  # this is already filtered for missing values of OPN
  # filter(!is.na(OPN)) %>%
  # now we do not have zeros anymore, therefore use the log instead of the log1p
  mutate(NfL = log(NfL),
         GFAP = log(GFAP),
         CHIT = log(CHIT),
         OPN = log(OPN))

# confirm the dimension of the dataet
data %>%  
  dim()

# 2. Diagnostic: age-marker correlations ------------------------------------
# Confirm age confounding before correcting for it
cor_df <- data %>%
  pivot_longer(all_of(markers), names_to = "marker", values_to = "value")

ggplot(cor_df, aes(x = age, y = value)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  facet_wrap(~marker, scales = "free_y") +
  labs(title = "Age vs. Marker concentration",
       x = "Age (years)", y = "log1 concentration") +
  theme_bw() +
  theme(strip.background = element_blank())

# ggsave("out/plot/age_marker_correlation_UPDATE.pdf", width = 8, height = 8)

# add also a variable to variable correlation
# lower: scatter points + lm trend with SE band
# upper: Pearson r with significance stars
# diag:  density curve
ggpairs(
  data,
  lower = list(
    continuous = function(data, mapping, ...) {
      ggplot(data = data, mapping = mapping) +
        geom_point(alpha = 0.45, size = 1.2) +
        geom_smooth(method = "lm", se = TRUE,
                    linewidth = 0.7, alpha = 0.1) +
        theme_bw()
    }
  ),
  upper = list(
    continuous = wrap("cor", method = "pearson", size = 3.2)
  ),
  diag = list(
    continuous = wrap("densityDiag", alpha = 0.4, color = NA,fill = "gray25")
  )
) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 8))
# ggsave("out/plot/cross_correlation_UPDATE.pdf", width = 8, height = 8)

# 2b. Model fit + residuals per point ---------------------------------------
# For each marker: fit lm(marker ~ age), collect observed, fitted, and residual.
# Vertical segments connect each observed point to its model estimate, making explicit what the regression "explains" vs what remains (the residual).
resid_df <- map_dfr(markers, function(m) {
  fit <- lm(reformulate("age", response = m), data = data)
  data.frame(
    age = data$age,
    observed = data[[m]],
    fitted = fitted(fit),
    residual = residuals(fit),
    marker = m
  )
})

ggplot(resid_df, aes(x = age)) +
  # residual segments: observed fitted value on the model line
  geom_segment(
    aes(xend = age, y = observed, yend = fitted, color = residual > 0),
    linewidth = 0.55, alpha = 0.85
  ) +
  # model regression line
  geom_line(aes(y = fitted), linewidth = 1) +
  # observed points
  geom_point(aes(y = observed),col="gray25", size = 1.8, alpha = 0.85) +
  # fitted values on the line (crosses)
  geom_point(aes(y = fitted), shape = 4,col = "gray25", size = 1.6, stroke = 0.7,alpha=0.5) +
  scale_color_manual(
    values = c("TRUE"  = "red",   # above the line
               "FALSE" = "blue"),  # below the line
    labels = c("TRUE"  = "Positive",
               "FALSE" = "Negative"),
    name   = "Residual direction"
  ) +
  facet_wrap(~marker, scales = "free_y", nrow = 2) +
  labs(
    title    = "Linear Model Fit and Residuals per Sample",
    # subtitle = "Dot = observed | Cross = model estimate | Segment = residual kept after correction",
    x        = "Age (years)",
    y        = "log1p concentration"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",strip.background = element_blank())
# ggsave("out/plot/age_model_residuals_UPDATE.pdf", width = 8, height = 8)

# 3. Age correction ---------------------------------------------------------
# For each marker: fit lm(marker ~ age), replace with residuals.
# Residuals = marker expression with the linear age effect removed.
# This ensures PCA and clustering capture disease signal, not age variation.
corrected <- data
for (m in markers) {
  fit <- lm(reformulate("age", response = m), data = data)
  corrected[[m]] <- residuals(fit)
}
corrected_markers <- corrected[, markers]

# save the table to share with michele
# corrected_markers %>%
#   rownames_to_column("sample_id") %>%
#   write_tsv("out/table/corrected_markers_UPDATE.tsv")

# 4. PCA on age-corrected markers -------------------------------------------
pca_res <- prcomp(corrected_markers, center = TRUE, scale. = TRUE)
var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc_labels <- paste0("PC", 1:length(var_exp), " (", round(var_exp * 100, 1), "%)")

print(data.frame(PC = pc_labels, Cumulative = round(cumsum(var_exp) * 100, 1)))

# Scree plot
scree_df <- data.frame(
  PC = factor(pc_labels, levels = pc_labels),
  Variance = var_exp * 100
)
ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_col(alpha = 0.8) +
  geom_line(aes(group = 1), linewidth = 0.8) +
  geom_point(size = 2.5) +
  labs(title = "Scree Plot (age-corrected markers)",
       x = NULL, y = "Variance Explained (%)") +
  theme_bw()
# ggsave("out/plot/scree_biomarkers_UPDATE.pdf", width = 5, height = 4)

# Use PC1 + PC2 for clustering (typically captures most variance with 4 vars)
scores <- as.data.frame(pca_res$x[, 1:2])

# 5. Optimal k: elbow + silhouette ------------------------------------------
# k_range <- 2:7
# 
# wss <- sapply(k_range, function(k)
#   kmeans(scores, centers = k, nstart = 50, iter.max = 100)$tot.withinss)
# 
# sil <- sapply(k_range, function(k) {
#   km  <- kmeans(scores, centers = k, nstart = 50, iter.max = 100)
#   mean(silhouette(km$cluster, dist(scores))[, 3])
# })
# 
# opt_df <- data.frame(k = k_range, WSS = wss, Silhouette = sil)
# 
# ggplot(opt_df, aes(x = k, y = WSS)) +
#   geom_line(color = "steelblue", linewidth = 0.9) +
#   geom_point(size = 3, color = "steelblue") +
#   geom_vline(xintercept = 3, linetype = "dashed", color = "darkred") +
#   scale_x_continuous(breaks = k_range) +
#   labs(title = "Elbow Plot", x = "Number of clusters (k)", y = "Total within-cluster SS") +
#   theme_bw()
# 
# ggplot(opt_df, aes(x = k, y = Silhouette)) +
#   geom_line(color = "steelblue", linewidth = 0.9) +
#   geom_point(size = 3, color = "steelblue") +
#   geom_vline(xintercept = 3, linetype = "dashed", color = "darkred") +
#   scale_x_continuous(breaks = k_range) +
#   labs(title = "Silhouette Score", x = "Number of clusters (k)", y = "Mean silhouette width") +
#   theme_bw()

# ggsave("out/plot/elbow_plot.pdf",     p_elbow, width = 5, height = 4)
# ggsave("out/plot/silhouette_plot.pdf", p_sil,   width = 5, height = 4)

# martina suggested to force it to 3 groups

# 6. k-means (k = 4) --------------------------------------------------------
set.seed(42)
km <- kmeans(scores, centers = 4, nstart = 100, iter.max = 200)

# Attach cluster labels to full data
data$cluster <- factor(km$cluster)

# Re-label clusters by increasing mean GFAP so Group 1 = lowest, 3 = highest
cluster_gfap_rank <- data %>%
  group_by(cluster) %>%
  summarise(mean_GFAP = mean(GFAP)) %>%
  arrange(mean_GFAP) %>%
  mutate(label = factor(paste("Group", 1:n())))

data <- data %>%
  rownames_to_column("sample_id") %>%
  left_join(cluster_gfap_rank[, c("cluster", "label")], by = "cluster") %>%
  mutate(cluster = label) %>%
  dplyr::select(-label)

print(table(data$cluster))
print(data %>% group_by(cluster) %>% summarise(across(all_of(markers), mean)))

# 7. Biplot: scores + loading vectors --------------------------------------
scores$cluster <- data$cluster
scores$sample_id <- data$sample_id

# Scale loadings into the score space so arrows are readable.
# Factor = 80% of the max absolute score range / max absolute loading.
scale_factor <- 0.8 * max(abs(scores[, 1:2])) / max(abs(pca_res$rotation[, 1:2]))

loadings_df <- as.data.frame(pca_res$rotation[, 1:2]) %>%
  mutate(
    variable = rownames(pca_res$rotation),
    x_end = PC1 * scale_factor,
    y_end = PC2 * scale_factor
  )

ggplot(scores, aes(x = PC1, y = PC2, color = cluster)) +
  geom_text_repel(aes(label = sample_id), size = 1.9, vjust = -0.7,
                  color = "grey40", show.legend = FALSE,bg.color = "white",max.overlaps = 3,
                  point.padding = unit(0.1, 'lines'), 
                  min.segment.length = unit(0.5, 'lines'),segment.alpha=0.5) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(linewidth = 0.7, linetype = "dashed") +
  # loading arrows
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = x_end, yend = y_end),
    arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
    color = "gray25", linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  # variable labels slightly beyond arrow tip
  geom_text_repel(data = loadings_df,
                  aes(x = x_end * 1.12, y = y_end * 1.12, label = variable),
                  color = "gray25", size = 3.2, fontface = "bold",
                  inherit.aes = FALSE) +
  scale_color_manual(values = c("Group 1" = "#2196F3",
                                "Group 2" = "#FF9800",
                                "Group 3" = "#F44336",
                                "Group 4" = "#4CAF50"))+
  labs(title = "PCA Biplot k-means clusters (age-corrected)",
       x = pc_labels[1],
       y = pc_labels[2],
       color = "Cluster"
  ) +
  theme_bw()
ggsave("out/plot/pca_biplot_kmeans_UPDATE_km4.pdf", width = 7, height = 6)

# How to read the biplot:                                                                              
# - Arrow direction — which way a marker pulls samples in PC space
# - Arrow length — how much that marker contributes to the displayed PCs (longer = stronger loading)     
# - Arrows pointing the same way — those markers co-vary (e.g. GFAP and NfL likely point similarly) 
# - Sample position relative to an arrow — samples in the direction of an arrow have high expression of that marker

# 8. Heatmap (ComplexHeatmap) -----------------------------------------------
# Show ORIGINAL (not residual) z-scored expression
# Annotate rows by cluster and age group — confirms clusters are not age-driven.
mat <- data %>%
  dplyr::select(all_of(markers)) %>%
  as.matrix()

# z-score per marker; ComplexHeatmap expects markers as rows
mat_z <- t(scale(mat))    
colnames(mat_z) <- data$sample_id

# Compute silhouette on original row order before any reordering.
# scores[, 1:2] drops the cluster/sample_id columns added in step 7.
sil_obj <- silhouette(as.integer(data$cluster), dist(scores[, 1:2]))
data$sil_width <- sil_obj[, 3]

# Sort: cluster first (keeps column_split blocks intact),
# then descending silhouette so the most confidently assigned samples
# appear at the left of each block and borderline cases at the right.
col_order <- order(data$cluster, -data$sil_width)
mat_z <- mat_z[, col_order]
data_ord  <- data[col_order, ]

# Color scale: blue → white → red
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

# Age as a continuous bar (gradient)
# age_col_fun <- colorRamp2(c(min(data_ord$age),max(data_ord$age)),c("#FFF9C4", "#E65100"))
show_col(viridis(option = "turbo",20))
age_col_fun <- colorRamp2(seq(from = min(data_ord$age), to = max(data_ord$age),length=20),viridis(option = "turbo",20))

# Top annotations: cluster (block) + age (gradient bar)
top_anno <- HeatmapAnnotation(
  Cluster = data_ord$cluster,
  Age = data_ord$age,
  col = list(
    Cluster = c("Group 1" = "#2196F3",
                "Group 2" = "#FF9800",
                "Group 3" = "#F44336",
                "Group 4" = "#4CAF50"),
    Age = age_col_fun
  ),
  annotation_label  = c("Cluster", "Age (years)"),
  annotation_height = unit(c(4, 3), "mm"),
  show_legend       = TRUE
)

# Bottom annotation: silhouette width per sample (quality of cluster assignment)
# Samples are columns → silhouette bar must be a column (bottom) annotation, not a row annotation. right_annotation matches nrow(mat_z) = 4 markers, not samples.
# sil_width already computed and attached to data above; just reorder.
sil_scores <- data_ord$sil_width

bottom_anno <- HeatmapAnnotation(
  Silhouette = anno_barplot(
    sil_scores,
    gp = gpar(fill = ifelse(sil_scores > 0, "steelblue", "tomato")),
    axis_param = list(gp = gpar(fontsize = 7)),
    height = unit(2, "cm")
  ),
  annotation_label = "Silhouette width"
)

ht <- Heatmap(
  mat_z,
  name = "z-score",
  col = col_fun,
  top_annotation = top_anno,
  bottom_annotation = bottom_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names  = TRUE,
  column_names_gp = gpar(fontsize = 3),
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  row_dend_side = "left",
  column_title = "Biomarker Expression by Cluster (z-scored)\nAge-corrected PCA + k-means (k = 4)",
  column_title_gp = gpar(fontsize = 11),
  column_split = data_ord$cluster,
  column_gap = unit(2, "mm"),
  border = TRUE,
  heatmap_legend_param = list(title = "z-score", direction = "horizontal")
)

pdf("out/plot/heatmap_clusters_UPDATE_km4.pdf", width = 12, height = 6)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# How to interpret the silhouette plot                                                                   
# The silhouette width for a single sample measures how well it fits its assigned cluster versus the nearest alternative cluster:                                                                           
# s(i) = (b - a) / max(a, b)                                                                             
# a = mean distance to all other samples in the SAME cluster                                             
# b = mean distance to all samples in the NEAREST other cluster                                        

# The value ranges from -1 to +1:                                                                        
# Close to +1 │ Sample sits deep inside its cluster, far from neighbours - well assigned
# Close to 0  │ Sample lies near the boundary between two clusters - ambiguous
# Negative │ Sample is closer to another cluster than its own - likely misclassified

# In the heatmap bottom bar specifically:                                                                
# - A blue bar (positive) = that sample is confidently placed in its group
# - A red bar (negative) = that sample may belong to a different group; worth inspecting as a borderline case                                                                                                  
# - The mean silhouette across all samples is the overall quality score for your chosen k - the elbow/silhouette plot shows this mean for each k, and you pick the k where it peaks                  

# A mean silhouette > 0.5 is generally considered a reasonable separation; > 0.7 is strong.

# 9. Covariate balance across clusters --------------------------------------
# read in the full dataset
data_full <- read_csv("data/data_new_full.csv") %>%
  # filter(CTR0_MS1 == 1) %>%
  left_join(data %>%
              dplyr::select(sample_id, cluster), by = "sample_id")

# cat_vars  <- c("Sex_F1_M2", "future_PIRA", "previous_PIRA",
#                "RELAPSE_act", "MRI_act", "Inflammatory_act")

cat_vars  <- c("future_PIRA",
               "previous_PIRA",
               "RELAPSE_act",
               "MRI_act",
               "Inflammatory_act",
               "pira_group",
               "phenot")

cluster_colors <- c("Group 1" = "#2196F3",
                    "Group 2" = "#FF9800",
                    "Group 3" = "#F44336",
                    "Group 4" = "#4CAF50")

# Recode all categorical variables to labelled factors
data_full <- data_full %>%
  mutate(
    phenot = factor(phenot, levels = c(0, 1),labels = c("0", "1")),
    pira_group = factor(pira_group, levels = c(0, 1,2),labels = c("0", "1","2")),
    
    Sex_F1_M2 = factor(Sex_F1_M2, levels = c(1, 2),labels = c("Female", "Male")),
    future_PIRA = factor(future_PIRA, levels = c(0, 1),labels = c("0", "1")),
    previous_PIRA = factor(previous_PIRA, levels = c(0, 1),labels = c("0", "1")),
    RELAPSE_act = factor(RELAPSE_act, levels = c(0, 1),labels = c("0", "1")),
    MRI_act = factor(MRI_act, levels = c(0, 1),labels = c("0", "1")),
    Inflammatory_act = factor(Inflammatory_act, levels = c(0, 1), labels = c("0", "1")))

# 9a. Categorical variables: proportional stacked bars + Fisher exact -------
cat_props <- data_full %>%
  filter(!is.na(cluster)) %>%
  dplyr::select(cluster, all_of(cat_vars)) %>%
  pivot_longer(-cluster, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(value = factor(value)) %>%
  group_by(variable, cluster, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(variable, cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

cat_props %>%
  arrange(n)

# Fisher exact test per variable (more robust than chi-sq for small n)
set.seed(123)
fisher_pvals <- map_dfr(cat_vars, function(v) {
  df <- data_full %>% filter(!is.na(cluster), !is.na(.data[[v]]))
  tbl <- table(df$cluster, df[[v]])
  p <- tryCatch(
    fisher.test(tbl, simulate.p.value = TRUE, B = 10000)$p.value,
    error = function(e) NA_real_
  )
  tibble(
    variable = v,
    label    = case_when(
      is.na(p) ~ "n.s.",
      p < 0.001 ~ "p < 0.001",
      p < 0.05 ~ paste0("p = ", round(p, 3)),
      TRUE ~ paste0("p = ", round(p, 2))
    )
  )
})

ggplot(cat_props, aes(x = cluster, y = prop, fill = value)) +
  geom_col(position = "stack", width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(prop > 0.08, scales::percent(prop, accuracy = 1), "")),
    position  = position_stack(vjust = 0.5),
    size = 2.8, color = "white", fontface = "bold"
  ) +
  geom_text(
    data        = fisher_pvals,
    aes(x = 2, y = 1.07, label = label),
    inherit.aes = FALSE, size = 2.6, color = "grey30"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.13)) +
  # scale_fill_brewer(palette = "Set2") +
  facet_wrap(~variable, nrow = 2) +
  labs(
    title    = "Categorical covariate distribution across clusters",
    subtitle = "Fisher exact test",
    x = "Cluster", y = "Proportion", fill = NULL
  ) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text       = element_text(size = 9),
        legend.position  = "bottom")

ggsave("out/plot/covariate_categorical_UPDATE_km4_stacked.pdf", width = 15, height = 7)

ggplot(cat_props, aes(x = cluster, y = prop, fill = value)) +
  geom_col(position = "dodge", width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = scales::percent(prop, accuracy = 1)),
    position = position_dodge(width = 0.65),
    vjust = -1, # Positive values move text down inside the bar; negative values move it above
    size = 2, fontface = "bold"
  ) +
  geom_text(
    data        = fisher_pvals,
    aes(x = 2, y = 1.07, label = label),
    inherit.aes = FALSE, size = 2.6, color = "grey30"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.13)) +
  # scale_fill_brewer(palette = "Set2") +
  facet_wrap(~variable, nrow = 2) +
  labs(
    title    = "Categorical covariate distribution across clusters",
    subtitle = "Fisher exact test",
    x = "Cluster", y = "Proportion", fill = NULL
  ) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text       = element_text(size = 9),
        legend.position  = "bottom")

ggsave("out/plot/covariate_categorical_UPDATE_km4_dodge.pdf", width = 15, height = 8)

# A significant result tells you only that at least one group differs — it does not say which one. This is the same logic as ANOVA: significant F-test → something is different, but you need post-hoc comparisons to find where.                                                                             
