# AIM ---------------------------------------------------------------------
# sample PCA on Clinical Dataset

# lirbaries ---------------------------------------------------------------
library(tidyverse)
library(skimr)

# 1. Load data ------------------------------------------------------------
# data <- read.csv("data/clinical_data.csv", row.names = "patient_id")
data_full <- read.csv("data/dataset_michele_short.csv", row.names = "MS_ID")
data_full <- read.csv("data/dataset_short_all.csv", row.names = "MS_ID")

previous_PIRA
Disease_duration
future_PIRA
# trimm down only the var of interest
data <- data_full %>%
  filter(CTR0_MS1 == 1) %>%
  mutate(gender = case_when(Sex_F1_M2==2~1,T~0)) %>%
  dplyr::select(Age_at_prelievo,NfL,GFAP,CHIT,OPN) %>%
  filter(OPN > 0) %>%
  # mutate(NfL = log1p(NfL),
  #        GFAP = log1p(GFAP),
  #        OPN = log1p(OPN),
  #        CHIT = log1p(CHIT)) %>%
  mutate(NfL = NfL,
         GFAP = GFAP,
         OPN = OPN,
         CHIT = CHIT) %>%
  # mutate(NfL = log1p(NfL/Age_at_prelievo),
  #        GFAP = log1p(GFAP/Age_at_prelievo),
  #        OPN = log1p(OPN/Age_at_prelievo),
  #        CHIT = log1p(CHIT/Age_at_prelievo)) %>%
  # mutate(NfL = log1p(NfL/Age_at_prelievo),
  #        GFAP = log1p(GFAP/Age_at_prelievo),
  #        OPN = log1p(OPN),
  #        CHIT = log1p(CHIT/Age_at_prelievo)) %>%
  # dplyr::select(NfL,GFAP,CHIT,OPN,Age_at_prelievo)
  dplyr::select(NfL,GFAP,CHIT,OPN,Age_at_prelievo)

skim(data)

# -------------------------------------------------------------------------
data %>%
  mutate(NfL = log1p(NfL),
         GFAP = log1p(GFAP),
         OPN = log1p(OPN),
         CHIT = log1p(CHIT)) %>%
  ggpairs()

# 2. Scale and run PCA ----------------------------------------------------
pca_result <- prcomp(data, center = TRUE, scale. = TRUE)
print(pca_result$sdev)

# Variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cumvar <- cumsum(var_explained)

pct_table <- data.frame(
  PC = paste0("PC", seq_along(var_explained)),
  Variance_pct = round(var_explained * 100, 2),
  Cumulative_pct = round(cumvar * 100, 2)
)
print(pct_table)

# 3. Loadings (variable contributions) ------------------------------------
pca_result$rotation[,]

# 4. Scree plot -----------------------------------------------------------
scree_df <- data.frame(
  PC = factor(paste0("PC", 1:length(var_explained)), levels = paste0("PC", 1:length(var_explained))),
  Variance = var_explained * 100
)

ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_line(aes(group = 1), color = "darkred", linewidth = 0.8) +
  geom_point(color = "darkred", size = 2.5) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained (%)") +
  theme_bw()

# ggsave("out/plot/scree_plot_michele.pdf", p_scree, width = 6, height = 4)

# 5. Biplot (PC1 vs PC2) --------------------------------------------------
scores_df <- as.data.frame(pca_result$x)
scores_df$patient_id <- rownames(scores_df)

loadings_df <- as.data.frame(pca_result$rotation[, 1:2])
loadings_df$variable <- rownames(loadings_df)
scale_factor <- 3   # scale arrows for visibility

scores_df %>%
  left_join(data_full %>% rownames_to_column("patient_id"),by = "patient_id") %>%
  # mutate(PIRA)
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(col= factor(previous_PIRA)), size = 2.5, alpha = 0.8) +
  geom_text(aes(label = patient_id), size = 2.2, vjust = -0.6, color = "grey40") +
  geom_segment(
    data = loadings_df,
    aes(x = 0, y = 0, xend = PC1 * scale_factor, yend = PC2 * scale_factor),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "darkred", linewidth = 0.7
  ) +
  geom_text(
    data = loadings_df,
    aes(x = PC1 * scale_factor * 1.15, y = PC2 * scale_factor * 1.15, label = variable),
    color = "darkred", size = 3, fontface = "bold"
  ) +
  labs(
    title = "PCA Biplot — PC1 vs PC2",
    x = paste0("PC1 (", round(var_explained[1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2] * 100, 1), "% variance)")
  ) +
  theme_bw()

# ggsave("out/plot/pca_biplot_michele.pdf", p_biplot, width = 7, height = 6)

# 6. All pairwise PC combinations -----------------------------------------
scores_all <- as.data.frame(pca_result$x) %>%
  rownames_to_column("patient_id") %>%
  left_join(data_full %>% rownames_to_column("patient_id"),by = "patient_id")

ggpairs(
  scores_all,
  columns = c(2:5),
  mapping = aes(color = factor(Sex_F1_M2)), 
  upper = "blank",
  lower = list(
    continuous = wrap("smooth", alpha = 0.7, se = FALSE, linewidth = 0.5)
  ),
  diag = list(
    continuous = wrap("densityDiag", alpha = 0.4)
  )
) + 
  theme_bw() +
  # scale_color_manual(values = c("steelblue", "orange", "darkgreen")) +
  # scale_fill_manual(values = c("steelblue", "orange", "darkgreen"))
  theme(strip.background = element_blank())

ggsave("out/plot/pairs_all_pcs_gg.pdf", p_pairs_gg, width = 12, height = 12)
