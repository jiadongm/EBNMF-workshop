# 01_PCA.R - PCA baseline for COVID-19 proteomics
#
# Input:  results/X_prot.rds        240 × 481 matrix, NPX log2 (samples × proteins)
#         results/sample_meta.rds   sample-level clinical metadata
# Output: results/pca_prot.rds      prcomp object
#         results/pca_scores.rds    data frame: PC scores + metadata
#         results/figs/prot_pca_*   figures
#
# Scaling: z-score per protein (center=TRUE, scale.=TRUE), i.e. subtract per-protein
# mean and divide by per-protein SD before PCA.  Each protein contributes equally
# regardless of dynamic range.
#
# Multi-PC visualisation: vizOmics::matrixPlot() - pairwise scatter + diagonal density
# for PC1–5, coloured by one metadata variable per figure.

library(ggplot2)
library(patchwork)
library(dplyr)
library(vizOmics)

# ── Paths ─────────────────────────────────────────────────────────────────────

results_dir <- "results"
fig_dir     <- file.path(results_dir, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load data ──────────────────────────────────────────────────────────────

X_prot <- readRDS(file.path(results_dir, "X_prot.rds"))   # 240 × 481
meta   <- readRDS(file.path(results_dir, "sample_meta.rds"))

# Ensure row order of meta matches X_prot
meta <- meta[match(rownames(X_prot), meta$sample), ]
stopifnot(all(meta$sample == rownames(X_prot)))

cat(sprintf("Data: %d samples × %d proteins\n", nrow(X_prot), ncol(X_prot)))

# ── 2. Compute WHO delta (T2 - T1) per patient ────────────────────────────────
# Positive = worsening; negative = improvement; 0 = stable.
# This is a patient-level variable - attach to both T1 and T2 rows.

who_wide <- meta |>
  filter(!who_ambiguous) |>                      # exclude ambiguous "1 or 2" scores
  select(patient_id, time_point, who_score_num) |>
  tidyr::pivot_wider(names_from = time_point, values_from = who_score_num) |>
  mutate(who_delta = T2 - T1)                    # positive = got worse

meta <- meta |>
  left_join(who_wide |> select(patient_id, who_delta), by = "patient_id")

cat(sprintf("WHO delta computed for %d patients (excluded %d with ambiguous scores)\n",
            sum(!is.na(who_wide$who_delta)), as.integer(sum(meta$who_ambiguous) / 2)))

# ── 3. PCA ────────────────────────────────────────────────────────────────────

cat("Running PCA (center=TRUE, scale.=TRUE)...\n")
pca <- prcomp(X_prot, center = TRUE, scale. = TRUE)

pve      <- pca$sdev^2 / sum(pca$sdev^2)
cum_pve  <- cumsum(pve)

cat(sprintf("  PC1: %.1f%%,  PC2: %.1f%%,  PC3: %.1f%%\n",
            100*pve[1], 100*pve[2], 100*pve[3]))
cat(sprintf("  PCs to reach 50%% variance: %d\n", which(cum_pve >= 0.50)[1]))
cat(sprintf("  PCs to reach 80%% variance: %d\n", which(cum_pve >= 0.80)[1]))

# Score data frame: PC1–10 + metadata
scores_df <- as.data.frame(pca$x[, 1:10])
scores_df <- cbind(meta, scores_df)

# ── 4. Scree plot ─────────────────────────────────────────────────────────────

scree_df <- data.frame(
  PC      = 1:30,
  pve     = pve[1:30] * 100,
  cum_pve = cum_pve[1:30] * 100
)

p_scree <- ggplot(scree_df, aes(x = PC)) +
  geom_col(aes(y = pve), fill = "steelblue", alpha = 0.8, width = 0.7) +
  geom_line(aes(y = cum_pve / 4), colour = "firebrick", linewidth = 0.8) +
  geom_point(aes(y = cum_pve / 4), colour = "firebrick", size = 1.5) +
  scale_y_continuous(
    name     = "% variance explained (bar)",
    sec.axis = sec_axis(~ . * 4, name = "Cumulative % variance (line)")
  ) +
  scale_x_continuous(breaks = seq(2, 30, 2)) +
  labs(title = "PCA scree plot - proteomics (481 proteins, 240 samples)",
       x = "Principal component") +
  theme_bw()

ggsave(file.path(fig_dir, "prot_pca_scree.pdf"), p_scree, width = 8, height = 4)
cat("Saved: prot_pca_scree.pdf\n")

# ── 5. matrixPlot panels coloured by different metadata ───────────────────────
# Each call plots PC1–5 as a pairwise matrix (5×5 grid with diagonal densities).

N_COMP <- 5   # number of PCs to show in the matrix plot
scores_mat <- pca$x[, 1:N_COMP]   # 240 × 5 matrix passed to matrixPlot()

# Helper: save a matrixPlot to PDF.
# matrixPlot() internally calls gridExtra::grid.arrange(), which draws immediately
# to the current graphics device. The call must therefore happen INSIDE pdf()/dev.off().
save_matplot <- function(colBy, title, fname, width = 12, height = 12,
                         manualCol = NULL, legendTitle = "") {
  pdf(file.path(fig_dir, fname), width = width, height = height)
  matrixPlot(
    scores      = scores_mat,
    max_ncomp   = N_COMP,
    colBy       = colBy,
    pointSize   = 1.5,
    pointAlpha  = 0.6,
    manualCol   = manualCol,
    legendTitle = legendTitle,
    compName    = "PC"
  )
  dev.off()
  cat(sprintf("Saved: %s\n", fname))
}

# 5a. Blood draw time point (T1 = early ~days 1-4, T2 = later ~days 7-10)
save_matplot(
  colBy       = meta$time_point,
  title       = "PCA coloured by blood draw time point (T1 / T2)",
  fname       = "prot_pca_timepoint.pdf",
  manualCol   = c(T1 = "steelblue", T2 = "firebrick"),
  legendTitle = "Time point"
)

# 5b. WHO severity score (continuous, 1–7)
save_matplot(
  colBy       = meta$who_score_num,
  title       = "PCA coloured by WHO ordinal scale (1=mild, 7=critical)",
  fname       = "prot_pca_severity.pdf",
  legendTitle = "WHO score"
)

# 5c. Age
save_matplot(
  colBy       = meta$age,
  title       = "PCA coloured by patient age",
  fname       = "prot_pca_age.pdf",
  legendTitle = "Age (years)"
)

# 5d. Sex
save_matplot(
  colBy       = meta$sex,
  title       = "PCA coloured by sex",
  fname       = "prot_pca_sex.pdf",
  manualCol   = c(Female = "#E07B8A", Male = "#5B8DB8"),
  legendTitle = "Sex"
)

# 5e. Patient location (Home / Clinic / Hospital / ICU) - proxy for severity tier
save_matplot(
  colBy       = meta$location,
  title       = "PCA coloured by patient location (Home / Clinic / Hospital / ICU)",
  fname       = "prot_pca_location.pdf",
  legendTitle = "Location"
)

# 5f. BMI
save_matplot(
  colBy       = meta$bmi,
  title       = "PCA coloured by BMI",
  fname       = "prot_pca_bmi.pdf",
  legendTitle = "BMI"
)

# 5g. WHO delta T2 - T1 (improvement vs worsening)
save_matplot(
  colBy       = scores_df$who_delta,
  title       = "PCA coloured by WHO score change T2-T1 (negative=improvement)",
  fname       = "prot_pca_who_delta.pdf",
  legendTitle = "dWHO"
)

# ── 6. Trajectory plot: PC1 vs PC2, arrows T1 -> T2 per patient ───────────────
# Arrow base = T1 position; arrowhead = T2 position.
# Coloured by the patient's T1 WHO severity score.

traj <- scores_df |>
  select(patient_id, time_point, PC1, PC2, who_score_num) |>
  tidyr::pivot_wider(
    names_from  = time_point,
    values_from = c(PC1, PC2, who_score_num)
  ) |>
  rename(
    x0        = PC1_T1,  y0 = PC2_T1,
    x1        = PC1_T2,  y1 = PC2_T2,
    severity  = who_score_num_T1
  )

p_traj <- ggplot(traj, aes(x = x0, y = y0, xend = x1, yend = y1, colour = severity)) +
  geom_segment(
    arrow     = arrow(length = unit(0.15, "cm"), type = "closed"),
    linewidth = 0.5,
    alpha     = 0.75
  ) +
  scale_colour_gradient2(
    low      = "steelblue",
    mid      = "gold",
    high     = "firebrick",
    midpoint = 4,
    name     = "WHO score\n(T1 baseline)"
  ) +
  labs(
    title    = "Patient trajectories in PC1/PC2 space (T1 -> T2)",
    subtitle = sprintf("PC1: %.1f%% var,  PC2: %.1f%% var",
                       100*pve[1], 100*pve[2]),
    x = sprintf("PC1 (%.1f%%)", 100*pve[1]),
    y = sprintf("PC2 (%.1f%%)", 100*pve[2])
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "prot_pca_trajectories.pdf"),
       p_traj, width = 7, height = 5.5)
cat("Saved: prot_pca_trajectories.pdf\n")

# ── 7. PC loadings: top proteins driving PC1 and PC2 ─────────────────────────

top_n <- 15

make_loading_plot <- function(pc_idx) {
  ld   <- pca$rotation[, pc_idx]
  df   <- data.frame(protein = names(ld), loading = ld) |>
    arrange(desc(abs(loading))) |>
    slice_head(n = top_n) |>
    mutate(
      protein = factor(protein, levels = rev(protein)),
      sign    = ifelse(loading > 0, "positive", "negative")
    )
  ggplot(df, aes(x = loading, y = protein, fill = sign)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(positive = "steelblue", negative = "firebrick"),
                      guide  = "none") +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(title = sprintf("Top %d protein loadings - PC%d", top_n, pc_idx),
         x = "Loading", y = NULL) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))
}

p_load <- make_loading_plot(1) | make_loading_plot(2)

ggsave(file.path(fig_dir, "prot_pca_loadings.pdf"),
       p_load, width = 10, height = 5)
cat("Saved: prot_pca_loadings.pdf\n")

# ── 8. Save objects ───────────────────────────────────────────────────────────

saveRDS(pca,       file.path(results_dir, "pca_prot.rds"))
saveRDS(scores_df, file.path(results_dir, "pca_scores.rds"))

cat("\nSaved:\n")
cat("  results/pca_prot.rds    - prcomp object\n")
cat("  results/pca_scores.rds  - PC scores + metadata data frame\n")
