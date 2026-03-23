# 07_delta_factor_plots.R — Flashier-vignette-style factor visualisations for delta analysis
#
# Mirrors 05_factor_plots.R but for the delta (T2-T1) matrix.
# Five variants: D1 (NMF), D2 (EBMF Laplace), D3 (EBMF Normal), D4 (GBCD), D5 (Asym EBMF).
#
# A. F-matrix pheatmaps (protein × factor) for all five NMF variants
#    A1  D1 standard NMF  (481×6,  non-negative, col-normalised [0,1])
#    A2  D2 EBMF Laplace  (481×20, signed, diverging blue→white→red, PVE order)
#    A3  D3 EBMF Normal   (481×20, signed, diverging, PVE order)
#    A4  D4 GBCD          (481×18, signed LFC, diverging)
#    A5  D5 Asym EBMF     (481×20, signed F, diverging, PVE order)
#
# B. flashier single-cell vignette style (fit_D2 — signed L, limited plot types)
#    B1  Scree / PVE bar chart
#    B2  Histograms × 2 colourings (who_delta_cat / T1_severity)
#    B3  SKIPPED — structure requires non-negative L
#    B4  SKIPPED — L heatmap requires non-negative L (see 06 output)
#    B5  Top-protein F heatmap (D2 factors 17, 5, 2, 4)
#    B6  F scatter plots for D2 factors 17, 5, 2
#
# C. flashier_intro bar plots for D2 (bar loadings works even with signed L)
#
# D. flashier plots for D5 (non-negative L → ALL plot types enabled)
#    D1  Scree / PVE bar chart
#    D2  Histograms × 2 colourings
#    D3  Structure plots × 2 colourings  ← enabled by non-neg L
#    D4  L heatmaps × 2 colourings       ← enabled by non-neg L
#    D5  Top-protein F heatmap (D5 most-associated factors)
#    D6  F scatter plots
#    D7  Bar plots × 2 colourings
#
# Patient colourings:
#   who_delta_cat : Improved (<0) / Stable (=0) / Worsened (>0) / Ambiguous (NA)
#   T1_severity   : WHO score at admission, levels 1, 1.5, 3, 4, 5, 6, 7
#
# Run from CovidCase/:   Rscript 07_delta_factor_plots.R

suppressPackageStartupMessages({
  library(NMF)
  library(flashier)
  library(ebnm)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(grid)
})

# ── Paths ──────────────────────────────────────────────────────────────────────

results_dir <- "results"
nmf_dir     <- file.path(results_dir, "nmf")
fig_dir     <- file.path(results_dir, "figs", "nmf")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load data ───────────────────────────────────────────────────────────────

X_prot     <- readRDS(file.path(results_dir, "X_prot.rds"))       # 240 × 481
pca_scores <- readRDS(file.path(results_dir, "pca_scores.rds"))
assoc      <- read.csv(file.path(nmf_dir, "delta_association.csv"))

fit_D1 <- readRDS(file.path(nmf_dir, "delta_D1_nmf.rds"))
fit_D2 <- readRDS(file.path(nmf_dir, "delta_D2_laplace.rds"))
fit_D3 <- readRDS(file.path(nmf_dir, "delta_D3_normal.rds"))
fit_D4 <- readRDS(file.path(nmf_dir, "delta_D4_gbcd.rds"))
fit_D5 <- readRDS(file.path(nmf_dir, "delta_D5_asym_ebmf.rds"))

# Ensure F_pm rownames (protein names) are set — flashier inherits from colnames
# of the training matrix; set explicitly as a safeguard.
if (is.null(rownames(fit_D2$F_pm))) rownames(fit_D2$F_pm) <- colnames(X_prot)
if (is.null(rownames(fit_D3$F_pm))) rownames(fit_D3$F_pm) <- colnames(X_prot)
if (is.null(rownames(fit_D5$F_pm))) rownames(fit_D5$F_pm) <- colnames(X_prot)

# ── 2. Reconstruct delta metadata (120 patients) ───────────────────────────────

meta   <- pca_scores[match(rownames(X_prot), pca_scores$sample), ]
T1_idx <- which(meta$time_point == "T1")
T2_idx <- which(meta$time_point == "T2")
pid_T1 <- meta$patient_id[T1_idx]
pid_T2 <- meta$patient_id[T2_idx]
shared_patients <- intersect(pid_T1, pid_T2)    # 120 patients
T1_order   <- T1_idx[match(shared_patients, pid_T1)]
meta_delta <- meta[T1_order, ]
rownames(meta_delta) <- shared_patients

who_delta_vec <- meta_delta$who_delta
who_score_T1  <- meta_delta$who_score_num   # severity at admission (T1)

# Verify L_pm row alignment (if rownames were preserved by flashier)
if (!is.null(rownames(fit_D2$L_pm))) {
  if (!all(rownames(fit_D2$L_pm) == rownames(meta_delta))) {
    # Re-order meta_delta to match L_pm row order
    cat("Note: reordering meta_delta to match fit_D2$L_pm row order\n")
    idx <- match(rownames(fit_D2$L_pm), rownames(meta_delta))
    meta_delta    <- meta_delta[idx, ]
    who_delta_vec <- who_delta_vec[idx]
    who_score_T1  <- who_score_T1[idx]
  }
}

# ── 3. Patient groupings for flashier plots ────────────────────────────────────

who_delta_cat <- ifelse(
  is.na(who_delta_vec), "Ambiguous",
  ifelse(who_delta_vec < 0, "Improved",
  ifelse(who_delta_vec == 0, "Stable", "Worsened"))
)

sev_t1_chr    <- as.character(who_score_T1)
sev_t1_levels <- sort(unique(sev_t1_chr))   # "1","1.5","3","4","5","6","7"

# flashier's plot() re-factors pm_groups via factor() and sorts levels
# alphabetically; pm_colors must be in that same sorted order.
delta_cat_levels <- sort(unique(who_delta_cat))  # Ambiguous,Improved,Stable,Worsened
delta_groups     <- factor(who_delta_cat, levels = delta_cat_levels)
sev_groups       <- factor(sev_t1_chr, levels = sev_t1_levels)

# ── 4. Colour palettes ─────────────────────────────────────────────────────────

# who_delta_cat: named colours, then extract in sorted alphabetical level order
delta_named <- c(Ambiguous = "gray70", Improved = "steelblue",
                 Stable    = "gold3",  Worsened = "firebrick")
delta_cols_vec <- unname(delta_named[delta_cat_levels])

# T1 severity gradient: light → dark blue (one per sorted level)
sev_colors <- colorRampPalette(c("#F7FBFF", "#08306B"))(length(sev_t1_levels))

# ── 5. PVE ordering (D2, D3) ──────────────────────────────────────────────────

D2_col_ord <- order(fit_D2$pve, decreasing = TRUE)
D3_col_ord <- order(fit_D3$pve, decreasing = TRUE)
D5_col_ord <- order(fit_D5$pve, decreasing = TRUE)

cat(sprintf("D2 K=%d, D3 K=%d, D5 K=%d\n",
            length(fit_D2$pve), length(fit_D3$pve), length(fit_D5$pve)))

# ── 6. Association rho helper ──────────────────────────────────────────────────

get_rho_delta <- function(variant_str, col_ord = NULL) {
  df  <- assoc[assoc$variant == variant_str, ]
  df  <- df[order(df$factor), ]
  rho <- abs(df$spearman_rho)
  rho[is.na(rho)] <- 0
  if (!is.null(col_ord)) rho <- rho[col_ord]
  rho
}

# ══════════════════════════════════════════════════════════════════════════════
# A.  F-MATRIX PHEATMAPS  (protein × factor, all four variants)
# ══════════════════════════════════════════════════════════════════════════════

col_norm <- function(M) {
  cm <- apply(M, 2, max)
  cm[cm == 0] <- 1
  sweep(M, 2, cm, "/")
}

div_scale_F <- function(M) {
  lim  <- ceiling(max(abs(M), na.rm = TRUE) * 10) / 10
  cols <- colorRampPalette(c("steelblue", "white", "red"))(101)
  brks <- seq(-lim, lim, length.out = 102)
  list(cols = cols, brks = brks)
}

make_col_anno_F <- function(rho_vec, K) {
  data.frame(abs_rho_who_delta = round(rho_vec, 3),
             row.names         = paste0("F", seq_len(K)))
}

heat_cols_pos   <- colorRampPalette(c("gray96", "red"))(50)
heat_brks_pos   <- seq(0, 1, length.out = 51)
rho_col_fn      <- colorRampPalette(c("white", "#E6550D"))(100)
anno_col_colors <- list(abs_rho_who_delta = rho_col_fn)

# ── A1: D1 Standard NMF — F heatmap (non-negative, col-normalised) ────────────

F_D1_raw <- NMF::basis(fit_D1)              # 481 × 6
rownames(F_D1_raw) <- colnames(X_prot)
colnames(F_D1_raw) <- paste0("F", seq_len(ncol(F_D1_raw)))
F_D1 <- col_norm(F_D1_raw)

anno_col_F_D1 <- make_col_anno_F(get_rho_delta("D1_NMF"), ncol(F_D1))

pdf(file.path(fig_dir, "heatmap_F_D1.pdf"), width = 6, height = 9)
pheatmap(
  F_D1,
  cluster_rows      = TRUE,
  clustering_method = "ward.D2",
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  annotation_col    = anno_col_F_D1,
  annotation_colors = anno_col_colors,
  angle_col         = 45,
  fontsize          = 9,
  color             = heat_cols_pos,
  breaks            = heat_brks_pos,
  main              = "D1 Standard NMF (K=6, shifted delta)\nF matrix — protein weights, col-normalised [0,1]",
  border_color      = NA
)
dev.off()
cat("Saved: heatmap_F_D1.pdf\n")

# ── A2: D2 EBMF Laplace — F heatmap (signed, PVE order) ─────────────────────

F_D2_raw <- fit_D2$F_pm                               # 481 × 20
rownames(F_D2_raw) <- colnames(X_prot)
F_D2_ord <- F_D2_raw[, D2_col_ord, drop = FALSE]
colnames(F_D2_ord) <- paste0("F", seq_len(ncol(F_D2_ord)))

sc_D2_F <- div_scale_F(F_D2_ord)
anno_col_F_D2 <- make_col_anno_F(
  get_rho_delta("D2_EBMF_Laplace", col_ord = D2_col_ord), ncol(F_D2_ord))

pdf(file.path(fig_dir, "heatmap_F_D2.pdf"), width = 10, height = 9)
pheatmap(
  F_D2_ord,
  cluster_rows      = TRUE,
  clustering_method = "ward.D2",
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  annotation_col    = anno_col_F_D2,
  annotation_colors = anno_col_colors,
  angle_col         = 45,
  fontsize          = 9,
  color             = sc_D2_F$cols,
  breaks            = sc_D2_F$brks,
  main              = "D2 EBMF point-Laplace (K=20, PVE order)\nF matrix — signed protein weights (blue=neg, red=pos)",
  border_color      = NA
)
dev.off()
cat("Saved: heatmap_F_D2.pdf\n")

# ── A3: D3 EBMF Normal — F heatmap (signed, PVE order) ──────────────────────

F_D3_raw <- fit_D3$F_pm                               # 481 × 20
rownames(F_D3_raw) <- colnames(X_prot)
F_D3_ord <- F_D3_raw[, D3_col_ord, drop = FALSE]
colnames(F_D3_ord) <- paste0("F", seq_len(ncol(F_D3_ord)))

sc_D3_F <- div_scale_F(F_D3_ord)
anno_col_F_D3 <- make_col_anno_F(
  get_rho_delta("D3_EBMF_Normal", col_ord = D3_col_ord), ncol(F_D3_ord))

pdf(file.path(fig_dir, "heatmap_F_D3.pdf"), width = 10, height = 9)
pheatmap(
  F_D3_ord,
  cluster_rows      = TRUE,
  clustering_method = "ward.D2",
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  annotation_col    = anno_col_F_D3,
  annotation_colors = anno_col_colors,
  angle_col         = 45,
  fontsize          = 9,
  color             = sc_D3_F$cols,
  breaks            = sc_D3_F$brks,
  main              = "D3 EBMF point-normal (K=20, PVE order)\nF matrix — signed protein weights (blue=neg, red=pos)",
  border_color      = NA
)
dev.off()
cat("Saved: heatmap_F_D3.pdf\n")

# ── A4: D4 GBCD — F heatmap (signed LFC, diverging) ─────────────────────────

F_D4_lfc <- fit_D4$F$lfc                             # 481 × 18
rownames(F_D4_lfc) <- colnames(X_prot)
colnames(F_D4_lfc) <- paste0("F", seq_len(ncol(F_D4_lfc)))

sc_D4_F <- div_scale_F(F_D4_lfc)
anno_col_F_D4 <- make_col_anno_F(get_rho_delta("D4_GBCD"), ncol(F_D4_lfc))

pdf(file.path(fig_dir, "heatmap_F_D4.pdf"), width = 9, height = 9)
pheatmap(
  F_D4_lfc,
  cluster_rows      = TRUE,
  clustering_method = "ward.D2",
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  annotation_col    = anno_col_F_D4,
  annotation_colors = anno_col_colors,
  angle_col         = 45,
  fontsize          = 9,
  color             = sc_D4_F$cols,
  breaks            = sc_D4_F$brks,
  main              = "D4 GBCD (K=18)\nF matrix — signed LFC (blue=decrease, red=increase)",
  border_color      = NA
)
dev.off()
cat("Saved: heatmap_F_D4.pdf\n")

# ── A5: D5 Asymmetric EBMF — F heatmap (signed F, PVE order) ─────────────────

F_D5_raw <- fit_D5$F_pm                               # 481 × 20 (signed)
rownames(F_D5_raw) <- colnames(X_prot)
F_D5_ord <- F_D5_raw[, D5_col_ord, drop = FALSE]
colnames(F_D5_ord) <- paste0("F", seq_len(ncol(F_D5_ord)))

sc_D5_F <- div_scale_F(F_D5_ord)
anno_col_F_D5 <- make_col_anno_F(
  get_rho_delta("D5_Asym_EBMF", col_ord = D5_col_ord), ncol(F_D5_ord))

pdf(file.path(fig_dir, "heatmap_F_D5.pdf"), width = 10, height = 9)
pheatmap(
  F_D5_ord,
  cluster_rows      = TRUE,
  clustering_method = "ward.D2",
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  annotation_col    = anno_col_F_D5,
  annotation_colors = anno_col_colors,
  angle_col         = 45,
  fontsize          = 9,
  color             = sc_D5_F$cols,
  breaks            = sc_D5_F$brks,
  main              = "D5 Asymmetric EBMF (K=20, PVE order)\nF matrix — signed protein weights (blue=neg, red=pos)",
  border_color      = NA
)
dev.off()
cat("Saved: heatmap_F_D5.pdf\n")

cat("\n[A] F-matrix pheatmaps done.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# B.  FLASHIER SINGLE-CELL VIGNETTE STYLE  (fit_D2 only — PRIMARY)
# ══════════════════════════════════════════════════════════════════════════════
# All calls use plot(fit_D2, plot_type = ..., pm_which = ..., pm_groups = ...,
# pm_colors = ...).  pm_colors must be an unnamed vector of length =
# length(unique(pm_groups)), in SORTED level order (flashier re-factors
# pm_groups internally via factor() before indexing pm_colors).
#
# Patient colourings:
#   delta_groups  : Improved / Stable / Worsened / Ambiguous (sorted: Amb,Imp,Sta,Wor)
#   sev_groups    : T1 severity 1, 1.5, 3, 4, 5, 6, 7 (sorted numerically as chr)

# ── B1: Scree / PVE bar chart ─────────────────────────────────────────────────

p_scree <- plot(fit_D2, plot_type = "scree")
ggsave(file.path(fig_dir, "d2_B1_scree.pdf"), p_scree, width = 6, height = 4)
cat("Saved: d2_B1_scree.pdf\n")

# ── B2: Histograms of loadings ─────────────────────────────────────────────────

p_hist_dc <- plot(fit_D2,
                  plot_type = "histogram",
                  pm_which  = "loadings",
                  pm_groups = delta_groups,
                  pm_colors = delta_cols_vec,
                  bins      = 20)
ggsave(file.path(fig_dir, "d2_B2a_histogram_delta_cat.pdf"),
       p_hist_dc, width = 14, height = 8)
cat("Saved: d2_B2a_histogram_delta_cat.pdf\n")

p_hist_sv <- plot(fit_D2,
                  plot_type = "histogram",
                  pm_which  = "loadings",
                  pm_groups = sev_groups,
                  pm_colors = sev_colors,
                  bins      = 20)
ggsave(file.path(fig_dir, "d2_B2b_histogram_severity.pdf"),
       p_hist_sv, width = 14, height = 8)
cat("Saved: d2_B2b_histogram_severity.pdf\n")

# ── B3: Structure plots ────────────────────────────────────────────────────────
# SKIPPED for fit_D2: flashier's structure_plot() calls verify.nonnegative.matrix()
# internally and rejects signed L matrices (point-Laplace loadings can be negative).
# Structure plots are only valid for non-negative factorizations (NMF/point-exponential).
# See 05_factor_plots.R (fit_S2, point-exponential) for the stacked-analysis equivalent.
cat("B3: structure plots skipped — requires non-negative L (point-Laplace is signed)\n")

# ── B4: L heatmaps ────────────────────────────────────────────────────────────
# SKIPPED for fit_D2: plot(, plot_type="heatmap", pm_which="loadings") calls
# select_loadings() → verify.fit() → verify.nonnegative.matrix(), which rejects
# signed L matrices. Use 06_delta_loading_heatmaps.R outputs instead:
#   results/figs/nmf/heatmap_delta_D2_loadings.pdf (pheatmap, diverging scale)
cat("B4: L heatmap skipped — signed L not supported; see 06_delta_loading_heatmaps.R output\n")

# ── B5: Top-protein F heatmap ─────────────────────────────────────────────────
# D2 most-associated factors (ORIGINAL indices, from delta_association.csv):
#   F17 (rho=+0.279, FDR=0.213), F5 (rho=+0.259), F2 (rho=+0.250), F4 (rho=-0.245)
# Top 5 proteins per factor by |F_pm| (signed — use abs() for ranking).

kset_key <- c(17, 5, 2, 4)
top_idx  <- unique(unlist(lapply(kset_key, function(k) {
  order(abs(fit_D2$F_pm[, k]), decreasing = TRUE)[1:5]
})))
top_prot   <- rownames(fit_D2$F_pm)[top_idx]
# Note: pm_groups cannot be used here — even with pm_which="factors", passing pm_groups
# triggers select_loadings() which calls verify.nonnegative.matrix() and rejects signed L.
# Use basic call; flashier uses t-SNE for protein ordering.

p_heat_F_top <- plot(fit_D2,
                     plot_type = "heatmap",
                     pm_which  = "factors",
                     pm_subset = top_prot,
                     kset      = kset_key)
ggsave(file.path(fig_dir, "d2_B5_heatmap_F_top_proteins.pdf"),
       p_heat_F_top, width = 8, height = 6)
cat("Saved: d2_B5_heatmap_F_top_proteins.pdf\n")

# ── B6: F scatter plots for key factors ───────────────────────────────────────
# x = protein weight (F_pm[,k], signed for point-Laplace)
# y = column mean of training data = mean (T2-T1) per protein across 120 patients

for (k in c(17, 5, 2)) {
  p_scatter <- plot(fit_D2,
                    plot_type  = "scatter",
                    pm_which   = "factors",
                    kset       = k,
                    labels     = TRUE,
                    n_labels   = 10,
                    label_size = 2.5) +
    labs(x     = "protein weight (F_pm, signed)",
         y     = "mean protein delta (T2-T1) across 120 patients",
         title = sprintf("D2 EBMF Laplace — Factor %d: protein weights vs mean delta", k))
  ggsave(file.path(fig_dir, sprintf("d2_B6_scatter_F%02d.pdf", k)),
         p_scatter, width = 6, height = 5)
  cat(sprintf("Saved: d2_B6_scatter_F%02d.pdf\n", k))
}

cat("\n[B] flashier single-cell vignette plots done.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# C.  FLASHIER_INTRO VIGNETTE STYLE — BAR PLOTS OF LOADINGS
# ══════════════════════════════════════════════════════════════════════════════

p_bar_dc <- plot(fit_D2,
                 plot_type = "bar",
                 pm_which  = "loadings",
                 pm_groups = delta_groups,
                 pm_colors = delta_cols_vec)
ggsave(file.path(fig_dir, "d2_C1_bar_loadings_delta_cat.pdf"),
       p_bar_dc, width = 16, height = 8)
cat("Saved: d2_C1_bar_loadings_delta_cat.pdf\n")

p_bar_sv <- plot(fit_D2,
                 plot_type = "bar",
                 pm_which  = "loadings",
                 pm_groups = sev_groups,
                 pm_colors = sev_colors)
ggsave(file.path(fig_dir, "d2_C2_bar_loadings_severity.pdf"),
       p_bar_sv, width = 16, height = 8)
cat("Saved: d2_C2_bar_loadings_severity.pdf\n")

cat("\n[C] flashier_intro bar plots done.\n")

# ══════════════════════════════════════════════════════════════════════════════
# D.  FLASHIER PLOTS FOR D5 — ALL TYPES ENABLED (non-negative L)
# ══════════════════════════════════════════════════════════════════════════════
# D5 has point-exponential prior on L → L >= 0 → structure, L heatmap, bar all work.
# D5 most-associated factors: F1 (rho=-0.256), F12 (rho=-0.218), F18 (rho=-0.214)

# ── D1: Scree / PVE bar chart ─────────────────────────────────────────────────

p_d5_scree <- plot(fit_D5, plot_type = "scree")
ggsave(file.path(fig_dir, "d5_D1_scree.pdf"), p_d5_scree, width = 6, height = 4)
cat("Saved: d5_D1_scree.pdf\n")

# ── D2: Histograms of loadings ─────────────────────────────────────────────────

p_d5_hist_dc <- plot(fit_D5,
                     plot_type = "histogram",
                     pm_which  = "loadings",
                     pm_groups = delta_groups,
                     pm_colors = delta_cols_vec,
                     bins      = 20)
ggsave(file.path(fig_dir, "d5_D2a_histogram_delta_cat.pdf"),
       p_d5_hist_dc, width = 14, height = 8)
cat("Saved: d5_D2a_histogram_delta_cat.pdf\n")

p_d5_hist_sv <- plot(fit_D5,
                     plot_type = "histogram",
                     pm_which  = "loadings",
                     pm_groups = sev_groups,
                     pm_colors = sev_colors,
                     bins      = 20)
ggsave(file.path(fig_dir, "d5_D2b_histogram_severity.pdf"),
       p_d5_hist_sv, width = 14, height = 8)
cat("Saved: d5_D2b_histogram_severity.pdf\n")

# ── D3: Structure plots (enabled by non-neg L) ────────────────────────────────

p_d5_str_dc <- plot(fit_D5,
                    plot_type = "structure",
                    pm_which  = "loadings",
                    pm_groups = delta_groups,
                    gap       = 5)
ggsave(file.path(fig_dir, "d5_D3a_structure_delta_cat.pdf"),
       p_d5_str_dc, width = 14, height = 5)
cat("Saved: d5_D3a_structure_delta_cat.pdf\n")

p_d5_str_sv <- plot(fit_D5,
                    plot_type = "structure",
                    pm_which  = "loadings",
                    pm_groups = sev_groups,
                    gap       = 3)
ggsave(file.path(fig_dir, "d5_D3b_structure_severity.pdf"),
       p_d5_str_sv, width = 14, height = 5)
cat("Saved: d5_D3b_structure_severity.pdf\n")

# ── D4: L heatmaps (enabled by non-neg L) ────────────────────────────────────

p_d5_heat_L_dc <- plot(fit_D5,
                       plot_type = "heatmap",
                       pm_which  = "loadings",
                       pm_groups = delta_groups,
                       gap       = 5)
ggsave(file.path(fig_dir, "d5_D4a_heatmap_L_delta_cat.pdf"),
       p_d5_heat_L_dc, width = 10, height = 6)
cat("Saved: d5_D4a_heatmap_L_delta_cat.pdf\n")

p_d5_heat_L_sv <- plot(fit_D5,
                       plot_type = "heatmap",
                       pm_which  = "loadings",
                       pm_groups = sev_groups,
                       gap       = 3)
ggsave(file.path(fig_dir, "d5_D4b_heatmap_L_severity.pdf"),
       p_d5_heat_L_sv, width = 10, height = 6)
cat("Saved: d5_D4b_heatmap_L_severity.pdf\n")

# ── D5: Top-protein F heatmap ─────────────────────────────────────────────────
# D5 most-associated factors (original indices): F1, F12, F18, F4

kset_d5 <- c(1, 12, 18, 4)
top_idx_d5 <- unique(unlist(lapply(kset_d5, function(k) {
  order(abs(fit_D5$F_pm[, k]), decreasing = TRUE)[1:5]
})))
top_prot_d5 <- rownames(fit_D5$F_pm)[top_idx_d5]

p_d5_heat_F <- plot(fit_D5,
                    plot_type = "heatmap",
                    pm_which  = "factors",
                    pm_subset = top_prot_d5,
                    kset      = kset_d5)
ggsave(file.path(fig_dir, "d5_D5_heatmap_F_top_proteins.pdf"),
       p_d5_heat_F, width = 8, height = 6)
cat("Saved: d5_D5_heatmap_F_top_proteins.pdf\n")

# ── D6: F scatter plots for key factors ───────────────────────────────────────

for (k in c(1, 12, 18)) {
  p_scatter_d5 <- plot(fit_D5,
                       plot_type  = "scatter",
                       pm_which   = "factors",
                       kset       = k,
                       labels     = TRUE,
                       n_labels   = 10,
                       label_size = 2.5) +
    labs(x     = "protein weight (F_pm, signed)",
         y     = "mean protein delta (T2-T1) across 120 patients",
         title = sprintf("D5 Asym EBMF — Factor %d: protein weights vs mean delta", k))
  ggsave(file.path(fig_dir, sprintf("d5_D6_scatter_F%02d.pdf", k)),
         p_scatter_d5, width = 6, height = 5)
  cat(sprintf("Saved: d5_D6_scatter_F%02d.pdf\n", k))
}

# ── D7: Bar plots of loadings (non-neg L) ────────────────────────────────────

p_d5_bar_dc <- plot(fit_D5,
                    plot_type = "bar",
                    pm_which  = "loadings",
                    pm_groups = delta_groups,
                    pm_colors = delta_cols_vec)
ggsave(file.path(fig_dir, "d5_D7a_bar_loadings_delta_cat.pdf"),
       p_d5_bar_dc, width = 16, height = 8)
cat("Saved: d5_D7a_bar_loadings_delta_cat.pdf\n")

p_d5_bar_sv <- plot(fit_D5,
                    plot_type = "bar",
                    pm_which  = "loadings",
                    pm_groups = sev_groups,
                    pm_colors = sev_colors)
ggsave(file.path(fig_dir, "d5_D7b_bar_loadings_severity.pdf"),
       p_d5_bar_sv, width = 16, height = 8)
cat("Saved: d5_D7b_bar_loadings_severity.pdf\n")

cat("\n[D] D5 flashier plots done (all plot types — non-neg L).\n")

# ── Summary ────────────────────────────────────────────────────────────────────

cat("\n--- 07_delta_factor_plots.R complete ---\n")
cat("\nOutputs written to:", fig_dir, "\n")
cat("A. F-matrix pheatmaps:\n")
cat("   heatmap_F_D1.pdf  heatmap_F_D2.pdf  heatmap_F_D3.pdf\n")
cat("   heatmap_F_D4.pdf  heatmap_F_D5.pdf\n")
cat("B. flashier vignette style (D2, signed L — limited types):\n")
cat("   d2_B1_scree.pdf\n")
cat("   d2_B2a_histogram_delta_cat.pdf, d2_B2b_histogram_severity.pdf\n")
cat("   d2_B3/B4: SKIPPED (structure/L-heatmap require non-neg L)\n")
cat("   d2_B5_heatmap_F_top_proteins.pdf\n")
cat("   d2_B6_scatter_F17.pdf, d2_B6_scatter_F05.pdf, d2_B6_scatter_F02.pdf\n")
cat("C. flashier_intro bar plots (D2):\n")
cat("   d2_C1_bar_loadings_delta_cat.pdf, d2_C2_bar_loadings_severity.pdf\n")
cat("D. flashier plots (D5, non-neg L — ALL types enabled):\n")
cat("   d5_D1_scree.pdf\n")
cat("   d5_D2a_histogram_delta_cat.pdf, d5_D2b_histogram_severity.pdf\n")
cat("   d5_D3a_structure_delta_cat.pdf, d5_D3b_structure_severity.pdf\n")
cat("   d5_D4a_heatmap_L_delta_cat.pdf, d5_D4b_heatmap_L_severity.pdf\n")
cat("   d5_D5_heatmap_F_top_proteins.pdf\n")
cat("   d5_D6_scatter_F01.pdf, d5_D6_scatter_F12.pdf, d5_D6_scatter_F18.pdf\n")
cat("   d5_D7a_bar_loadings_delta_cat.pdf, d5_D7b_bar_loadings_severity.pdf\n")
