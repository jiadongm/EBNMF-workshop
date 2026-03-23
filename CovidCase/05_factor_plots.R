# 05_factor_plots.R — Flashier-vignette-style factor visualisations
#
# Three blocks of plots, all for the stacked 240×481 proteomics analysis.
#
# A. F-matrix pheatmaps (protein × factor) for all three NMF variants
#    A1  S1 standard NMF  (481×8,  col-normalised [0,1])
#    A2  S2 EBNMF         (481×14, col-normalised [0,1], PVE order)
#    A3  S3 GBCD          (481×10, signed LFC, diverging blue-white-red)
#
# B. flashier single-cell vignette style (fit_S2 only)
#    B1  Scree / PVE bar chart
#    B2  Histograms of loadings    × 2 colourings (time point / severity)
#    B3  Structure plots           × 2 colourings
#    B4  L heatmaps                × 2 colourings
#    B5  Top-protein F heatmap (factors 5, 1, 8 — most clinically associated)
#    B6  F scatter plots for factors 5, 1, 8
#
# C. flashier_intro vignette style — bar plots of L × 2 colourings
#
# References:
#   https://willwerscheid.github.io/flashier/articles/flashier_single_cell.html
#   https://willwerscheid.github.io/flashier/articles/flashier_intro.html
#
# Run from CovidCase/:   Rscript 05_factor_plots.R

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

X_prot     <- readRDS(file.path(results_dir, "X_prot.rds"))        # 240 × 481
pca_scores <- readRDS(file.path(results_dir, "pca_scores.rds"))
assoc      <- read.csv(file.path(nmf_dir, "stacked_association.csv"))

fit_S1 <- readRDS(file.path(nmf_dir, "stacked_S1_nmf.rds"))
fit_S2 <- readRDS(file.path(nmf_dir, "stacked_S2_ebnmf.rds"))
fit_S3 <- readRDS(file.path(nmf_dir, "stacked_S3_gbcd.rds"))

# ── 2. Metadata ────────────────────────────────────────────────────────────────

meta <- pca_scores[match(rownames(X_prot), pca_scores$sample), ]
stopifnot(all(meta$sample == rownames(X_prot)))

# Sorted severity levels present in the data
sev_raw    <- as.character(meta$who_score_num)  # "1","1.5","2","3","4","5","6","7"
sev_uniq   <- sort(unique(sev_raw))             # sorted unique levels
sev_groups <- factor(sev_raw, levels = sev_uniq)

# Time-point grouping
tp_groups  <- factor(meta$time_point, levels = c("T1", "T2"))

# ── 3. Colour palettes ─────────────────────────────────────────────────────────

# Time point
tp_colors  <- c("steelblue", "firebrick")    # T1 first (level order), T2 second

# Severity: gradient light → dark blue, one colour per unique level in sorted order
sev_colors <- colorRampPalette(c("#F7FBFF", "#08306B"))(length(sev_uniq))

# ── 4. PVE ordering (S2) ──────────────────────────────────────────────────────

pve_S2     <- fit_S2$pve
S2_col_ord <- order(pve_S2, decreasing = TRUE)

# ── 5. Association rho helper ─────────────────────────────────────────────────

get_rho <- function(variant_str, col_ord = NULL) {
  df  <- assoc[assoc$variant == variant_str, ]
  df  <- df[order(df$factor), ]
  rho <- abs(df$spearman_rho)
  rho[is.na(rho)] <- 0
  if (!is.null(col_ord)) rho <- rho[col_ord]
  rho
}

# ══════════════════════════════════════════════════════════════════════════════
# A.  F-MATRIX PHEATMAPS  (protein × factor, all three variants)
# ══════════════════════════════════════════════════════════════════════════════

# ── Shared pheatmap helpers ────────────────────────────────────────────────────

col_norm <- function(M) {
  cm <- apply(M, 2, max)
  cm[cm == 0] <- 1
  sweep(M, 2, cm, "/")
}

make_col_anno_F <- function(rho_vec, K) {
  data.frame(abs_rho_who_delta = round(rho_vec, 3),
             row.names         = paste0("F", seq_len(K)))
}

heat_cols_pos   <- colorRampPalette(c("gray96", "red"))(50)
heat_brks_pos   <- seq(0, 1, length.out = 51)
rho_col_fn      <- colorRampPalette(c("white", "#E6550D"))(100)
anno_col_colors <- list(abs_rho_who_delta = rho_col_fn)

# ── A1: S1 Standard NMF — F heatmap ──────────────────────────────────────────

F_S1_raw <- NMF::basis(fit_S1)                   # 481 × 8
rownames(F_S1_raw) <- colnames(X_prot)
colnames(F_S1_raw) <- paste0("F", seq_len(ncol(F_S1_raw)))
F_S1 <- col_norm(F_S1_raw)

anno_col_F_S1 <- make_col_anno_F(get_rho("S1_NMF"), ncol(F_S1))

pdf(file.path(fig_dir, "heatmap_F_S1.pdf"), width = 6, height = 9)
pheatmap(
  F_S1,
  cluster_rows      = TRUE,
  clustering_method = "ward.D2",
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  annotation_col    = anno_col_F_S1,
  annotation_colors = anno_col_colors,
  angle_col         = 45,
  fontsize          = 9,
  color             = heat_cols_pos,
  breaks            = heat_brks_pos,
  main              = "S1 Standard NMF (K=8)\nF matrix — protein weights, col-normalised [0,1]",
  border_color      = NA
)
dev.off()
cat("Saved: heatmap_F_S1.pdf\n")

# ── A2: S2 EBNMF — F heatmap (columns in PVE order) ─────────────────────────

F_S2_raw <- fit_S2$F_pm                          # 481 × 14 (rownames already set)
F_S2_ord <- F_S2_raw[, S2_col_ord, drop = FALSE] # reorder columns by PVE
colnames(F_S2_ord) <- paste0("F", seq_len(ncol(F_S2_ord)))
F_S2 <- col_norm(F_S2_ord)

anno_col_F_S2 <- make_col_anno_F(get_rho("S2_EBNMF", col_ord = S2_col_ord), ncol(F_S2))

pdf(file.path(fig_dir, "heatmap_F_S2.pdf"), width = 8, height = 9)
pheatmap(
  F_S2,
  cluster_rows      = TRUE,
  clustering_method = "ward.D2",
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  annotation_col    = anno_col_F_S2,
  annotation_colors = anno_col_colors,
  angle_col         = 45,
  fontsize          = 9,
  color             = heat_cols_pos,
  breaks            = heat_brks_pos,
  main              = "S2 EBNMF point-exponential (K=14, PVE order)\nF matrix — protein weights, col-normalised [0,1]",
  border_color      = NA
)
dev.off()
cat("Saved: heatmap_F_S2.pdf\n")

# ── A3: S3 GBCD — F heatmap (signed LFC, diverging colour) ───────────────────

F_S3_lfc <- fit_S3$F$lfc                         # 481 × 10 (rownames already set)
colnames(F_S3_lfc) <- paste0("F", seq_len(ncol(F_S3_lfc)))

lfc_lim        <- ceiling(max(abs(F_S3_lfc), na.rm = TRUE) * 10) / 10
heat_cols_div  <- colorRampPalette(c("steelblue", "white", "red"))(101)
heat_brks_div  <- seq(-lfc_lim, lfc_lim, length.out = 102)

anno_col_F_S3 <- make_col_anno_F(get_rho("S3_GBCD"), ncol(F_S3_lfc))

pdf(file.path(fig_dir, "heatmap_F_S3.pdf"), width = 7, height = 9)
pheatmap(
  F_S3_lfc,
  cluster_rows      = TRUE,
  clustering_method = "ward.D2",
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  annotation_col    = anno_col_F_S3,
  annotation_colors = anno_col_colors,
  angle_col         = 45,
  fontsize          = 9,
  color             = heat_cols_div,
  breaks            = heat_brks_div,
  main              = "S3 GBCD (K=10)\nF matrix — signed LFC (blue = decrease, red = increase)",
  border_color      = NA
)
dev.off()
cat("Saved: heatmap_F_S3.pdf\n")

cat("\n[A] F-matrix pheatmaps done.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# B.  FLASHIER SINGLE-CELL VIGNETTE STYLE  (fit_S2 only)
# ══════════════════════════════════════════════════════════════════════════════
# All calls use plot(fit_S2, plot_type = ..., pm_which = ..., pm_groups = ...,
# pm_colors = ...).  pm_colors must be an unnamed vector of length =
# length(unique(pm_groups)), in SORTED level order (flashier re-factors
# pm_groups internally via factor() before indexing pm_colors).

# ── B1: Scree / PVE bar chart ─────────────────────────────────────────────────

p_scree <- plot(fit_S2, plot_type = "scree")
ggsave(file.path(fig_dir, "s2_B1_scree.pdf"), p_scree, width = 6, height = 4)
cat("Saved: s2_B1_scree.pdf\n")

# ── B2: Histograms of loadings ─────────────────────────────────────────────────

p_hist_tp <- plot(fit_S2,
                  plot_type = "histogram",
                  pm_which  = "loadings",
                  pm_groups = tp_groups,
                  pm_colors = tp_colors,
                  bins      = 20)
ggsave(file.path(fig_dir, "s2_B2a_histogram_timepoint.pdf"),
       p_hist_tp, width = 14, height = 8)
cat("Saved: s2_B2a_histogram_timepoint.pdf\n")

p_hist_sv <- plot(fit_S2,
                  plot_type = "histogram",
                  pm_which  = "loadings",
                  pm_groups = sev_groups,
                  pm_colors = sev_colors,
                  bins      = 20)
ggsave(file.path(fig_dir, "s2_B2b_histogram_severity.pdf"),
       p_hist_sv, width = 14, height = 8)
cat("Saved: s2_B2b_histogram_severity.pdf\n")

# ── B3: Structure plots ────────────────────────────────────────────────────────
# Structure plot: each patient is a stacked bar; row-normalised L is shown.
# gap = gap between groups (number of samples wide).

p_str_tp <- plot(fit_S2,
                 plot_type = "structure",
                 pm_which  = "loadings",
                 pm_groups = tp_groups,
                 # pm_colors not used: structure plot uses pm_colors for FACTOR
                 # fill colours (K=14 needed), not for patient groups. Groups
                 # are distinguished by gaps and x-axis labels.
                 gap       = 10)
ggsave(file.path(fig_dir, "s2_B3a_structure_timepoint.pdf"),
       p_str_tp, width = 14, height = 5)
cat("Saved: s2_B3a_structure_timepoint.pdf\n")

p_str_sv <- plot(fit_S2,
                 plot_type = "structure",
                 pm_which  = "loadings",
                 pm_groups = sev_groups,
                 gap       = 5)
ggsave(file.path(fig_dir, "s2_B3b_structure_severity.pdf"),
       p_str_sv, width = 14, height = 5)
cat("Saved: s2_B3b_structure_severity.pdf\n")

# ── B4: L heatmaps ────────────────────────────────────────────────────────────
# flashier heatmap: samples on y-axis, factors on x-axis, fill = loading value.
# pm_colors here sets the fill gradient (low/mid/high), not the group colours.
# Leave pm_colors as default (darkred / white / darkblue) — or supply a 2- or
# 3-element vector.  Group colouring appears via y-axis tick labels coloured by
# pm_groups (handled internally).

p_heat_L_tp <- plot(fit_S2,
                    plot_type = "heatmap",
                    pm_which  = "loadings",
                    pm_groups = tp_groups,
                    gap       = 10)
ggsave(file.path(fig_dir, "s2_B4a_heatmap_L_timepoint.pdf"),
       p_heat_L_tp, width = 10, height = 6)
cat("Saved: s2_B4a_heatmap_L_timepoint.pdf\n")

p_heat_L_sv <- plot(fit_S2,
                    plot_type = "heatmap",
                    pm_which  = "loadings",
                    pm_groups = sev_groups,
                    gap       = 5)
ggsave(file.path(fig_dir, "s2_B4b_heatmap_L_severity.pdf"),
       p_heat_L_sv, width = 10, height = 6)
cat("Saved: s2_B4b_heatmap_L_severity.pdf\n")

# ── B5: Top-protein F heatmap ─────────────────────────────────────────────────
# Show top 5 proteins (by F_pm weight) for the most clinically associated
# factors, using the ORIGINAL S2 factor indices:
#   F5 (rho=-0.344, FDR=0.054), F1 (rho=-0.273), F8 (rho=+0.226),
#   F2 (rho=+0.209), F7 (rho=+0.205)

kset_key   <- c(5, 1, 8, 2, 7)
top_idx    <- unique(unlist(lapply(kset_key, function(k) {
  order(fit_S2$F_pm[, k], decreasing = TRUE)[1:5]
})))
top_prot   <- rownames(fit_S2$F_pm)[top_idx]   # character vector of protein names

# pm_groups = protein names, one per protein (so each protein is its own group)
# pm_colors must have length = length(unique(pm_groups)) = length(top_prot)
# Use gray for all so the group colouring doesn't dominate (the heatmap fill matters)
top_colors <- rep("gray50", length(top_prot))

p_heat_F_top <- plot(fit_S2,
                     plot_type = "heatmap",
                     pm_which  = "factors",
                     pm_subset = top_prot,
                     pm_groups = factor(top_prot, levels = rev(top_prot)),
                     pm_colors = top_colors,
                     kset      = kset_key,
                     gap       = 0.2)
ggsave(file.path(fig_dir, "s2_B5_heatmap_F_top_proteins.pdf"),
       p_heat_F_top, width = 8, height = 6)
cat("Saved: s2_B5_heatmap_F_top_proteins.pdf\n")

# ── B6: F scatter plots — one per key factor ──────────────────────────────────
# x = protein weight (F_pm[,k])
# y = data column mean (mean NPX per protein across all 240 samples — computed
#     from the data stored inside the flashier fit)
# Top 10 proteins by |F_pm| are labelled.

for (k in c(5, 1, 8)) {
  p_scatter <- plot(fit_S2,
                    plot_type  = "scatter",
                    pm_which   = "factors",
                    kset       = k,
                    labels     = TRUE,
                    n_labels   = 10,
                    label_size = 2.5) +
    labs(x     = "protein weight (F_pm)",
         y     = "mean protein NPX across all 240 samples",
         title = sprintf("S2 EBNMF — Factor %d: protein weights vs mean expression", k))
  ggsave(file.path(fig_dir, sprintf("s2_B6_scatter_F%02d.pdf", k)),
         p_scatter, width = 6, height = 5)
  cat(sprintf("Saved: s2_B6_scatter_F%02d.pdf\n", k))
}

cat("\n[B] flashier single-cell vignette plots done.\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# C.  FLASHIER_INTRO VIGNETTE STYLE — BAR PLOTS OF LOADINGS
# ══════════════════════════════════════════════════════════════════════════════
# Each panel = one factor.  Each bar = one patient sample.
# Bars coloured by patient group (time point or severity).
# Mirrors:  plot(fl, pm_which="factors", pm_colors=gtex_colors, plot_type="bar")
# but applied to the loadings (L) side so patients are coloured.

p_bar_tp <- plot(fit_S2,
                 plot_type = "bar",
                 pm_which  = "loadings",
                 pm_groups = tp_groups,
                 pm_colors = tp_colors)
ggsave(file.path(fig_dir, "s2_C1_bar_loadings_timepoint.pdf"),
       p_bar_tp, width = 16, height = 8)
cat("Saved: s2_C1_bar_loadings_timepoint.pdf\n")

p_bar_sv <- plot(fit_S2,
                 plot_type = "bar",
                 pm_which  = "loadings",
                 pm_groups = sev_groups,
                 pm_colors = sev_colors)
ggsave(file.path(fig_dir, "s2_C2_bar_loadings_severity.pdf"),
       p_bar_sv, width = 16, height = 8)
cat("Saved: s2_C2_bar_loadings_severity.pdf\n")

cat("\n[C] flashier_intro bar plots done.\n")

# ── Summary ────────────────────────────────────────────────────────────────────

cat("\n--- 05_factor_plots.R complete ---\n")
cat("\nOutputs written to:", fig_dir, "\n")
cat("A. F-matrix pheatmaps:\n")
cat("   heatmap_F_S1.pdf, heatmap_F_S2.pdf, heatmap_F_S3.pdf\n")
cat("B. flashier single-cell vignette style (S2):\n")
cat("   s2_B1_scree.pdf\n")
cat("   s2_B2a_histogram_timepoint.pdf, s2_B2b_histogram_severity.pdf\n")
cat("   s2_B3a_structure_timepoint.pdf, s2_B3b_structure_severity.pdf\n")
cat("   s2_B4a_heatmap_L_timepoint.pdf, s2_B4b_heatmap_L_severity.pdf\n")
cat("   s2_B5_heatmap_F_top_proteins.pdf\n")
cat("   s2_B6_scatter_F05.pdf, s2_B6_scatter_F01.pdf, s2_B6_scatter_F08.pdf\n")
cat("C. flashier_intro style bar plots:\n")
cat("   s2_C1_bar_loadings_timepoint.pdf, s2_C2_bar_loadings_severity.pdf\n")
