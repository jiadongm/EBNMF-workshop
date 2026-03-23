# 04_loading_heatmaps.R — Pheatmap visualisation of NMF loading matrices
#
# Reproduces the style of the GBCD HNSCC vignette
# (https://stephenslab.github.io/gbcd/articles/hnscc.html):
#   - Rows = samples (240 patients × time points), ordered by time_point then WHO severity
#   - Columns = programs (factors), ordered by variance explained (S2) or as-fit (S1, S3)
#   - Color scale: gray96 → red, range [0, 1] (L columns normalised per-column)
#   - Row annotation: time point (T1/T2) + WHO severity score
#   - Column annotation: Spearman |rho| with who_delta (highlights associated factors)
#
# Inputs:
#   results/nmf/stacked_S1_nmf.rds
#   results/nmf/stacked_S2_ebnmf.rds
#   results/nmf/stacked_S3_gbcd.rds
#   results/nmf/stacked_association.csv
#   results/pca_scores.rds
#
# Outputs:
#   results/figs/nmf/heatmap_S1_loadings.pdf
#   results/figs/nmf/heatmap_S2_loadings.pdf
#   results/figs/nmf/heatmap_S3_loadings.pdf
#   results/figs/nmf/heatmap_all3_loadings.pdf   (combined 3-panel)

suppressPackageStartupMessages({
  library(NMF)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(gridExtra)
})

# ── Paths ──────────────────────────────────────────────────────────────────────

results_dir <- "results"
nmf_dir     <- file.path(results_dir, "nmf")
fig_dir     <- file.path(results_dir, "figs", "nmf")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load data ───────────────────────────────────────────────────────────────

pca_scores <- readRDS(file.path(results_dir, "pca_scores.rds"))
X_prot     <- readRDS(file.path(results_dir, "X_prot.rds"))
assoc      <- read.csv(file.path(nmf_dir, "stacked_association.csv"))

fit_S1 <- readRDS(file.path(nmf_dir, "stacked_S1_nmf.rds"))
fit_S2 <- readRDS(file.path(nmf_dir, "stacked_S2_ebnmf.rds"))
fit_S3 <- readRDS(file.path(nmf_dir, "stacked_S3_gbcd.rds"))

# ── 2. Extract and align L matrices ───────────────────────────────────────────

meta <- pca_scores[match(rownames(X_prot), pca_scores$sample), ]
stopifnot(all(meta$sample == rownames(X_prot)))

L_S1_raw <- t(NMF::coef(fit_S1))   # 240 × 8
L_S2_raw <- fit_S2$L_pm            # 240 × 14
L_S3_raw <- fit_S3$L               # 240 × 10  (already [0,1])

rownames(L_S1_raw) <- rownames(X_prot)
rownames(L_S2_raw) <- rownames(X_prot)
rownames(L_S3_raw) <- rownames(X_prot)

# ── 3. Per-column normalise to [0, 1] ─────────────────────────────────────────
# Follows the vignette: makes loadings comparable across factors on a common scale.
# S3 GBCD L is already in [0,1] by construction (GB prior); normalise anyway for
# consistency (column max may be < 1 for some factors).

col_norm <- function(L) {
  col_max <- apply(L, 2, max)
  col_max[col_max == 0] <- 1  # avoid division by zero (should not occur)
  sweep(L, 2, col_max, "/")
}

L_S1 <- col_norm(L_S1_raw)   # 240 × 8,  values in [0, 1]
L_S2 <- col_norm(L_S2_raw)   # 240 × 14, values in [0, 1]
L_S3 <- col_norm(L_S3_raw)   # 240 × 10, values in [0, 1]

# ── 4. Column ordering ────────────────────────────────────────────────────────
# S2: order factors by PVE (decreasing) so the dominant programs appear first.
# S1, S3: keep as-fit order (no PVE metric directly comparable).

pve_S2    <- fit_S2$pve
S2_col_ord <- order(pve_S2, decreasing = TRUE)
L_S2      <- L_S2[, S2_col_ord, drop = FALSE]

# Relabel columns
colnames(L_S1) <- paste0("F", seq_len(ncol(L_S1)))
colnames(L_S2) <- paste0("F", seq_len(ncol(L_S2)))   # now sorted by PVE
colnames(L_S3) <- paste0("F", seq_len(ncol(L_S3)))

# ── 5. Row ordering ───────────────────────────────────────────────────────────
# Mirror the vignette: order rows by a meaningful biological grouping.
# Here: time_point (T1 first, T2 second), then WHO severity ascending within each.

row_ord <- order(
  ifelse(meta$time_point == "T1", 0L, 1L),
  meta$who_score_num
)

# ── 6. Row annotation ─────────────────────────────────────────────────────────

anno_row <- data.frame(
  Time_point   = meta$time_point,
  WHO_severity = meta$who_score_num,
  row.names    = rownames(X_prot)
)

anno_colors <- list(
  Time_point   = c(T1 = "steelblue", T2 = "firebrick"),
  WHO_severity = colorRampPalette(c("#F7FBFF", "#08306B"))(7)  # light → dark blue
)
names(anno_colors$WHO_severity) <- as.character(1:7)

# ── 7. Column annotation — |rho| with who_delta ───────────────────────────────
# Re-apply the PVE ordering for S2 before extracting rho values.

make_col_anno <- function(assoc_df, variant_str, col_ord = NULL) {
  df <- assoc_df[assoc_df$variant == variant_str, ]
  df <- df[order(df$factor), ]   # ensure factor order 1…K
  rho_vals <- abs(df$spearman_rho)
  rho_vals[is.na(rho_vals)] <- 0
  if (!is.null(col_ord)) rho_vals <- rho_vals[col_ord]
  data.frame(
    abs_rho_who_delta = round(rho_vals, 3),
    row.names         = paste0("F", seq_along(rho_vals))
  )
}

anno_col_S1 <- make_col_anno(assoc, "S1_NMF")
anno_col_S2 <- make_col_anno(assoc, "S2_EBNMF", col_ord = S2_col_ord)
anno_col_S3 <- make_col_anno(assoc, "S3_GBCD")

# Shared column annotation color: white → orange gradient
rho_col_fn <- colorRampPalette(c("white", "#E6550D"))(100)

anno_col_colors <- list(
  abs_rho_who_delta = rho_col_fn
)

# ── 8. Heatmap color scale ────────────────────────────────────────────────────

heat_cols <- colorRampPalette(c("gray96", "red"))(50)
heat_brks <- seq(0, 1, length.out = 51)

# ── 9. Heatmap function ───────────────────────────────────────────────────────

draw_heatmap <- function(L_mat, row_ord, anno_row, anno_col, anno_colors_row,
                         anno_col_colors, title_str, filename,
                         width = 7, height = 9) {
  pdf(filename, width = width, height = height)
  pheatmap(
    L_mat[row_ord, ],
    cluster_rows        = FALSE,
    cluster_cols        = FALSE,
    show_rownames       = FALSE,
    annotation_row      = anno_row,
    annotation_col      = anno_col,
    annotation_colors   = c(anno_colors_row, anno_col_colors),
    annotation_names_row = FALSE,
    angle_col           = 45,
    fontsize            = 9,
    color               = heat_cols,
    breaks              = heat_brks,
    main                = title_str,
    border_color        = NA
  )
  dev.off()
  cat(sprintf("Saved: %s\n", basename(filename)))
}

# ── 10. Draw individual heatmaps ──────────────────────────────────────────────

draw_heatmap(
  L_mat          = L_S1,
  row_ord        = row_ord,
  anno_row       = anno_row,
  anno_col       = anno_col_S1,
  anno_colors_row = anno_colors,
  anno_col_colors = anno_col_colors,
  title_str      = "S1 Standard NMF (K=8)\nLoadings normalised per column",
  filename       = file.path(fig_dir, "heatmap_S1_loadings.pdf"),
  width = 6, height = 9
)

draw_heatmap(
  L_mat          = L_S2,
  row_ord        = row_ord,
  anno_row       = anno_row,
  anno_col       = anno_col_S2,
  anno_colors_row = anno_colors,
  anno_col_colors = anno_col_colors,
  title_str      = "S2 EBNMF point-exponential (K=14, ordered by PVE)\nLoadings normalised per column",
  filename       = file.path(fig_dir, "heatmap_S2_loadings.pdf"),
  width = 8, height = 9
)

draw_heatmap(
  L_mat          = L_S3,
  row_ord        = row_ord,
  anno_row       = anno_row,
  anno_col       = anno_col_S3,
  anno_colors_row = anno_colors,
  anno_col_colors = anno_col_colors,
  title_str      = "S3 GBCD (K=10)\nLoadings normalised per column",
  filename       = file.path(fig_dir, "heatmap_S3_loadings.pdf"),
  width = 7, height = 9
)

# ── 11. Combined 3-panel PDF ──────────────────────────────────────────────────
# Draw all three heatmaps onto a single wide PDF page using grid capture.

capture_pheatmap <- function(L_mat, row_ord, anno_row, anno_col,
                              anno_colors_row, anno_col_colors, title_str) {
  pheatmap(
    L_mat[row_ord, ],
    cluster_rows        = FALSE,
    cluster_cols        = FALSE,
    show_rownames       = FALSE,
    annotation_row      = anno_row,
    annotation_col      = anno_col,
    annotation_colors   = c(anno_colors_row, anno_col_colors),
    annotation_names_row = FALSE,
    angle_col           = 45,
    fontsize            = 8,
    color               = heat_cols,
    breaks              = heat_brks,
    main                = title_str,
    border_color        = NA,
    silent              = TRUE   # return grob without drawing
  )
}

g1 <- capture_pheatmap(L_S1, row_ord, anno_row, anno_col_S1,
                        anno_colors, anno_col_colors,
                        "S1 Standard NMF (K=8)")
g2 <- capture_pheatmap(L_S2, row_ord, anno_row, anno_col_S2,
                        anno_colors, anno_col_colors,
                        "S2 EBNMF (K=14, by PVE)")
g3 <- capture_pheatmap(L_S3, row_ord, anno_row, anno_col_S3,
                        anno_colors, anno_col_colors,
                        "S3 GBCD (K=10)")

pdf(file.path(fig_dir, "heatmap_all3_loadings.pdf"), width = 22, height = 9)
grid.arrange(g1$gtable, g2$gtable, g3$gtable, ncol = 3)
dev.off()
cat("Saved: heatmap_all3_loadings.pdf\n")

cat("\n--- 04_loading_heatmaps.R complete ---\n")
