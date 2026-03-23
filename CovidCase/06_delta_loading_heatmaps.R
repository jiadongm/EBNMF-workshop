# 06_delta_loading_heatmaps.R — Pheatmap of L matrices for the delta NMF analysis
#
# Mirrors 04_loading_heatmaps.R but for the delta (T2-T1) matrix.
# Four variants: D1 (NMF), D2 (EBMF Laplace), D3 (EBMF Normal), D4 (GBCD).
#
# Rows  = 120 patients (one per patient, ordered by who_delta then T1 severity)
# Cols  = programs (K per variant; D2/D3 ordered by PVE)
# Color = gray96→red for D1 [0,1]; blue→white→red (diverging) for D2/D3/D4 (signed)
#
# Row annotation:
#   WHO_delta   — Improved (<0) / Stable (=0) / Worsened (>0) / Ambiguous (NA)
#   T1_severity — WHO score at admission (T1), continuous 1–7
#
# Column annotation:
#   abs_rho_who_delta — |Spearman rho| with who_delta (white → orange)
#
# Run from CovidCase/:   Rscript 06_delta_loading_heatmaps.R

suppressPackageStartupMessages({
  library(NMF)
  library(flashier)
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

X_prot     <- readRDS(file.path(results_dir, "X_prot.rds"))       # 240 × 481
pca_scores <- readRDS(file.path(results_dir, "pca_scores.rds"))
assoc      <- read.csv(file.path(nmf_dir, "delta_association.csv"))

fit_D1 <- readRDS(file.path(nmf_dir, "delta_D1_nmf.rds"))
fit_D2 <- readRDS(file.path(nmf_dir, "delta_D2_laplace.rds"))
fit_D3 <- readRDS(file.path(nmf_dir, "delta_D3_normal.rds"))
fit_D4 <- readRDS(file.path(nmf_dir, "delta_D4_gbcd.rds"))
fit_D5 <- readRDS(file.path(nmf_dir, "delta_D5_asym_ebmf.rds"))

# ── 2. Reconstruct delta metadata (120 patients) ───────────────────────────────

meta   <- pca_scores[match(rownames(X_prot), pca_scores$sample), ]
T1_idx <- which(meta$time_point == "T1")
T2_idx <- which(meta$time_point == "T2")
pid_T1 <- meta$patient_id[T1_idx]
pid_T2 <- meta$patient_id[T2_idx]
shared_patients <- intersect(pid_T1, pid_T2)    # 120 patients
T1_order <- T1_idx[match(shared_patients, pid_T1)]
meta_delta <- meta[T1_order, ]
rownames(meta_delta) <- shared_patients

who_delta_vec <- meta_delta$who_delta
who_score_T1  <- meta_delta$who_score_num   # severity at admission

# ── 3. Extract and align L matrices ───────────────────────────────────────────

L_D1_raw <- t(NMF::coef(fit_D1))   # 120 × 6  (non-negative)
L_D2_raw <- fit_D2$L_pm             # 120 × 20 (signed Laplace)
L_D3_raw <- fit_D3$L_pm             # 120 × 20 (signed normal)
L_D4_raw <- fit_D4$L                # 120 × 18 (binary-like, slight negatives)
L_D5_raw <- fit_D5$L_pm             # 120 × 20 (non-negative: point-exponential prior)

# Ensure all have correct rownames (patient IDs)
rownames(L_D1_raw) <- shared_patients
rownames(L_D2_raw) <- shared_patients
rownames(L_D3_raw) <- shared_patients
rownames(L_D4_raw) <- shared_patients
rownames(L_D5_raw) <- shared_patients

cat(sprintf("L dims — D1:%dx%d  D2:%dx%d  D3:%dx%d  D4:%dx%d  D5:%dx%d\n",
            nrow(L_D1_raw), ncol(L_D1_raw),
            nrow(L_D2_raw), ncol(L_D2_raw),
            nrow(L_D3_raw), ncol(L_D3_raw),
            nrow(L_D4_raw), ncol(L_D4_raw),
            nrow(L_D5_raw), ncol(L_D5_raw)))

# ── 4. Column ordering ────────────────────────────────────────────────────────
# D1, D4: keep as-fit.  D2, D3, D5: sort by PVE (decreasing).

pve_D2     <- fit_D2$pve
pve_D3     <- fit_D3$pve
pve_D5     <- fit_D5$pve
D2_col_ord <- order(pve_D2, decreasing = TRUE)
D3_col_ord <- order(pve_D3, decreasing = TRUE)
D5_col_ord <- order(pve_D5, decreasing = TRUE)

L_D2 <- L_D2_raw[, D2_col_ord, drop = FALSE]
L_D3 <- L_D3_raw[, D3_col_ord, drop = FALSE]
L_D5 <- L_D5_raw[, D5_col_ord, drop = FALSE]

colnames(L_D1_raw) <- paste0("F", seq_len(ncol(L_D1_raw)))
colnames(L_D2)     <- paste0("F", seq_len(ncol(L_D2)))
colnames(L_D3)     <- paste0("F", seq_len(ncol(L_D3)))
colnames(L_D4_raw) <- paste0("F", seq_len(ncol(L_D4_raw)))
colnames(L_D5)     <- paste0("F", seq_len(ncol(L_D5)))

# ── 5. Row ordering ───────────────────────────────────────────────────────────
# Primary: ascending who_delta (improved first, NA patients last)
# Secondary: ascending T1 severity within each who_delta level

row_ord <- order(
  is.na(who_delta_vec),          # 0 = non-NA first, 1 = NA (ambiguous) last
  who_delta_vec,                 # ascending: -4, -3, -2, -1, 0, +1, +2
  who_score_T1,                  # then by T1 severity
  na.last = TRUE
)

# ── 6. Row annotation ─────────────────────────────────────────────────────────

who_delta_cat <- ifelse(
  is.na(who_delta_vec), "Ambiguous",
  ifelse(who_delta_vec < 0, "Improved",
  ifelse(who_delta_vec == 0, "Stable", "Worsened"))
)

# T1 severity: treat as character for discrete pheatmap colors
sev_t1_chr    <- as.character(who_score_T1)
sev_t1_levels <- sort(unique(sev_t1_chr))   # "1","1.5","3","4","5","6","7"
sev_t1_pal    <- colorRampPalette(c("#F7FBFF", "#08306B"))(length(sev_t1_levels))
names(sev_t1_pal) <- sev_t1_levels

anno_row <- data.frame(
  WHO_delta   = who_delta_cat,
  T1_severity = sev_t1_chr,
  row.names   = shared_patients
)

anno_colors <- list(
  WHO_delta   = c(Ambiguous = "gray70",
                  Improved  = "steelblue",
                  Stable    = "gold3",
                  Worsened  = "firebrick"),
  T1_severity = sev_t1_pal
)

# ── 7. Column annotation — |rho| with who_delta ───────────────────────────────

rho_col_fn      <- colorRampPalette(c("white", "#E6550D"))(100)
anno_col_colors <- list(abs_rho_who_delta = rho_col_fn)

make_col_anno_delta <- function(variant_str, col_ord = NULL) {
  df  <- assoc[assoc$variant == variant_str, ]
  df  <- df[order(df$factor), ]
  rho <- abs(df$spearman_rho)
  rho[is.na(rho)] <- 0
  if (!is.null(col_ord)) rho <- rho[col_ord]
  data.frame(abs_rho_who_delta = round(rho, 3),
             row.names         = paste0("F", seq_along(rho)))
}

anno_col_D1 <- make_col_anno_delta("D1_NMF")
anno_col_D2 <- make_col_anno_delta("D2_EBMF_Laplace", col_ord = D2_col_ord)
anno_col_D3 <- make_col_anno_delta("D3_EBMF_Normal",  col_ord = D3_col_ord)
anno_col_D4 <- make_col_anno_delta("D4_GBCD")
anno_col_D5 <- make_col_anno_delta("D5_Asym_EBMF",    col_ord = D5_col_ord)

# ── 8. Color scales ───────────────────────────────────────────────────────────

# D1: non-negative (col-normalised [0,1])
col_norm <- function(M) {
  cm <- apply(M, 2, max)
  cm[cm == 0] <- 1
  sweep(M, 2, cm, "/")
}
heat_cols_pos  <- colorRampPalette(c("gray96", "red"))(50)
heat_brks_pos  <- seq(0, 1, length.out = 51)

# D2/D3/D4: signed, diverging (computed per matrix)
div_scale <- function(M) {
  lim  <- ceiling(max(abs(M), na.rm = TRUE) * 10) / 10
  cols <- colorRampPalette(c("steelblue", "white", "red"))(101)
  brks <- seq(-lim, lim, length.out = 102)
  list(cols = cols, brks = brks)
}

# ── 9. Heatmap function ───────────────────────────────────────────────────────

draw_delta_heatmap <- function(L_mat, row_ord, anno_row, anno_col,
                                anno_colors_all, heat_cols, heat_brks,
                                title_str, filename, width = 8, height = 8) {
  pdf(filename, width = width, height = height)
  pheatmap(
    L_mat[row_ord, ],
    cluster_rows        = FALSE,
    cluster_cols        = FALSE,
    show_rownames       = FALSE,
    annotation_row      = anno_row,
    annotation_col      = anno_col,
    annotation_colors   = anno_colors_all,
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

# ── 10. Individual heatmaps ───────────────────────────────────────────────────

# D1 — Standard NMF (non-negative, col-normalised)
L_D1_n <- col_norm(L_D1_raw)
draw_delta_heatmap(
  L_mat          = L_D1_n,
  row_ord        = row_ord,
  anno_row       = anno_row,
  anno_col       = anno_col_D1,
  anno_colors_all = c(anno_colors, anno_col_colors),
  heat_cols      = heat_cols_pos,
  heat_brks      = heat_brks_pos,
  title_str      = "D1 Standard NMF (K=6, shifted delta)\nLoadings col-normalised [0,1]",
  filename       = file.path(fig_dir, "heatmap_delta_D1_loadings.pdf"),
  width = 6, height = 8
)

# D2 — EBMF Laplace (signed, PVE order)
sc_D2 <- div_scale(L_D2)
draw_delta_heatmap(
  L_mat          = L_D2,
  row_ord        = row_ord,
  anno_row       = anno_row,
  anno_col       = anno_col_D2,
  anno_colors_all = c(anno_colors, anno_col_colors),
  heat_cols      = sc_D2$cols,
  heat_brks      = sc_D2$brks,
  title_str      = "D2 EBMF point-Laplace (K=20, PVE order)\nSigned loadings (blue=neg, red=pos)",
  filename       = file.path(fig_dir, "heatmap_delta_D2_loadings.pdf"),
  width = 10, height = 8
)

# D3 — EBMF Normal (signed, PVE order)
sc_D3 <- div_scale(L_D3)
draw_delta_heatmap(
  L_mat          = L_D3,
  row_ord        = row_ord,
  anno_row       = anno_row,
  anno_col       = anno_col_D3,
  anno_colors_all = c(anno_colors, anno_col_colors),
  heat_cols      = sc_D3$cols,
  heat_brks      = sc_D3$brks,
  title_str      = "D3 EBMF point-normal (K=20, PVE order)\nSigned loadings (blue=neg, red=pos)",
  filename       = file.path(fig_dir, "heatmap_delta_D3_loadings.pdf"),
  width = 10, height = 8
)

# D4 — GBCD (binary-like, slight negatives)
sc_D4 <- div_scale(L_D4_raw)
draw_delta_heatmap(
  L_mat          = L_D4_raw,
  row_ord        = row_ord,
  anno_row       = anno_row,
  anno_col       = anno_col_D4,
  anno_colors_all = c(anno_colors, anno_col_colors),
  heat_cols      = sc_D4$cols,
  heat_brks      = sc_D4$brks,
  title_str      = "D4 GBCD (K=18, binary-like)\nLoadings (mostly [0,1])",
  filename       = file.path(fig_dir, "heatmap_delta_D4_loadings.pdf"),
  width = 9, height = 8
)

# D5 — Asymmetric EBMF (non-negative L: col-normalised [0,1], gray96→red)
L_D5_n <- col_norm(L_D5)
draw_delta_heatmap(
  L_mat          = L_D5_n,
  row_ord        = row_ord,
  anno_row       = anno_row,
  anno_col       = anno_col_D5,
  anno_colors_all = c(anno_colors, anno_col_colors),
  heat_cols      = heat_cols_pos,
  heat_brks      = heat_brks_pos,
  title_str      = "D5 Asymmetric EBMF (K=20, PVE order)\nLoadings col-normalised [0,1] (non-negative)",
  filename       = file.path(fig_dir, "heatmap_delta_D5_loadings.pdf"),
  width = 10, height = 8
)

# ── 11. Combined figure — all five variants ────────────────────────────────────
# Capture each as a grob and arrange side-by-side (D1|D2|D3|D4|D5)

capture_pheatmap_delta <- function(L_mat, row_ord, anno_row, anno_col,
                                    anno_colors_all, heat_cols, heat_brks,
                                    title_str) {
  pheatmap(
    L_mat[row_ord, ],
    cluster_rows        = FALSE,
    cluster_cols        = FALSE,
    show_rownames       = FALSE,
    annotation_row      = anno_row,
    annotation_col      = anno_col,
    annotation_colors   = anno_colors_all,
    annotation_names_row = FALSE,
    angle_col           = 45,
    fontsize            = 7,
    color               = heat_cols,
    breaks              = heat_brks,
    main                = title_str,
    border_color        = NA,
    silent              = TRUE
  )
}

g1 <- capture_pheatmap_delta(L_D1_n,    row_ord, anno_row, anno_col_D1,
                              c(anno_colors, anno_col_colors),
                              heat_cols_pos, heat_brks_pos, "D1 NMF (K=6)")
g2 <- capture_pheatmap_delta(L_D2,      row_ord, anno_row, anno_col_D2,
                              c(anno_colors, anno_col_colors),
                              sc_D2$cols, sc_D2$brks, "D2 EBMF Laplace (K=20, PVE)")
g3 <- capture_pheatmap_delta(L_D3,      row_ord, anno_row, anno_col_D3,
                              c(anno_colors, anno_col_colors),
                              sc_D3$cols, sc_D3$brks, "D3 EBMF Normal (K=20, PVE)")
g4 <- capture_pheatmap_delta(L_D4_raw,  row_ord, anno_row, anno_col_D4,
                              c(anno_colors, anno_col_colors),
                              sc_D4$cols, sc_D4$brks, "D4 GBCD (K=18)")
g5 <- capture_pheatmap_delta(L_D5_n,    row_ord, anno_row, anno_col_D5,
                              c(anno_colors, anno_col_colors),
                              heat_cols_pos, heat_brks_pos, "D5 Asym EBMF (K=20, PVE)")

pdf(file.path(fig_dir, "heatmap_delta_all5_loadings.pdf"), width = 30, height = 9)
grid.arrange(g1$gtable, g2$gtable, g3$gtable, g4$gtable, g5$gtable, ncol = 5)
dev.off()
cat("Saved: heatmap_delta_all5_loadings.pdf\n")

cat("\n--- 06_delta_loading_heatmaps.R complete ---\n")
