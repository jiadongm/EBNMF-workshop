# 03_NMF_delta.R — NMF on T2-T1 difference matrix (PRIMARY)
#
# Biological question:
#   What protein programs CHANGED between T1 and T2?
#   Which programs are linked to WHO severity change (who_delta = T2 - T1)?
#
# Matrix: delta_mat (120 × 481), one row per patient, values = (T2 NPX - T1 NPX).
#   Range: approximately [-6.19, 7.28] — signed.
#
# Five variants:
#   D1 — Standard NMF on min-shifted delta (removes directionality, baseline)
#   D2 — EBMF, point-Laplace prior on BOTH L and F (signed sparse)
#   D3 — EBMF, point-normal prior on BOTH L and F  (signed dense — comparison)
#   D4 — GBCD (generalised binary, binary-like patient memberships + signed LFC)
#   D5 — Asymmetric EBMF: point-exponential on L (non-neg), point-Laplace on F (signed)
#        Direction of protein change → F (signed); patient participation → L (≥ 0)
#        Best of both: non-neg L enables all flashier plots; signed F captures direction
#
# Outputs:
#   results/nmf/delta_D1_nmf.rds         — NMFfit
#   results/nmf/delta_D2_laplace.rds     — flashier fit (point-Laplace, both)
#   results/nmf/delta_D3_normal.rds      — flashier fit (point-normal, both)
#   results/nmf/delta_D4_gbcd.rds        — gbcd fit (L, F$lfc, F$lfsr)
#   results/nmf/delta_D5_asym_ebmf.rds  — flashier fit (asymmetric: exp-L / Laplace-F)
#   results/nmf/delta_association.csv    — Spearman corr(L_k, who_delta) per factor
#   results/figs/nmf/delta_*.pdf         — figures

suppressPackageStartupMessages({
  library(NMF)
  library(flashier)
  library(ebnm)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# gbcd: install from GitHub if not present (first run only; updates flashier)
if (!requireNamespace("gbcd", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("stephenslab/gbcd")
}
suppressPackageStartupMessages(library(gbcd))

# ── Paths ──────────────────────────────────────────────────────────────────────

results_dir <- "results"
nmf_dir     <- file.path(results_dir, "nmf")
fig_dir     <- file.path(results_dir, "figs", "nmf")
dir.create(nmf_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load data ───────────────────────────────────────────────────────────────

X_prot     <- readRDS(file.path(results_dir, "X_prot.rds"))       # 240 × 481
pca_scores <- readRDS(file.path(results_dir, "pca_scores.rds"))   # PC scores + metadata

cat(sprintf("X_prot dimensions: %d samples × %d proteins\n", nrow(X_prot), ncol(X_prot)))

# ── 2. Build delta matrix (T2 - T1 per patient) ────────────────────────────────

meta <- pca_scores[match(rownames(X_prot), pca_scores$sample), ]
stopifnot(all(meta$sample == rownames(X_prot)))

# Extract T1 and T2 sub-matrices, match by patient_id
T1_idx <- which(meta$time_point == "T1")
T2_idx <- which(meta$time_point == "T2")
pid_T1 <- meta$patient_id[T1_idx]
pid_T2 <- meta$patient_id[T2_idx]
shared_patients <- intersect(pid_T1, pid_T2)
cat(sprintf("Patients with both T1 and T2: %d\n", length(shared_patients)))

T1_order <- T1_idx[match(shared_patients, pid_T1)]
T2_order <- T2_idx[match(shared_patients, pid_T2)]

X_T1 <- X_prot[T1_order, ]   # 120 × 481
X_T2 <- X_prot[T2_order, ]   # 120 × 481

# Delta matrix: T2 - T1 (change in NPX log2 units)
delta_mat <- X_T2 - X_T1     # 120 × 481
rownames(delta_mat) <- shared_patients

cat(sprintf("delta_mat: %d patients × %d proteins\n", nrow(delta_mat), ncol(delta_mat)))
cat(sprintf("delta_mat range: [%.3f, %.3f]\n", min(delta_mat), max(delta_mat)))

# Patient metadata aligned to delta_mat rows (120 patients)
meta_delta <- meta[T1_order, ]   # T1 rows carry who_delta computed in 01_PCA.R
rownames(meta_delta) <- shared_patients
stopifnot(all(meta_delta$patient_id == shared_patients))

who_delta_vec <- meta_delta$who_delta
cat(sprintf("who_delta: %d non-NA; range [%d, %d]; %d stable (=0)\n",
            sum(!is.na(who_delta_vec)),
            min(who_delta_vec, na.rm = TRUE),
            max(who_delta_vec, na.rm = TRUE),
            sum(who_delta_vec == 0, na.rm = TRUE)))

# ── 3. Min-shifted delta for standard NMF (D1) ─────────────────────────────────

prot_min_delta <- apply(delta_mat, 2, min)
delta_shifted  <- sweep(delta_mat, 2, prot_min_delta, "-")  # per-column min-shift
cat(sprintf("delta_shifted range: [%.3f, %.3f]\n", min(delta_shifted), max(delta_shifted)))
stopifnot(min(delta_shifted) >= 0)

# ── Helper: association analysis ───────────────────────────────────────────────
# For delta script: L has 120 rows (one per patient), associate directly with who_delta

compute_associations_delta <- function(L_mat, who_delta_vec, variant_label) {
  K     <- ncol(L_mat)
  non_na <- !is.na(who_delta_vec)
  cat(sprintf("\n--- Association analysis: %s (K=%d, n=%d non-NA who_delta) ---\n",
              variant_label, K, sum(non_na)))

  results <- lapply(seq_len(K), function(k) {
    if (sum(non_na) < 5) return(data.frame(variant = variant_label, factor = k,
                                            spearman_rho = NA, p_value = NA))
    ct <- cor.test(L_mat[non_na, k], who_delta_vec[non_na], method = "spearman",
                   exact = FALSE)
    data.frame(
      variant      = variant_label,
      factor       = k,
      spearman_rho = round(ct$estimate, 4),
      p_value      = signif(ct$p.value, 3)
    )
  })

  assoc_df <- do.call(rbind, results)
  assoc_df$bh_fdr <- p.adjust(assoc_df$p_value, method = "BH")
  assoc_df <- assoc_df[order(abs(assoc_df$spearman_rho), decreasing = TRUE), ]

  top <- assoc_df[!is.na(assoc_df$spearman_rho) & abs(assoc_df$spearman_rho) > 0.2, ]
  cat(sprintf("  Factors with |rho| > 0.2 and BH FDR < 0.2:\n"))
  top_sig <- top[!is.na(top$bh_fdr) & top$bh_fdr < 0.2, ]
  if (nrow(top_sig) > 0) {
    print(top_sig, row.names = FALSE)
  } else {
    cat("  None.\n")
    if (nrow(top) > 0) {
      cat("  All factors with |rho| > 0.2 (any FDR):\n")
      print(top, row.names = FALSE)
    }
  }
  return(assoc_df)
}

# ── Helper: PVE bar chart ──────────────────────────────────────────────────────

plot_pve <- function(pve_vec, title_str) {
  df <- data.frame(factor = seq_along(pve_vec), pve = pve_vec * 100)
  ggplot(df, aes(x = factor(factor), y = pve)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pve)), vjust = -0.3, size = 2.8) +
    labs(title = title_str, x = "Factor", y = "% variance explained") +
    theme_bw(base_size = 10)
}

# ── Helper: factor-protein plots (signed or unsigned) ─────────────────────────

plot_factor_proteins_signed <- function(F_mat, top_n = 25, variant_label = "",
                                         max_factors = 6) {
  K      <- min(ncol(F_mat), max_factors)
  pnames <- rownames(F_mat)
  plots  <- vector("list", K)
  for (k in seq_len(K)) {
    f_k <- F_mat[, k]
    ord <- order(abs(f_k), decreasing = TRUE)[1:min(top_n, length(f_k))]
    df  <- data.frame(
      protein = pnames[ord],
      weight  = f_k[ord],
      sign    = ifelse(f_k[ord] >= 0, "positive", "negative")
    ) |>
      dplyr::arrange(dplyr::desc(weight)) |>
      dplyr::mutate(protein = factor(protein, levels = rev(protein)))

    plots[[k]] <- ggplot(df, aes(x = weight, y = protein, fill = sign)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = c(positive = "steelblue", negative = "firebrick"),
                        guide = "none") +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      labs(title = sprintf("Factor %d", k), x = "Weight", y = NULL) +
      theme_bw(base_size = 8) +
      theme(axis.text.y = element_text(size = 6))
  }
  wrap_plots(plots, ncol = 3) +
    plot_annotation(title = sprintf("%s — top proteins per factor (|weight|)", variant_label))
}

# ── Helper: loading dot plot sorted by L[,k], coloured by who_delta ───────────

plot_loading_dotplot <- function(L_mat, who_delta_vec, variant_label, max_factors = 6) {
  K <- min(ncol(L_mat), max_factors)
  plots <- vector("list", K)
  for (k in seq_len(K)) {
    df <- data.frame(
      patient   = seq_len(nrow(L_mat)),
      score     = L_mat[, k],
      who_delta = who_delta_vec
    ) |>
      dplyr::arrange(score) |>
      dplyr::mutate(rank = seq_len(nrow(L_mat)))

    plots[[k]] <- ggplot(df, aes(x = rank, y = score, colour = who_delta)) +
      geom_point(size = 1.8, alpha = 0.85) +
      scale_colour_gradient2(low = "steelblue", mid = "gold", high = "firebrick",
                             midpoint = 0, name = "dWHO", na.value = "grey70") +
      labs(title = sprintf("Factor %d", k),
           x = "Patient (sorted by score)", y = "Loading") +
      theme_bw(base_size = 9)
  }
  plots <- Filter(Negate(is.null), plots)
  wrap_plots(plots, ncol = 3) +
    plot_annotation(title = sprintf("%s — loading dot plot (sorted by score, coloured by who_delta)",
                                    variant_label))
}

# ══════════════════════════════════════════════════════════════════════════════
# VARIANT D1 — Standard NMF on min-shifted delta (K=6)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== D1: Standard NMF on min-shifted delta (K=6) ==========\n")

K_D1 <- 6

if (!file.exists(file.path(nmf_dir, "delta_D1_nmf.rds"))) {
  set.seed(123)
  t_start <- proc.time()
  fit_D1 <- NMF::nmf(
    x      = t(delta_shifted),   # NMF package: features × samples (481 × 120)
    rank   = K_D1,
    method = "brunet",
    nrun   = 10,
    seed   = 123
  )
  t_D1 <- proc.time() - t_start
  saveRDS(fit_D1, file.path(nmf_dir, "delta_D1_nmf.rds"))
  cat(sprintf("D1 NMF fitted and saved.  [%.1f s]\n", t_D1["elapsed"]))
} else {
  fit_D1 <- readRDS(file.path(nmf_dir, "delta_D1_nmf.rds"))
  cat("D1 NMF loaded from cache.\n")
}

W_D1 <- NMF::basis(fit_D1)   # 481 × 6
H_D1 <- NMF::coef(fit_D1)    # 6 × 120
L_D1 <- t(H_D1)              # 120 × 6
F_D1 <- W_D1                  # 481 × 6

cat(sprintf("D1: L %d×%d, F %d×%d\n", nrow(L_D1), ncol(L_D1), nrow(F_D1), ncol(F_D1)))

assoc_D1 <- compute_associations_delta(L_D1, who_delta_vec, "D1_NMF")

p_D1_factors <- plot_factor_proteins_signed(F_D1, variant_label = "D1 Standard NMF (shifted delta)")
ggsave(file.path(fig_dir, "delta_D1_factor_proteins.pdf"),
       p_D1_factors, width = 14, height = 10)
cat("Saved: delta_D1_factor_proteins.pdf\n")

p_D1_dots <- plot_loading_dotplot(L_D1, who_delta_vec, "D1_NMF")
ggsave(file.path(fig_dir, "delta_D1_loading_dots.pdf"),
       p_D1_dots, width = 12, height = 8)
cat("Saved: delta_D1_loading_dots.pdf\n")

# ══════════════════════════════════════════════════════════════════════════════
# VARIANT D2 — EBMF, point-Laplace (signed sparse, PRIMARY)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== D2: EBMF point-Laplace (signed sparse, PRIMARY) ==========\n")

if (!file.exists(file.path(nmf_dir, "delta_D2_laplace.rds"))) {
  set.seed(42)
  t_start <- proc.time()
  fit_D2 <- flash(
    data        = delta_mat,          # 120 × 481 (signed, no shift)
    ebnm_fn     = ebnm_point_laplace, # symmetric, sparse
    greedy_Kmax = 20,
    backfit     = TRUE,
    nullcheck   = TRUE,
    verbose     = 1L
  )
  t_D2 <- proc.time() - t_start
  saveRDS(fit_D2, file.path(nmf_dir, "delta_D2_laplace.rds"))
  cat(sprintf("D2 EBMF (point-Laplace) fitted and saved.  [%.1f s]\n", t_D2["elapsed"]))
} else {
  fit_D2 <- readRDS(file.path(nmf_dir, "delta_D2_laplace.rds"))
  cat("D2 EBMF (point-Laplace) loaded from cache.\n")
}

K_D2   <- fit_D2$n_factors
L_D2   <- fit_D2$L_pm    # 120 × K
F_D2   <- fit_D2$F_pm    # 481 × K
pve_D2 <- fit_D2$pve

cat(sprintf("D2 auto K=%d; PVE sum=%.3f\n", K_D2, sum(pve_D2)))
cat(sprintf("D2 PVE per factor: %s\n", paste(round(pve_D2 * 100, 1), collapse = ", ")))

assoc_D2 <- compute_associations_delta(L_D2, who_delta_vec, "D2_EBMF_Laplace")

# PVE plot
p_D2_pve <- plot_pve(pve_D2, "D2 EBMF point-Laplace — PVE per factor")
ggsave(file.path(fig_dir, "delta_D2_pve.pdf"), p_D2_pve, width = 7, height = 4)
cat("Saved: delta_D2_pve.pdf\n")

# Factor-protein plots (signed weights)
p_D2_factors <- plot_factor_proteins_signed(F_D2, variant_label = "D2 EBMF point-Laplace",
                                              max_factors = K_D2)
ggsave(file.path(fig_dir, "delta_D2_factor_proteins.pdf"),
       p_D2_factors, width = 14, height = 10)
cat("Saved: delta_D2_factor_proteins.pdf\n")

# Loading dot plots (PRIMARY — coloured by who_delta)
p_D2_dots <- plot_loading_dotplot(L_D2, who_delta_vec, "D2_EBMF_Laplace (PRIMARY)",
                                   max_factors = K_D2)
ggsave(file.path(fig_dir, "delta_D2_loading_dots.pdf"),
       p_D2_dots, width = 12, height = 8)
cat("Saved: delta_D2_loading_dots.pdf\n")

# ══════════════════════════════════════════════════════════════════════════════
# VARIANT D3 — EBMF, point-normal (signed dense)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== D3: EBMF point-normal (signed dense) ==========\n")

if (!file.exists(file.path(nmf_dir, "delta_D3_normal.rds"))) {
  set.seed(77)
  t_start <- proc.time()
  fit_D3 <- flash(
    data        = delta_mat,         # 120 × 481 (signed)
    ebnm_fn     = ebnm_point_normal, # symmetric, less sparse
    greedy_Kmax = 20,
    backfit     = TRUE,
    nullcheck   = TRUE,
    verbose     = 1L
  )
  t_D3 <- proc.time() - t_start
  saveRDS(fit_D3, file.path(nmf_dir, "delta_D3_normal.rds"))
  cat(sprintf("D3 EBMF (point-normal) fitted and saved.  [%.1f s]\n", t_D3["elapsed"]))
} else {
  fit_D3 <- readRDS(file.path(nmf_dir, "delta_D3_normal.rds"))
  cat("D3 EBMF (point-normal) loaded from cache.\n")
}

K_D3   <- fit_D3$n_factors
L_D3   <- fit_D3$L_pm    # 120 × K
F_D3   <- fit_D3$F_pm    # 481 × K
pve_D3 <- fit_D3$pve

cat(sprintf("D3 auto K=%d; PVE sum=%.3f\n", K_D3, sum(pve_D3)))
cat(sprintf("D3 PVE per factor: %s\n", paste(round(pve_D3 * 100, 1), collapse = ", ")))

assoc_D3 <- compute_associations_delta(L_D3, who_delta_vec, "D3_EBMF_Normal")

# PVE plot
p_D3_pve <- plot_pve(pve_D3, "D3 EBMF point-normal — PVE per factor")
ggsave(file.path(fig_dir, "delta_D3_pve.pdf"), p_D3_pve, width = 7, height = 4)
cat("Saved: delta_D3_pve.pdf\n")

# Factor-protein plots (signed weights)
p_D3_factors <- plot_factor_proteins_signed(F_D3, variant_label = "D3 EBMF point-normal",
                                              max_factors = K_D3)
ggsave(file.path(fig_dir, "delta_D3_factor_proteins.pdf"),
       p_D3_factors, width = 14, height = 10)
cat("Saved: delta_D3_factor_proteins.pdf\n")

# Loading dot plots
p_D3_dots <- plot_loading_dotplot(L_D3, who_delta_vec, "D3_EBMF_Normal",
                                   max_factors = K_D3)
ggsave(file.path(fig_dir, "delta_D3_loading_dots.pdf"),
       p_D3_dots, width = 12, height = 8)
cat("Saved: delta_D3_loading_dots.pdf\n")

# ══════════════════════════════════════════════════════════════════════════════
# VARIANT D4 — GBCD (generalised binary, binary-like patient memberships)
# ══════════════════════════════════════════════════════════════════════════════
#
# Input: delta_mat (120 × 481, signed T2-T1 differences in log2 NPX units).
# GBCD fits X ≈ L F^T where L has binary-like (0/1) patient memberships and
# F$lfc gives signed log-fold change of each protein per program.
# The generalised binary prior enforces "in/out" membership structure, which
# is well-suited for identifying patient subgroups with distinct change profiles.

cat("\n========== D4: GBCD on delta matrix (binary-like memberships) ==========\n")

if (!file.exists(file.path(nmf_dir, "delta_D4_gbcd.rds"))) {
  set.seed(99)
  t_start <- proc.time()
  fit_D4 <- fit_gbcd(
    Y        = delta_mat,                      # 120 × 481, signed
    Kmax     = 10,
    prior    = ebnm::ebnm_generalized_binary,
    maxiter1 = 500, maxiter2 = 200, maxiter3 = 500,
    verbose  = 1
  )
  t_D4 <- proc.time() - t_start
  saveRDS(fit_D4, file.path(nmf_dir, "delta_D4_gbcd.rds"))
  cat(sprintf("D4 GBCD fitted and saved.  [%.1f s]\n", t_D4["elapsed"]))
} else {
  fit_D4 <- readRDS(file.path(nmf_dir, "delta_D4_gbcd.rds"))
  cat("D4 GBCD loaded from cache.\n")
}

K_D4   <- ncol(fit_D4$L)
L_D4   <- fit_D4$L           # 120 × K, binary-like memberships in [0,1]
F_D4_lfc  <- fit_D4$F$lfc   # 481 × K, signed LFC per protein per program
F_D4_lfsr <- fit_D4$F$lfsr  # 481 × K, local false sign rate

# Set rownames if not already present (gbcd does set them from the input matrix)
if (is.null(rownames(L_D4)))   rownames(L_D4)      <- rownames(delta_mat)
if (is.null(rownames(F_D4_lfc))) rownames(F_D4_lfc) <- colnames(delta_mat)

cat(sprintf("D4 auto K=%d\n", K_D4))
cat(sprintf("D4 L range: [%.3f, %.3f]  (expected ~ [0,1])\n",
            min(L_D4), max(L_D4)))

assoc_D4 <- compute_associations_delta(L_D4, who_delta_vec, "D4_GBCD")

# ── Factor-protein LFC bar chart (lfsr shading) ───────────────────────────────
# Proteins with lfsr >= 0.05 are shown semi-transparent (grayed signal).

plot_factor_lfc_gbcd <- function(F_lfc, F_lfsr, top_n = 25, variant_label = "",
                                  max_factors = 6, lfsr_thresh = 0.05) {
  K      <- min(ncol(F_lfc), max_factors)
  pnames <- rownames(F_lfc)
  plots  <- vector("list", K)
  for (k in seq_len(K)) {
    lfc_k  <- F_lfc[, k]
    lfsr_k <- F_lfsr[, k]
    ord    <- order(abs(lfc_k), decreasing = TRUE)[1:min(top_n, length(lfc_k))]
    df     <- data.frame(
      protein = pnames[ord],
      lfc     = lfc_k[ord],
      lfsr    = lfsr_k[ord],
      sign    = ifelse(lfc_k[ord] >= 0, "positive", "negative"),
      sig     = lfsr_k[ord] < lfsr_thresh
    ) |>
      dplyr::arrange(dplyr::desc(lfc)) |>
      dplyr::mutate(protein = factor(protein, levels = rev(protein)))

    plots[[k]] <- ggplot(df, aes(x = lfc, y = protein, fill = sign,
                                  alpha = ifelse(sig, 1, 0.25))) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = c(positive = "steelblue", negative = "firebrick"),
                        guide = "none") +
      scale_alpha_identity() +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      labs(title = sprintf("Factor %d", k), x = "LFC", y = NULL) +
      theme_bw(base_size = 8) +
      theme(axis.text.y = element_text(size = 6))
  }
  wrap_plots(plots, ncol = 3) +
    plot_annotation(title = sprintf("%s — top proteins per factor (faded = lfsr >= %.2f)",
                                    variant_label, lfsr_thresh))
}

p_D4_factors <- plot_factor_lfc_gbcd(F_D4_lfc, F_D4_lfsr,
                                      variant_label = "D4 GBCD (delta matrix)",
                                      max_factors = K_D4)
ggsave(file.path(fig_dir, "delta_D4_factor_proteins.pdf"),
       p_D4_factors, width = 14, height = 10)
cat("Saved: delta_D4_factor_proteins.pdf\n")

p_D4_dots <- plot_loading_dotplot(L_D4, who_delta_vec, "D4_GBCD", max_factors = K_D4)
ggsave(file.path(fig_dir, "delta_D4_loading_dots.pdf"),
       p_D4_dots, width = 12, height = 8)
cat("Saved: delta_D4_loading_dots.pdf\n")

# ══════════════════════════════════════════════════════════════════════════════
# VARIANT D5 — Asymmetric EBMF: point-exponential on L, point-Laplace on F
# ══════════════════════════════════════════════════════════════════════════════
#
# ebnm_fn as a list: list(fn_for_L, fn_for_F)
#   L prior: ebnm_point_exponential — non-negative patient participation scores
#   F prior: ebnm_point_laplace     — signed sparse protein weights
#
# Biological interpretation:
#   L[patient, k] >= 0 — how strongly patient i participated in change program k
#   F[protein, k] signed — which proteins increased (>0) or decreased (<0) in program k
#
# This routes directionality entirely into F and keeps L ≥ 0, enabling all
# flashier plot types (structure, bar, heatmap) that require non-negative L.

cat("\n========== D5: Asymmetric EBMF (exp-L / Laplace-F) ==========\n")

if (!file.exists(file.path(nmf_dir, "delta_D5_asym_ebmf.rds"))) {
  set.seed(55)
  t_start <- proc.time()
  fit_D5 <- flash(
    data        = delta_mat,                           # 120 × 481 (signed, no shift)
    ebnm_fn     = list(ebnm_point_exponential,         # L: non-negative
                       ebnm_point_laplace),            # F: signed sparse
    greedy_Kmax = 20,
    backfit     = TRUE,
    nullcheck   = TRUE,
    verbose     = 1L
  )
  t_D5 <- proc.time() - t_start
  saveRDS(fit_D5, file.path(nmf_dir, "delta_D5_asym_ebmf.rds"))
  cat(sprintf("D5 asymmetric EBMF fitted and saved.  [%.1f s]\n", t_D5["elapsed"]))
} else {
  fit_D5 <- readRDS(file.path(nmf_dir, "delta_D5_asym_ebmf.rds"))
  cat("D5 asymmetric EBMF loaded from cache.\n")
}

K_D5   <- fit_D5$n_factors
L_D5   <- fit_D5$L_pm    # 120 × K  (non-negative)
F_D5   <- fit_D5$F_pm    # 481 × K  (signed)
pve_D5 <- fit_D5$pve

cat(sprintf("D5 auto K=%d; PVE sum=%.3f\n", K_D5, sum(pve_D5)))
cat(sprintf("D5 L range: [%.3f, %.3f]  (expected >= 0)\n", min(L_D5), max(L_D5)))
cat(sprintf("D5 F range: [%.3f, %.3f]  (signed)\n",        min(F_D5), max(F_D5)))
cat(sprintf("D5 PVE per factor: %s\n", paste(round(pve_D5 * 100, 1), collapse = ", ")))

assoc_D5 <- compute_associations_delta(L_D5, who_delta_vec, "D5_Asym_EBMF")

# PVE plot
p_D5_pve <- plot_pve(pve_D5, "D5 Asymmetric EBMF (exp-L / Laplace-F) — PVE per factor")
ggsave(file.path(fig_dir, "delta_D5_pve.pdf"), p_D5_pve, width = 7, height = 4)
cat("Saved: delta_D5_pve.pdf\n")

# Factor-protein plots (signed F weights)
p_D5_factors <- plot_factor_proteins_signed(F_D5, variant_label = "D5 Asymmetric EBMF",
                                              max_factors = K_D5)
ggsave(file.path(fig_dir, "delta_D5_factor_proteins.pdf"),
       p_D5_factors, width = 14, height = 10)
cat("Saved: delta_D5_factor_proteins.pdf\n")

# Loading dot plots (non-negative L, coloured by who_delta)
p_D5_dots <- plot_loading_dotplot(L_D5, who_delta_vec, "D5_Asym_EBMF", max_factors = K_D5)
ggsave(file.path(fig_dir, "delta_D5_loading_dots.pdf"),
       p_D5_dots, width = 12, height = 8)
cat("Saved: delta_D5_loading_dots.pdf\n")

# ══════════════════════════════════════════════════════════════════════════════
# D2 vs D3 comparison — do sparse and dense priors find different patterns?
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== D2 vs D3 comparison ==========\n")

# PVE comparison on the same axes
pve_compare_df <- rbind(
  data.frame(variant = "D2 point-Laplace (sparse)",
             factor = seq_len(K_D2), pve = pve_D2 * 100),
  data.frame(variant = "D3 point-normal (dense)",
             factor = seq_len(K_D3), pve = pve_D3 * 100)
)

p_pve_compare <- ggplot(pve_compare_df,
                         aes(x = factor, y = pve, fill = variant)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("D2 point-Laplace (sparse)" = "#F28E2B",
                                "D3 point-normal (dense)"  = "#76B7B2")) +
  labs(title = "Delta matrix — PVE comparison: D2 (Laplace) vs D3 (normal)",
       x = "Factor", y = "% variance explained", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "delta_D2D3_pve_comparison.pdf"),
       p_pve_compare, width = 9, height = 5)
cat("Saved: delta_D2D3_pve_comparison.pdf\n")

# ══════════════════════════════════════════════════════════════════════════════
# Combined association table + summary
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== Cross-variant summary (delta matrix) ==========\n")

summary_df <- data.frame(
  variant       = c("D1_NMF", "D2_EBMF_Laplace", "D3_EBMF_Normal",
                    "D4_GBCD", "D5_Asym_EBMF"),
  K             = c(K_D1, K_D2, K_D3, K_D4, K_D5),
  K_selection   = c("fixed (6)", "auto (ELBO)", "auto (ELBO)",
                    "auto (GBCD)", "auto (ELBO)"),
  prior_L       = c("non-neg (shifted)", "point-Laplace", "point-normal",
                    "gen. binary",       "point-exponential"),
  prior_F       = c("non-neg (shifted)", "point-Laplace", "point-normal",
                    "signed LFC",        "point-Laplace"),
  top_rho       = c(
    max(abs(assoc_D1$spearman_rho), na.rm = TRUE),
    max(abs(assoc_D2$spearman_rho), na.rm = TRUE),
    max(abs(assoc_D3$spearman_rho), na.rm = TRUE),
    max(abs(assoc_D4$spearman_rho), na.rm = TRUE),
    max(abs(assoc_D5$spearman_rho), na.rm = TRUE)
  )
)
cat("\nVariant summary:\n")
print(summary_df, row.names = FALSE)

# Combined association table
all_assoc <- rbind(assoc_D1, assoc_D2, assoc_D3, assoc_D4, assoc_D5)
write.csv(all_assoc,
          file.path(nmf_dir, "delta_association.csv"),
          row.names = FALSE)
cat("\nSaved: delta_association.csv\n")

# Waterfall plot of rho values for D2 (primary)
plot_assoc_waterfall <- function(assoc_df, variant_label, rho_thresh = 0.2, fdr_thresh = 0.2) {
  df <- assoc_df[!is.na(assoc_df$spearman_rho), ]
  df$factor  <- factor(df$factor, levels = df$factor[order(df$spearman_rho)])
  df$sig     <- df$bh_fdr < fdr_thresh
  ggplot(df, aes(x = factor, y = spearman_rho,
                  fill = ifelse(spearman_rho > 0, "positive", "negative"),
                  alpha = ifelse(sig, 1, 0.45))) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(positive = "firebrick", negative = "steelblue"),
                      guide = "none") +
    scale_alpha_identity() +
    geom_hline(yintercept = c(-rho_thresh, rho_thresh), linetype = "dashed",
               colour = "grey40") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    labs(title = sprintf("%s — Spearman rho(L_k, who_delta)", variant_label),
         subtitle = sprintf("Faded bars: BH FDR ≥ %.1f; dashed: |rho|=%.1f",
                            fdr_thresh, rho_thresh),
         x = "Factor (sorted by rho)", y = "Spearman rho") +
    theme_bw(base_size = 10)
}

p_D2_assoc <- plot_assoc_waterfall(assoc_D2, "D2 EBMF point-Laplace (PRIMARY)")
ggsave(file.path(fig_dir, "delta_D2_association_waterfall.pdf"),
       p_D2_assoc, width = 7, height = 4)
cat("Saved: delta_D2_association_waterfall.pdf\n")

p_D3_assoc <- plot_assoc_waterfall(assoc_D3, "D3 EBMF point-normal")
ggsave(file.path(fig_dir, "delta_D3_association_waterfall.pdf"),
       p_D3_assoc, width = 7, height = 4)
cat("Saved: delta_D3_association_waterfall.pdf\n")

p_D4_assoc <- plot_assoc_waterfall(assoc_D4, "D4 GBCD (binary, delta matrix)")
ggsave(file.path(fig_dir, "delta_D4_association_waterfall.pdf"),
       p_D4_assoc, width = 7, height = 4)
cat("Saved: delta_D4_association_waterfall.pdf\n")

p_D5_assoc <- plot_assoc_waterfall(assoc_D5, "D5 Asymmetric EBMF (exp-L / Laplace-F)")
ggsave(file.path(fig_dir, "delta_D5_association_waterfall.pdf"),
       p_D5_assoc, width = 7, height = 4)
cat("Saved: delta_D5_association_waterfall.pdf\n")

# Combined multi-panel summary figure (all five variants)
p_combined <- (p_D2_pve | p_D2_assoc) /
              (p_D3_pve | p_D3_assoc) /
              (p_D5_pve | p_D5_assoc) /
              (plot_spacer() | p_D4_assoc) +
  plot_annotation(
    title    = "Delta matrix NMF — PVE and who_delta association",
    subtitle = "D2: Laplace | D3: Normal | D5: Asym EBMF (exp-L/Laplace-F) | D4: GBCD"
  )
ggsave(file.path(fig_dir, "delta_combined_summary.pdf"),
       p_combined, width = 12, height = 16)
cat("Saved: delta_combined_summary.pdf\n")

cat("\n--- 03_NMF_delta.R complete ---\n")
cat("Key outputs:\n")
cat(sprintf("  D1 NMF (K=%d, shifted delta):            delta_D1_nmf.rds\n", K_D1))
cat(sprintf("  D2 EBMF Laplace (auto K=%d, signed L+F): delta_D2_laplace.rds\n", K_D2))
cat(sprintf("  D3 EBMF Normal  (auto K=%d, signed L+F): delta_D3_normal.rds\n", K_D3))
cat(sprintf("  D4 GBCD         (auto K=%d, binary L):   delta_D4_gbcd.rds\n", K_D4))
cat(sprintf("  D5 Asym EBMF    (auto K=%d, nonneg L):   delta_D5_asym_ebmf.rds\n", K_D5))
cat("  Association table: delta_association.csv\n")
cat("\n  D5 is the asymmetric variant: point-exponential prior on L enforces\n")
cat("  non-negative patient participation scores; point-Laplace on F captures\n")
cat("  signed sparse protein changes. Non-neg L enables all flashier plot types.\n")
cat("  D4 GBCD: binary-like patient memberships + signed LFC per protein.\n")
