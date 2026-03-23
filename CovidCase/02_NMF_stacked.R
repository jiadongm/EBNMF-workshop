# 02_NMF_stacked.R — NMF on stacked 240×481 proteomics matrix
#
# Biological question:
#   What protein co-expression modules exist across all samples (T1 + T2)?
#   Do factor loading scores change from T1 to T2, and does that change
#   correlate with WHO severity change (who_delta = T2 - T1)?
#
# Matrix: X_shifted (240 × 481), per-protein min-shift applied here.
#   X_shifted[,j] = X_prot[,j] - min(X_prot[,j])  → all values ≥ 0
#
# Three NMF variants:
#   S1 — Standard NMF (NMF package, brunet, K=8)
#   S2 — EBNMF (flashier + ebnm_point_exponential, auto K)
#   S3 — GBCD  (gbcd package, generalised binary prior, Kmax=15)
#
# Outputs:
#   results/nmf/stacked_S1_nmf.rds      — NMFfit object
#   results/nmf/stacked_S2_ebnmf.rds    — flashier fit
#   results/nmf/stacked_S3_gbcd.rds     — gbcd fit
#   results/nmf/stacked_association.csv — Spearman corr(ΔL_k, who_delta) per factor/variant
#   results/figs/nmf/stacked_*.pdf      — figures

suppressPackageStartupMessages({
  library(NMF)
  library(flashier)
  library(ebnm)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

# ── Paths ──────────────────────────────────────────────────────────────────────

results_dir <- "results"
nmf_dir     <- file.path(results_dir, "nmf")
fig_dir     <- file.path(results_dir, "figs", "nmf")
dir.create(nmf_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load data ───────────────────────────────────────────────────────────────

X_prot    <- readRDS(file.path(results_dir, "X_prot.rds"))        # 240 × 481
pca_scores <- readRDS(file.path(results_dir, "pca_scores.rds"))   # PC scores + metadata

cat(sprintf("X_prot dimensions: %d samples × %d proteins\n", nrow(X_prot), ncol(X_prot)))
cat(sprintf("X_prot range: [%.3f, %.3f]\n", min(X_prot), max(X_prot)))

# ── 2. Build X_shifted (per-protein min-shift) ─────────────────────────────────

prot_min  <- apply(X_prot, 2, min)     # per-protein minimum
X_shifted <- sweep(X_prot, 2, prot_min, "-")   # subtract per-protein min → all ≥ 0

cat(sprintf("X_shifted range: [%.3f, %.3f]\n", min(X_shifted), max(X_shifted)))
stopifnot(min(X_shifted) >= 0)

# ── 3. Metadata alignment ──────────────────────────────────────────────────────

# pca_scores has one row per sample (240 rows), columns include:
#   patient_id, time_point, who_delta, who_score_num, who_ambiguous
# Align to X_prot row order
meta <- pca_scores[match(rownames(X_prot), pca_scores$sample), ]
stopifnot(all(meta$sample == rownames(X_prot)))

cat(sprintf("who_delta range: [%d, %d]; NA count (ambiguous): %d\n",
            min(meta$who_delta, na.rm = TRUE),
            max(meta$who_delta, na.rm = TRUE),
            sum(is.na(meta$who_delta))))

# Patient-level metadata (T1 rows, one per patient for association analysis)
meta_T1 <- meta[meta$time_point == "T1", ]
cat(sprintf("T1 patients with non-NA who_delta: %d\n",
            sum(!is.na(meta_T1$who_delta))))

# ── Helper: association analysis ───────────────────────────────────────────────
# For each factor k, compute Spearman corr(L_T2[,k] - L_T1[,k], who_delta)
# L must be ordered to match X_prot rows (T1 and T2 interleaved or stacked).

compute_associations <- function(L_mat, meta_df, variant_label) {
  # meta_df: 240 rows aligned to L_mat rows
  # L_mat: 240 × K non-negative factor scores
  K <- ncol(L_mat)
  cat(sprintf("\n--- Association analysis: %s (K=%d) ---\n", variant_label, K))

  # Split into T1 and T2 (match by patient_id)
  # We need paired rows: same patient, T1 and T2
  T1_idx <- which(meta_df$time_point == "T1")
  T2_idx <- which(meta_df$time_point == "T2")

  # Match T2 to T1 by patient_id
  pid_T1 <- meta_df$patient_id[T1_idx]
  pid_T2 <- meta_df$patient_id[T2_idx]
  shared_patients <- intersect(pid_T1, pid_T2)

  T1_order <- T1_idx[match(shared_patients, pid_T1)]
  T2_order <- T2_idx[match(shared_patients, pid_T2)]

  L_T1   <- L_mat[T1_order, , drop = FALSE]
  L_T2   <- L_mat[T2_order, , drop = FALSE]
  delta_L <- L_T2 - L_T1

  who_delta_vec <- meta_df$who_delta[T1_order]
  non_na        <- !is.na(who_delta_vec)

  results <- lapply(seq_len(K), function(k) {
    if (sum(non_na) < 5) return(data.frame(variant = variant_label, factor = k,
                                            spearman_rho = NA, p_value = NA))
    ct <- cor.test(delta_L[non_na, k], who_delta_vec[non_na], method = "spearman",
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

  # Print top hits
  top <- assoc_df[!is.na(assoc_df$spearman_rho) & abs(assoc_df$spearman_rho) > 0.1, ]
  if (nrow(top) > 0) {
    cat("  Factors with |rho| > 0.1:\n")
    print(top, row.names = FALSE)
  } else {
    cat("  No factors with |rho| > 0.1 found.\n")
  }
  return(assoc_df)
}

# ── Helper: factor-level trajectory plot ──────────────────────────────────────

plot_factor_trajectories <- function(L_mat, meta_df, variant_label, assoc_df = NULL,
                                     max_factors = 9) {
  K <- min(ncol(L_mat), max_factors)
  plots <- vector("list", K)
  for (k in seq_len(K)) {
    df_k <- data.frame(
      patient_id  = meta_df$patient_id,
      time_point  = meta_df$time_point,
      score       = L_mat[, k],
      who_delta   = meta_df$who_delta
    )
    # Pivot to wide: T1 vs T2
    df_wide <- df_k |>
      dplyr::filter(!is.na(who_delta)) |>
      tidyr::pivot_wider(names_from = time_point, values_from = score,
                         id_cols = c(patient_id, who_delta)) |>
      dplyr::filter(!is.na(T1) & !is.na(T2))

    if (nrow(df_wide) < 3) next

    rho_label <- ""
    if (!is.null(assoc_df)) {
      rho_row <- assoc_df[assoc_df$factor == k & assoc_df$variant == variant_label, ]
      if (nrow(rho_row) > 0 && !is.na(rho_row$spearman_rho)) {
        rho_label <- sprintf("rho=%.2f, FDR=%.2f",
                             rho_row$spearman_rho, rho_row$bh_fdr)
      }
    }

    plots[[k]] <- ggplot(df_wide, aes(x = T1, y = T2, colour = who_delta)) +
      geom_point(size = 1.8, alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50",
                  linewidth = 0.4) +
      scale_colour_gradient2(low = "steelblue", mid = "gold", high = "firebrick",
                             midpoint = 0, name = "dWHO", na.value = "grey80") +
      labs(title = sprintf("Factor %d  %s", k, rho_label),
           x = "L score (T1)", y = "L score (T2)") +
      theme_bw(base_size = 9)
  }
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(invisible(NULL))
  wrap_plots(plots, ncol = 3) +
    plot_annotation(title = sprintf("%s — factor loading trajectories (T1 → T2)",
                                    variant_label))
}

# ── Helper: PVE bar chart ──────────────────────────────────────────────────────

plot_pve <- function(pve_vec, title_str) {
  df <- data.frame(factor = seq_along(pve_vec), pve = pve_vec * 100)
  ggplot(df, aes(x = factor(factor), y = pve)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", pve)), vjust = -0.3, size = 2.8) +
    labs(title = title_str, x = "Factor", y = "% variance explained") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(size = 8))
}

# ── Helper: factor-protein bar plots ──────────────────────────────────────────

plot_factor_proteins <- function(F_mat, top_n = 20, variant_label = "", max_factors = 6) {
  K       <- min(ncol(F_mat), max_factors)
  prot_names <- rownames(F_mat)
  plots   <- vector("list", K)
  for (k in seq_len(K)) {
    f_k <- F_mat[, k]
    ord <- order(f_k, decreasing = TRUE)
    idx <- c(ord[1:min(top_n, length(ord))])  # top positives only (NMF >= 0)
    df  <- data.frame(
      protein = prot_names[idx],
      weight  = f_k[idx]
    ) |>
      dplyr::arrange(dplyr::desc(weight)) |>
      dplyr::mutate(protein = factor(protein, levels = rev(protein)))

    plots[[k]] <- ggplot(df, aes(x = weight, y = protein)) +
      geom_col(fill = "steelblue", width = 0.7) +
      labs(title = sprintf("Factor %d", k), x = "Weight", y = NULL) +
      theme_bw(base_size = 8) +
      theme(axis.text.y = element_text(size = 6))
  }
  wrap_plots(plots, ncol = 3) +
    plot_annotation(title = sprintf("%s — top proteins per factor", variant_label))
}

# ══════════════════════════════════════════════════════════════════════════════
# VARIANT S1 — Standard NMF (brunet, K=8)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== S1: Standard NMF (K=8) ==========\n")

K_S1 <- 8

if (!file.exists(file.path(nmf_dir, "stacked_S1_nmf.rds"))) {
  set.seed(123)
  t_start <- proc.time()
  fit_S1 <- NMF::nmf(
    x      = t(X_shifted),   # NMF package expects features × samples (p × n)
    rank   = K_S1,
    method = "brunet",
    nrun   = 10,
    seed   = 123
  )
  t_S1 <- proc.time() - t_start
  saveRDS(fit_S1, file.path(nmf_dir, "stacked_S1_nmf.rds"))
  cat(sprintf("S1 NMF fitted and saved.  [%.1f s]\n", t_S1["elapsed"]))
} else {
  fit_S1 <- readRDS(file.path(nmf_dir, "stacked_S1_nmf.rds"))
  cat("S1 NMF loaded from cache.\n")
}

# Extract matrices — NMF package: basis() = W (p×K), coef() = H (K×n)
# Our convention: L (n×K loadings), F (p×K factors)
W_S1 <- NMF::basis(fit_S1)   # 481 × 8 (proteins × factors)
H_S1 <- NMF::coef(fit_S1)    # 8 × 240 (factors × samples)
L_S1 <- t(H_S1)              # 240 × 8 (samples × factors)
F_S1 <- W_S1                  # 481 × 8 (proteins × factors)

cat(sprintf("S1: W %d×%d, H %d×%d → L %d×%d, F %d×%d\n",
            nrow(W_S1), ncol(W_S1), nrow(H_S1), ncol(H_S1),
            nrow(L_S1), ncol(L_S1), nrow(F_S1), ncol(F_S1)))

assoc_S1 <- compute_associations(L_S1, meta, "S1_NMF")

# Factor-protein plots
p_S1_factors <- plot_factor_proteins(F_S1, variant_label = "S1 Standard NMF")
ggsave(file.path(fig_dir, "stacked_S1_factor_proteins.pdf"),
       p_S1_factors, width = 14, height = 10)
cat("Saved: stacked_S1_factor_proteins.pdf\n")

# Trajectory plots
p_S1_traj <- plot_factor_trajectories(L_S1, meta, "S1_NMF", assoc_S1)
if (!is.null(p_S1_traj)) {
  ggsave(file.path(fig_dir, "stacked_S1_trajectories.pdf"),
         p_S1_traj, width = 12, height = 10)
  cat("Saved: stacked_S1_trajectories.pdf\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# VARIANT S2 — EBNMF (flashier, point-exponential, auto K)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== S2: EBNMF (point-exponential, auto K) ==========\n")

if (!file.exists(file.path(nmf_dir, "stacked_S2_ebnmf.rds"))) {
  set.seed(42)
  t_start <- proc.time()
  fit_S2 <- flash(
    data        = X_shifted,           # 240 × 481 (n × p)
    ebnm_fn     = ebnm_point_exponential,
    greedy_Kmax = 30,
    backfit     = TRUE,
    nullcheck   = TRUE,
    verbose     = 1L
  )
  t_S2 <- proc.time() - t_start
  saveRDS(fit_S2, file.path(nmf_dir, "stacked_S2_ebnmf.rds"))
  cat(sprintf("S2 EBNMF fitted and saved.  [%.1f s]\n", t_S2["elapsed"]))
} else {
  fit_S2 <- readRDS(file.path(nmf_dir, "stacked_S2_ebnmf.rds"))
  cat("S2 EBNMF loaded from cache.\n")
}

K_S2   <- fit_S2$n_factors
L_S2   <- fit_S2$L_pm    # 240 × K
F_S2   <- fit_S2$F_pm    # 481 × K
pve_S2 <- fit_S2$pve

cat(sprintf("S2 auto K=%d; PVE sum=%.3f\n", K_S2, sum(pve_S2)))
cat(sprintf("S2 PVE per factor: %s\n",
            paste(round(pve_S2 * 100, 1), collapse = ", ")))
stopifnot(sum(pve_S2) < 1 || K_S2 == 0)  # automatic stopping check

assoc_S2 <- compute_associations(L_S2, meta, "S2_EBNMF")

# PVE plot
p_S2_pve <- plot_pve(pve_S2, "S2 EBNMF — PVE per factor (point-exponential)")
ggsave(file.path(fig_dir, "stacked_S2_pve.pdf"), p_S2_pve, width = 7, height = 4)
cat("Saved: stacked_S2_pve.pdf\n")

# Factor-protein plots
p_S2_factors <- plot_factor_proteins(F_S2, variant_label = "S2 EBNMF", max_factors = K_S2)
ggsave(file.path(fig_dir, "stacked_S2_factor_proteins.pdf"),
       p_S2_factors, width = 14, height = 10)
cat("Saved: stacked_S2_factor_proteins.pdf\n")

# Trajectory plots
p_S2_traj <- plot_factor_trajectories(L_S2, meta, "S2_EBNMF", assoc_S2)
if (!is.null(p_S2_traj)) {
  ggsave(file.path(fig_dir, "stacked_S2_trajectories.pdf"),
         p_S2_traj, width = 12, height = 10)
  cat("Saved: stacked_S2_trajectories.pdf\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# VARIANT S3 — GBCD (generalised binary prior)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== S3: GBCD (Kmax=15) ==========\n")

if (!requireNamespace("gbcd", quietly = TRUE)) {
  cat("gbcd not installed. Installing from GitHub...\n")
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("stephenslab/gbcd")
}
library(gbcd)

if (!file.exists(file.path(nmf_dir, "stacked_S3_gbcd.rds"))) {
  set.seed(99)
  # fit_gbcd expects Y: samples × features (n × p), log-normalised
  # X_shifted is 240 × 481 after min-shift — suitable input
  t_start <- proc.time()
  fit_S3 <- fit_gbcd(
    Y        = X_shifted,
    Kmax     = 15,
    prior    = ebnm::ebnm_generalized_binary,
    maxiter1 = 500,
    maxiter2 = 200,
    maxiter3 = 500,
    verbose  = 1
  )
  t_S3 <- proc.time() - t_start
  saveRDS(fit_S3, file.path(nmf_dir, "stacked_S3_gbcd.rds"))
  cat(sprintf("S3 GBCD fitted and saved.  [%.1f s]\n", t_S3["elapsed"]))
} else {
  fit_S3 <- readRDS(file.path(nmf_dir, "stacked_S3_gbcd.rds"))
  cat("S3 GBCD loaded from cache.\n")
}

L_S3   <- fit_S3$L          # 240 × K — patient memberships (non-negative, binary-like)
F_lfc  <- fit_S3$F$lfc      # 481 × K — signed log-fold change per protein per program
F_lfsr <- fit_S3$F$lfsr     # 481 × K — local false sign rates
K_S3   <- ncol(L_S3)

cat(sprintf("S3 GBCD K=%d (max 2*Kmax-1 = %d)\n", K_S3, 2*15-1))

assoc_S3 <- compute_associations(L_S3, meta, "S3_GBCD")

# Factor-protein plot using LFC (signed)
plot_gbcd_factor_proteins <- function(F_lfc_mat, F_lfsr_mat, top_n = 20,
                                       max_factors = 6) {
  K      <- min(ncol(F_lfc_mat), max_factors)
  pnames <- rownames(F_lfc_mat)
  plots  <- vector("list", K)
  for (k in seq_len(K)) {
    lfc_k  <- F_lfc_mat[, k]
    lfsr_k <- F_lfsr_mat[, k]
    # Top by absolute LFC
    ord <- order(abs(lfc_k), decreasing = TRUE)[1:min(top_n, length(lfc_k))]
    df  <- data.frame(
      protein = pnames[ord],
      lfc     = lfc_k[ord],
      lfsr    = lfsr_k[ord],
      sig     = lfsr_k[ord] < 0.05
    ) |>
      dplyr::arrange(dplyr::desc(lfc)) |>
      dplyr::mutate(protein = factor(protein, levels = rev(protein)))

    plots[[k]] <- ggplot(df, aes(x = lfc, y = protein,
                                  fill = ifelse(lfc > 0, "up", "down"),
                                  alpha = ifelse(sig, 1, 0.4))) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = c(up = "steelblue", down = "firebrick"),
                        guide = "none") +
      scale_alpha_identity() +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      labs(title = sprintf("Factor %d", k),
           subtitle = sprintf("lfsr<0.05: %d proteins", sum(df$sig)),
           x = "LFC (GBCD)", y = NULL) +
      theme_bw(base_size = 8) +
      theme(axis.text.y = element_text(size = 6))
  }
  wrap_plots(plots, ncol = 3) +
    plot_annotation(title = "S3 GBCD — top proteins per factor (LFC, shaded = lfsr≥0.05)")
}

p_S3_factors <- plot_gbcd_factor_proteins(F_lfc, F_lfsr, max_factors = K_S3)
ggsave(file.path(fig_dir, "stacked_S3_factor_proteins.pdf"),
       p_S3_factors, width = 14, height = 10)
cat("Saved: stacked_S3_factor_proteins.pdf\n")

# Trajectory plots
p_S3_traj <- plot_factor_trajectories(L_S3, meta, "S3_GBCD", assoc_S3)
if (!is.null(p_S3_traj)) {
  ggsave(file.path(fig_dir, "stacked_S3_trajectories.pdf"),
         p_S3_traj, width = 12, height = 10)
  cat("Saved: stacked_S3_trajectories.pdf\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# COMBINED: comparison overview + save association table
# ══════════════════════════════════════════════════════════════════════════════

cat("\n========== Cross-variant summary ==========\n")

# Number of factors summary
summary_df <- data.frame(
  variant       = c("S1_NMF", "S2_EBNMF", "S3_GBCD"),
  K             = c(K_S1, K_S2, K_S3),
  K_selection   = c("fixed (8)", "auto (ELBO)", "auto (GBCD)"),
  top_rho       = c(
    max(abs(assoc_S1$spearman_rho), na.rm = TRUE),
    max(abs(assoc_S2$spearman_rho), na.rm = TRUE),
    max(abs(assoc_S3$spearman_rho), na.rm = TRUE)
  )
)
cat("\nVariant summary:\n")
print(summary_df, row.names = FALSE)

# Combined association table
all_assoc <- rbind(assoc_S1, assoc_S2, assoc_S3)
write.csv(all_assoc,
          file.path(nmf_dir, "stacked_association.csv"),
          row.names = FALSE)
cat("\nSaved: stacked_association.csv\n")

# Side-by-side PVE comparison (S1 uses RSS reduction proxy; S2 has direct PVE)
p_summary <- ggplot(all_assoc[!is.na(all_assoc$spearman_rho), ],
                    aes(x = factor(factor), y = spearman_rho, fill = variant)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c(S1_NMF = "#4E79A7", S2_EBNMF = "#F28E2B",
                                S3_GBCD = "#59A14F")) +
  facet_wrap(~variant, scales = "free_x", ncol = 1) +
  labs(title = "Stacked matrix — Spearman rho(ΔL_k, who_delta) per factor",
       subtitle = "Dashed lines: |rho| = 0.2 threshold",
       x = "Factor", y = "Spearman rho") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "stacked_association_summary.pdf"),
       p_summary, width = 9, height = 10)
cat("Saved: stacked_association_summary.pdf\n")

cat("\n--- 02_NMF_stacked.R complete ---\n")
cat("Key outputs:\n")
cat(sprintf("  S1 NMF (K=%d): stacked_S1_nmf.rds\n", K_S1))
cat(sprintf("  S2 EBNMF (auto K=%d): stacked_S2_ebnmf.rds\n", K_S2))
cat(sprintf("  S3 GBCD (K=%d): stacked_S3_gbcd.rds\n", K_S3))
cat("  Association table: stacked_association.csv\n")
