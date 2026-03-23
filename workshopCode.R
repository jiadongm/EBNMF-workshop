library(ggplot2)
set.seed(42)
n_null <- 490; n_signal <- 10; s <- 1
theta  <- c(rep(0, n_null), rep(5, n_signal))
x      <- rnorm(500, mean = theta, sd = s)

df <- data.frame(
  x    = x,
  type = c(rep("null (theta=0)", n_null), rep("signal (theta=3)", n_signal))
)

ggplot(df, aes(x = x, fill = type)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("grey60", "steelblue")) +
  labs(x = "Observed x_i", y = "Count", fill = "Truth",
       title = "500 genes: 490 null, 10 signal",
       subtitle = "EB learns g from this histogram, then shrinks accordingly") +
  theme_classic()



# Block 2
library(ggplot2)
library(tidyr)

# Shrinkage function for normal prior: E[theta | x] = (sigma^2 / (sigma^2 + s^2)) * x
sigma2 <- 1; s2 <- 1     # equal prior variance and noise variance
x_vals <- seq(-4, 4, length.out = 200)

df_shrink <- data.frame(
  x            = x_vals,
  MLE          = x_vals,
  EB_posterior = (sigma2 / (sigma2 + s2)) * x_vals
)

df_long <- pivot_longer(df_shrink, -x, names_to = "Estimator", values_to = "Estimate")

ggplot(df_long, aes(x = x, y = Estimate, colour = Estimator, linetype = Estimator)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, colour = "grey70") +
  scale_colour_manual(values = c("steelblue", "tomato")) +
  labs(x = expression(paste("Observed ", x[i])),
       y = expression(paste("Estimate of ", theta[i])),
       title = "Shrinkage: EB posterior mean vs MLE (normal prior)",
       subtitle = "EB pulls all estimates toward zero; strong signals survive") +
  theme_classic()


## Block 3
library(ebnm)
fit_pe  <- ebnm_point_exponential(x, s)
fit_pl  <- ebnm_point_laplace(x, s)
fit_pn  <- ebnm_point_normal(x, s)

logLik(fit_pe)   # pick the family with the highest value
logLik(fit_pl)
logLik(fit_pn)


# Block 4
library(ggplot2)
library(patchwork)

x <- seq(-4, 4, length.out = 400)

# Panel 1: Normal prior
p1 <- ggplot(data.frame(x = x, y = dnorm(x, 0, 1)), aes(x, y)) +
  geom_line(colour = "steelblue", linewidth = 1) +
  labs(title = "Normal", x = expression(theta), y = "Density") +
  theme_classic()

# Panel 2: Point-normal (spike + Gaussian slab)
p2 <- ggplot(data.frame(x = x, y = 0.7 * dnorm(x, 0, 1.5)), aes(x, y)) +
  geom_line(colour = "darkorange", linewidth = 1) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1.2),
               colour = "darkorange", linewidth = 2) +
  labs(title = "Point-normal", x = expression(theta), y = "") +
  theme_classic()

# Panel 3: Point-exponential (spike + exponential slab, support >= 0)
x_pos <- seq(0, 4, length.out = 200)
p3 <- ggplot(data.frame(x = x_pos, y = 0.6 * dexp(x_pos, rate = 1)), aes(x, y)) +
  geom_line(colour = "forestgreen", linewidth = 1) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1.5),
               colour = "forestgreen", linewidth = 2) +
  labs(title = "Point-exponential (non-negative)", x = expression(theta), y = "Density") +
  theme_classic()

# Panel 4: Point-Laplace (spike + Laplace slab)
dlaplace <- function(x, b = 1) exp(-abs(x) / b) / (2 * b)
p4 <- ggplot(data.frame(x = x, y = 0.7 * dlaplace(x, b = 1.5)), aes(x, y)) +
  geom_line(colour = "purple", linewidth = 1) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1.2),
               colour = "purple", linewidth = 2) +
  labs(title = "Point-Laplace (heavy tail)", x = expression(theta), y = "") +
  theme_classic()

(p1 + p2) / (p3 + p4) +
  plot_annotation(
    title    = "Key prior shapes in ebnm",
    subtitle = "Vertical segment = point mass at 0; curve = continuous slab"
  )

# Block 5
library(flashier)
fit <- flash(
  data        = X,                        # n x p matrix (rows = cells, columns = genes)
  ebnm_fn     = ebnm_point_exponential,   # prior family; or list c(G_l, G_f) for asymmetric
  greedy_Kmax = 30,                       # maximum factors to add greedily
  backfit     = TRUE,
  nullcheck   = TRUE,
  verbose     = 1L
)




# Block 6
library(flashier); library(ebnm); library(gbcd); library(NMF)

# ── Data preparation ──────────────────────────────────────────────────────────
X_prot   <- readRDS("EBNMF/CovidCase/results/X_prot.rds")      # 240 × 481 (samples × proteins)
meta     <- readRDS("EBNMF/CovidCase/results/pca_scores.rds")  # sample metadata
col_name <- colnames(X_prot)                   # protein names

# Per-protein min-shift → non-negative (required for S1 and S2)
X_shifted <- sweep(X_prot, 2, apply(X_prot, 2, min), "-")

if(F){
  # ── S1: Standard NMF (brunet, K=8) ───────────────────────────────────────────
  set.seed(123)
  tic <- Sys.time()
  fit_S1 <- NMF::nmf(t(X_shifted), rank = 8, method = "brunet", nrun = 10, seed = 123)
  Sys.time() - tic
  
  # ── S2: EBNMF — point-exponential on both L and F ────────────────────────────
  set.seed(42)
  tic <- Sys.time()
  fit_S2 <- flash(
    data        = X_shifted,
    ebnm_fn     = ebnm_point_exponential,   # same prior for L and F
    greedy_Kmax = 15, backfit = TRUE, nullcheck = TRUE, verbose = 1L
  )
  Sys.time() - tic
  
  
  # ── S3: Asymmetric EBMF — point-exponential L, point-Laplace F ───────────────
  set.seed(55)
  tic <- Sys.time()
  fit_S3 <- flash(
    data        = X_prot,
    ebnm_fn     = list(ebnm_point_exponential,   # L: non-negative participation scores
                       ebnm_point_laplace),       # F: signed sparse protein weights
    greedy_Kmax = 15, backfit = TRUE, nullcheck = TRUE, verbose = 1L
  )
  Sys.time() - tic
  
  # ── S4: GBCD — generalised binary prior ──────────────────────────────────────
  set.seed(99)
  tic <- Sys.time()
  fit_S4 <- fit_gbcd(
    Y        = X_prot,
    Kmax     = 15,
    prior    = ebnm::ebnm_generalized_binary,
    maxiter1 = 500, maxiter2 = 200, maxiter3 = 500,
    verbose  = 1
  )
  Sys.time() - tic
  
  saveRDS(list(fit_S1 = fit_S1, fit_S2 = fit_S2, fit_S3 = fit_S3, fit_S4 = fit_S4),
          "EBNMF/CovidCase/results/fit_res_for_workshop.rds")
} else {
  
  fitRes <- readRDS("EBNMF/CovidCase/results/fit_res_for_workshop.rds")
  fit_S1 <- fitRes$fit_S1
  fit_S2 <- fitRes$fit_S2
  fit_S3 <- fitRes$fit_S3
  fit_S4 <- fitRes$fit_S4
}

# NMF package: basis() = W (p×K), coef() = H (K×n)
L_S1 <- t(NMF::coef(fit_S1))   # 240 × 8  (samples × factors)
F_S1 <- NMF::basis(fit_S1)     # 481 × 8  (proteins × factors)
rownames(L_S1) <- rownames(X_shifted)
rownames(F_S1) <- col_name

K_S2   <- fit_S2$n_factors
L_S2   <- fit_S2$L_pm    # 240 × K  (non-negative)
F_S2   <- fit_S2$F_pm    # 481 × K  (non-negative)
pve_S2 <- fit_S2$pve

K_S3   <- fit_S3$n_factors
L_S3   <- fit_S3$L_pm    # 240 × K  (non-negative)
F_S3   <- fit_S3$F_pm    # 481 × K  (signed)
pve_S3 <- fit_S3$pve

K_S4   <- ncol(fit_S4$L)
L_S4   <- fit_S4$L          # 240 × K  (binary-like, [0,1])
F_S4   <- fit_S4$F$lfc      # 481 × K  (signed log-fold change)
F_lfsr <- fit_S4$F$lfsr     # 481 × K  (local false sign rate)
rownames(L_S4) <- rownames(X_shifted)
rownames(F_S4) <- rownames(F_lfsr) <- col_name




# Block 7
library(pheatmap)

# ── Shared helpers ─────────────────────────────────────────────────────────────

# Row order: T1 first, then T2; within each, ascending WHO severity
row_ord <- order(meta$time_point, meta$who_score_num)

# Row annotation: time point + WHO severity
ann_row <- data.frame(
  time_point = meta$time_point,
  WHO        = meta$who_score_num,
  row.names  = rownames(X_shifted)
)

# Helper: per-column max-normalise to [0,1]
# col_norm <- function(M) sweep(M, 2, pmax(apply(M, 2, max), 1e-8), "/")

# Colour scales
heat_pos <- colorRampPalette(c("gray96", "red"))(50)        # non-negative L and F
heat_div <- colorRampPalette(c("#313695", "white", "#A50026"))(100)  # signed F

# ── Loading heatmaps (L matrices) ─────────────────────────────────────────────
# All four variants have non-negative L, so the same scale applies.
# Columns ordered by decreasing PVE for S2/S3; as-fit for S1 and S4.

# S1 — Standard NMF
pheatmap(L_S1[row_ord, ], # col_norm(L_S1)[row_ord, ], 
         annotation_row = ann_row, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat_pos, show_rownames = FALSE, main = "S1 Standard NMF — loadings")

# S2 — EBNMF (point-exponential)
pheatmap(L_S2[row_ord, order(pve_S2, decreasing = TRUE)], #col_norm(L_S2)[row_ord, order(pve_S2, decreasing = TRUE)],
         annotation_row = ann_row, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat_pos, show_rownames = FALSE, main = "S2 EBNMF (point-exp) — loadings")

# S3 — Asymmetric EBMF
pheatmap((L_S3)[row_ord, order(pve_S3, decreasing = TRUE)], #col_norm(L_S3)[row_ord, order(pve_S3, decreasing = TRUE)],
         annotation_row = ann_row, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat_pos, show_rownames = FALSE, main = "S3 Asymmetric EBMF — loadings")

# S4 — GBCD (L already in [0,1]; no normalisation needed)
pheatmap(L_S4[row_ord, ],
         annotation_row = ann_row, cluster_rows = FALSE, cluster_cols = FALSE,
         color = heat_pos, show_rownames = FALSE, main = "S4 GBCD — loadings")

# ── Factor heatmaps (F matrices — protein weights) ────────────────────────────
# S1 and S2: F non-negative → per-column max-normalise, gray96 → red.
# S3: F signed (point-Laplace) → diverging scale.
# S4: F signed LFC; faded where lfsr ≥ 0.05 (handled separately below).

# S1 — non-negative F
pheatmap(t(F_S1), #t(col_norm(F_S1)),
         cluster_rows = FALSE, cluster_cols = TRUE,
         color = heat_pos, show_colnames = FALSE, show_rownames = FALSE,
         main = "S1 Standard NMF — factor weights (proteins)")

# S2 — non-negative F
pheatmap(t((F_S2[, order(pve_S2, decreasing = TRUE)])), #t(col_norm(F_S2[, order(pve_S2, decreasing = TRUE)])),
         cluster_rows = FALSE, cluster_cols = TRUE,
         color = heat_pos, show_colnames = FALSE,
         main = "S2 EBNMF (point-exp) — factor weights (proteins)")

# S3 — signed F; symmetric diverging scale
F_S3_lim <- max(abs(F_S3))
pheatmap(t(F_S3[, order(pve_S3, decreasing = TRUE)]),
         cluster_rows = FALSE, cluster_cols = TRUE,
         color = heat_div, breaks = seq(-F_S3_lim, F_S3_lim, length.out = 101),
         show_colnames = FALSE,
         main = "S3 Asymmetric EBMF — factor weights (signed, proteins)")

# S4 — signed LFC; mask low-confidence entries (lfsr ≥ 0.05 → grey)
F_S4_plot <- F_S4
F_S4_plot[F_lfsr >= 0.05] <- 0   # zero out non-significant entries
F_S4_lim  <- max(abs(F_S4_plot))
pheatmap(t(F_S4_plot),
         cluster_rows = FALSE, cluster_cols = TRUE,
         color = heat_div, breaks = seq(-F_S4_lim, F_S4_lim, length.out = 101),
         show_colnames = FALSE,
         main = "S4 GBCD — factor LFC (lfsr ≥ 0.05 zeroed, proteins)")





library(ggplot2); library(patchwork)

# ── T1 rows and T1 severity ──────────────────────────────────────────────────
T1_idx <- which(meta$time_point == "T1")
who_T1 <- meta$who_score_num[T1_idx]   # 120 patients, all non-NA

# T1 loading for the best factor per variant (indices from the table above)
L_T1 <- list(
  "S1 Standard NMF\nF1 (\u03c1 = +0.750, FDR < 0.001)"       = L_S1[T1_idx, 1],
  "S2 EBNMF (point-exp)\nF1 (\u03c1 = \u22120.819, FDR < 0.001)" = L_S2[T1_idx, 1],
  "S3 Asymmetric EBMF\nF4 (\u03c1 = +0.737, FDR < 0.001)"    = L_S3[T1_idx, 4],
  "S4 GBCD\nF3 (\u03c1 = \u22120.690, FDR < 0.001)"              = L_S4[T1_idx, 3]
)

# ── One scatter panel per variant ───────────────────────────────────────────
plots <- lapply(names(L_T1), function(nm) {
  df <- data.frame(loading = L_T1[[nm]], severity = who_T1)
  ggplot(df, aes(y = severity, x = loading)) +
    geom_jitter(alpha = 0.7, size = 2, width = 0.1, colour = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.7) +
    labs(x = "T1 WHO severity", y = "T1 loading", title = nm) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 9))
})

wrap_plots(plots, nrow = 2)



# Block DIVAS
library(DIVAS)

# ── Build blocks: features × samples (DIVAS convention) ─────────────────────
# Prot_T1 and Prot_T2: 481 proteins × 120 patients, same column order
Prot_T1 <- X_prot[grepl("T1", rownames(X_prot)), ]
Prot_T2 <- X_prot[grepl("T2", rownames(X_prot)), ]
data_list <- list(T1 = t(Prot_T1), T2 = t(Prot_T2))

# ── Fit DIVAS ────────────────────────────────────────────────────────────────
if(F){
  set.seed(123)
  tic <- Sys.time()
  divas_res <- DIVASmain(
    datablock    = data_list,
    nsim         = 400,       # bootstrap simulations for rank estimation
    colCent      = T,        # centre proteins within each block
    rowCent      = F,      
    seed         = 123,
    ReturnDetail = TRUE       # needed for getTopFeatures()
  )
  Sys.time() - tic
  saveRDS(divas_res,  "EBNMF/CovidCase/results/divas_2block_forWorkshop.rds")
  
} else {
  divas_res <- readRDS("EBNMF/CovidCase/results/divas_2block_forWorkshop.rds")
}

# ── Extract results ──────────────────────────────────────────────────────────
scores_mat <- divas_res$Scores   # 120 patients × 24 components
comp_names <- colnames(scores_mat)
cat(sprintf("Components: %d\n", ncol(scores_mat)))
print(table(ifelse(grepl("\\+", comp_names), "Shared", "Individual")))


# ── Score heatmap ────────────────────────────────────────────────────────────
library(pheatmap)
ann_col <- data.frame(
  T1_severity = t1_severity,
  WHO_delta   = who_delta,
  row.names   = rownames(scores_mat)
)
pheatmap(
  t(scores_mat),
  annotation_col = ann_col,
  cluster_rows   = FALSE,
  cluster_cols   = TRUE,
  show_colnames  = FALSE,
  color          = colorRampPalette(c("#313695", "white", "#A50026"))(100),
  main           = "DIVAS scores: proteomics T1 vs T2"
)


# ── Associate with clinical metadata ────────────────────────────────────────
t1_meta <- meta[rownames(Prot_T1),]
t2_meta <- meta[rownames(Prot_T2),]
t1_severity <- t1_meta$who_score_num
t2_severity <- t2_meta$who_score_num

rho_sev <- apply(scores_mat, 2, function(x)
  cor(x, t1_severity, method = "spearman", use = "complete.obs"))
pval_sev <- sapply(seq_len(ncol(scores_mat)), function(k)
  cor.test(scores_mat[, k], t1_severity, method = "spearman")$p.value)
fdr_sev  <- p.adjust(pval_sev, method = "BH")
data.frame(
  rho = rho_sev, pval = pval_sev, pval_adj = fdr_sev
)

rho_sev <- apply(scores_mat, 2, function(x)
  cor(x, who_delta, method = "spearman", use = "complete.obs"))
pval_sev <- sapply(seq_len(ncol(scores_mat)), function(k)
  cor.test(scores_mat[, k], who_delta, method = "spearman")$p.value)
fdr_sev  <- p.adjust(pval_sev, method = "BH")
data.frame(
  rho = rho_sev, pval = pval_sev, pval_adj = fdr_sev
)


data.frame(
  T1T2_1 = scores_mat[,"T1+T2-1"], 
  T1T2_2 = scores_mat[,"T1+T2-2"],
  T1_severity = t1_severity
) |>
  ggplot(aes(T1_severity, T1T2_1)) +
  geom_point() +
  geom_smooth()



# ── Top proteins for a shared component ─────────────────────────────────────
# Shared components: compName = e.g. "T1+T2-1", modName = "T1" or "T2"
feats_T1 <- getTopFeatures(divas_res, compName = "T1+T2-1", modName = "T1",
                           n_top_pos = 10, n_top_neg = 10)
cat("T1+T2-1 top positive proteins (T1):", paste(feats_T1$top_positive, collapse = ", "), "\n")
feats_T2 <- getTopFeatures(divas_res, compName = "T1+T2-1", modName = "T1",
                           n_top_pos = 10, n_top_neg = 10)
cat("T1+T2-1 top positive proteins (T2):", paste(feats_T1$top_positive, collapse = ", "), "\n")

feats_T1 <- getTopFeatures(divas_res, compName = "T1+T2-6", modName = "T1",
                           n_top_pos = 10, n_top_neg = 10)
cat("T1+T2-6 top positive proteins (T1):", paste(feats_T1$top_positive, collapse = ", "), "\n")
feats_T2 <- getTopFeatures(divas_res, compName = "T1+T2-6", modName = "T1",
                           n_top_pos = 10, n_top_neg = 10)
cat("T1+T2-6 top positive proteins (T1):", paste(feats_T1$top_positive, collapse = ", "), "\n")





