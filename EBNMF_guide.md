# EBNMF Practical Guide: flashier and GBCD

> Reference for future Claude Code sessions.
> All patterns confirmed on the COVID-19 proteomics CovidCase project (2026-02/03).
> flashier v1.0.58, gbcd from stephenslab/gbcd.

---

## 1. Packages and installation

```r
# CRAN
install.packages("flashier")   # v1.0.58 — confirmed working version
install.packages("ebnm")

# GBCD from GitHub (also updates flashier)
remotes::install_github("stephenslab/gbcd")
# NOTE: gbcd install updates flashier. If flashier was already loaded in the
# current R session before gbcd install, restart R and reload. Run the script
# twice: first run installs; second run runs cleanly from cache.

library(flashier)
library(ebnm)
library(gbcd)
```

---

## 2. When to use which prior (decision table)

| Prior | Function | Sign constraint | Sparsity | Use case |
|-------|----------|-----------------|----------|----------|
| `ebnm_point_exponential` | L ≥ 0, F ≥ 0 | Non-negative only | Moderate | NMF on counts/expression (both matrices non-negative) |
| `ebnm_point_laplace` | Signed | Sparse (double-exponential) | Moderate–high | Sparse signed decomposition (delta matrices, T2−T1, LFC-like data) |
| `ebnm_point_normal` | Signed | Modest | Low | Dense signed decomposition (comparison baseline, PCA-like) |
| `ebnm_generalized_binary` | Binary-like ∈ [0,1] | Non-negative | Structured | Patient subgroup memberships (GBCD, binary in/out programs) |
| `list(ebnm_point_exponential, ebnm_point_laplace)` | L ≥ 0, F signed | Asymmetric | Moderate | Non-neg participation scores + signed gene/protein programs |

**Practical rule:** fit two candidate priors, compare PVE sum and association signal, pick
the one that finds interpretable structure. For bioinformatics NMF (counts, expression):
use `ebnm_point_exponential`. For change matrices (T2−T1, signed data): use
`ebnm_point_laplace` or asymmetric `list(ebnm_point_exponential, ebnm_point_laplace)`.

---

## 3. Data preparation

### Input orientation

```r
# flash() expects:  n_samples × p_features  (rows = samples/cells/patients)
# This is OPPOSITE to the NMF package (which expects features × samples)

# Example: 240 patients × 481 proteins
dim(X)   # should be [240, 481]
```

### Non-negative input for point-exponential (NMF)

```r
# Per-column (per-feature) min-shift: ensures all values >= 0
prot_min  <- apply(X, 2, min)        # per-feature minimum
X_shifted <- sweep(X, 2, prot_min, "-")   # subtract per-feature min
stopifnot(min(X_shifted) >= 0)

# Range check
cat(sprintf("X_shifted range: [%.3f, %.3f]\n", min(X_shifted), max(X_shifted)))
```

### Signed input for point-Laplace / point-normal (no shift needed)

```r
# delta_mat = T2 - T1 per patient (signed, negative = improvement)
delta_mat <- X_T2 - X_T1     # 120 × 481, range approximately [-6, 7]
# No shift required — point-Laplace handles signed data
```

### Asymmetric case (non-neg L, signed F)

```r
# Same signed input as point-Laplace
# The point-exponential prior on L will enforce non-negativity on patient scores
# The point-Laplace prior on F will capture signed protein changes
# Min-shift NOT needed — L non-negativity comes from the prior, not the data
```

### Set rownames before fitting (preserved in fit object)

```r
rownames(X_shifted) <- sample_ids   # preserved in fit$L_pm rownames
colnames(X_shifted) <- feature_ids  # preserved in fit$F_pm rownames
```

---

## 4. Fitting with flash()

### Full function signature

```r
fit <- flash(
  data        = X,                      # n × p matrix (required)
  ebnm_fn     = ebnm_point_exponential, # prior (see §7 for all options)
  greedy_Kmax = 30,                     # max factors (auto-stops earlier)
  backfit     = TRUE,                   # refine factors after greedy pass
  nullcheck   = TRUE,                   # remove factors that don't help
  verbose     = 1L                      # 0=silent, 1=progress, 2=detailed
)
```

### ebnm_fn: single function vs list

```r
# Single function: same prior applied to BOTH L and F
fit <- flash(data = X, ebnm_fn = ebnm_point_exponential, ...)
# → L and F are both non-negative

# List of two: first applies to L, second to F
fit <- flash(
  data    = X,
  ebnm_fn = list(ebnm_point_exponential,   # L: non-negative scores
                 ebnm_point_laplace),       # F: signed sparse features
  ...
)
# → L >= 0, F signed — the asymmetric (D5) variant
```

### Auto K via ELBO

flashier selects K automatically: the greedy pass adds factors until adding another
does not improve ELBO. `greedy_Kmax` is an upper bound, not a target. Nullcheck may
remove additional factors. Final K is almost always < `greedy_Kmax`.

### Runtimes (COVID-19 proteomics, 240×481)

| Variant | greedy_Kmax | K selected | Runtime |
|---------|-------------|------------|---------|
| point-exponential (S2) | 30 | 14 | 40 s |
| point-Laplace (D2) | 20 | 20 | 29 s |
| point-normal (D3) | 20 | 20 | 13 s |
| asymmetric exp-Laplace (D5) | 20 | 20 | 57 s |

---

## 5. Accessing results

```r
K   <- fit$n_factors   # auto-selected K (integer)
L   <- fit$L_pm        # n × K  posterior mean loadings
F   <- fit$F_pm        # p × K  posterior mean factors
pve <- fit$pve         # K-vector, proportion of variance explained

# DO NOT USE (do not exist in flashier v1.0.58):
# flash_get_ldf()    ← does not exist
# nfactors(fit)      ← does not exist

# PVE summary
cat(sprintf("K=%d; PVE sum=%.3f\n", K, sum(pve)))
cat(sprintf("PVE per factor: %s\n", paste(round(pve * 100, 1), collapse = ", ")))

# Rownames: preserved from input matrix
# fit$L_pm has rownames from rownames(data)
# fit$F_pm has rownames from colnames(data)

# Top features in factor k (unsigned)
top_genes_k <- order(F[, k], decreasing = TRUE)[1:20]
rownames(F)[top_genes_k]

# Top features in factor k (by absolute weight, for signed F)
top_genes_k <- order(abs(F[, k]), decreasing = TRUE)[1:20]
```

---

## 6. Visualisation with plot.flash()

### Which plot types work with signed vs non-negative L

| `plot_type` | `pm_which` | Non-neg L | Signed L | Notes |
|-------------|------------|-----------|----------|-------|
| `"scree"` | — | ✓ | ✓ | Always works |
| `"histogram"` | `"loadings"` | ✓ | ✓ | Works with either |
| `"bar"` | `"loadings"` | ✓ | ✓ | Works with either |
| `"scatter"` | `"factors"` | ✓ | ✓ | Works with either |
| `"heatmap"` | `"factors"` | ✓ | ✓ | Works WITHOUT pm_groups |
| `"structure"` | `"loadings"` | ✓ | **✗** | Calls `verify.nonnegative.matrix()` — FAILS with signed L |
| `"heatmap"` | `"loadings"` | ✓ | **✗** | Same check — FAILS with signed L |
| `"heatmap"` | `"factors"` + `pm_groups` | ✓ | **✗** | pm_groups triggers select_loadings — FAILS with signed L |

**Key rule:** Any plot type that uses `pm_groups` with loadings (L) requires non-negative L.
If L is signed (point-Laplace or point-normal prior), use the asymmetric list prior
`list(ebnm_point_exponential, ebnm_point_laplace)` on the same data to get non-neg L.

### pm_groups vs pm_colors semantics

```r
# STRUCTURE PLOT:
# pm_groups = sample labels (patient groupings, controls gaps + x-axis text)
# pm_colors = K-length vector of factor fill colours (NOT sample group colours!)

plot(fit, plot_type = "structure",
     pm_groups = time_point_vector,   # sample grouping (n-vector of labels)
     pm_colors = rainbow(K))          # K colours for K factors (not groups!)

# HISTOGRAM / BAR PLOT:
# pm_groups = sample labels (used to colour samples)
# pm_colors = length(unique(pm_groups)) colours, in SORTED level order
#             NOT the order you might expect — use sorted unique values

sev_uniq   <- sort(unique(as.character(sev_groups)))  # sorted unique levels
sev_colors <- colorRampPalette(c("#F7FBFF", "#08306B"))(length(sev_uniq))
# sev_colors[1] = colour for lowest severity, etc.

plot(fit, plot_type = "histogram", pm_which = "loadings",
     pm_groups = sev_groups,
     pm_colors = sev_colors)

# HEATMAP of factors (no pm_groups — safe with signed or non-neg L):
plot(fit, plot_type = "heatmap", pm_which = "factors",
     kset = c(5, 1, 8))  # subset of factors to show
```

---

## 7. Prior choice guide (with biological examples)

### point-exponential (NMF)

```r
# Both L and F non-negative.
# Best for: count/expression data shifted to [0, ∞).
# Example: stacked proteomics T1+T2 (after per-protein min-shift)
fit <- flash(data = X_shifted, ebnm_fn = ebnm_point_exponential,
             greedy_Kmax = 30, backfit = TRUE, nullcheck = TRUE)
```

### point-Laplace (sparse signed)

```r
# Both L and F signed, sparse (double-exponential = Laplace prior).
# Best for: signed change matrices (T2−T1), LFC-like data.
# Example: delta NPX matrix (T2 - T1 per patient, range ~ [-6, 7])
fit <- flash(data = delta_mat, ebnm_fn = ebnm_point_laplace,
             greedy_Kmax = 20, backfit = TRUE, nullcheck = TRUE)
# WARNING: signed L blocks structure/heatmap plots (see §6)
```

### point-normal (dense signed)

```r
# Both L and F signed, dense (Gaussian shrinkage).
# Best for: comparison baseline; PCA-like decomposition.
# Finds more factors than point-Laplace on same data (less aggressive shrinkage).
fit <- flash(data = delta_mat, ebnm_fn = ebnm_point_normal,
             greedy_Kmax = 20, backfit = TRUE, nullcheck = TRUE)
```

### Asymmetric (non-neg L, signed F) — recommended for change matrices

```r
# L: non-negative (patient participation >= 0)
# F: signed sparse (protein increases positive, decreases negative)
# Advantages: all flashier plots work (non-neg L); F captures direction explicitly.
# Trade-off: slightly weaker association signal than point-Laplace in some datasets.
fit_D5 <- flash(
  data        = delta_mat,                         # signed, no shift needed
  ebnm_fn     = list(ebnm_point_exponential,       # L prior: non-negative
                     ebnm_point_laplace),           # F prior: signed sparse
  greedy_Kmax = 20,
  backfit     = TRUE,
  nullcheck   = TRUE,
  verbose     = 1L
)
# D5 result (COVID-19, 120×481): K=20, PVE=0.840, L >= 0, F signed
```

---

## 8. GBCD

GBCD (Generalized Binary Covariance Decomposition) finds binary-like patient subgroup
memberships (L ∈ [0,1]) with signed LFC per feature (F$lfc) and uncertainty (F$lfsr).

### fit_gbcd() signature

```r
library(gbcd)

fit <- fit_gbcd(
  Y        = X,                              # n × p (samples × features), log-normalised
  Kmax     = 15,                             # upper bound; auto K ≤ 2*Kmax - 1
  prior    = ebnm::ebnm_generalized_binary,  # standard GBCD prior
  maxiter1 = 500,                            # greedy flash iterations
  maxiter2 = 200,                            # signature estimation iterations
  maxiter3 = 500,                            # backfit iterations
  verbose  = 1
)
```

### Accessing results

```r
K      <- ncol(fit$L)   # auto-selected K (up to 2*Kmax - 1)
L      <- fit$L          # n × K, binary-like memberships in [0,1]
F_lfc  <- fit$F$lfc     # p × K, signed log-fold change per feature per program
F_lfsr <- fit$F$lfsr    # p × K, local false sign rate per feature per program

# DO NOT attempt fit$n_factors or fit$pve — GBCD does not have these fields
# GBCD has its own result structure (not a flashier fit object)

# Set rownames if missing
if (is.null(rownames(L)))     rownames(L)     <- rownames(Y)
if (is.null(rownames(F_lfc))) rownames(F_lfc) <- colnames(Y)

# Top features for program k by absolute LFC
top_k  <- order(abs(F_lfc[, k]), decreasing = TRUE)[1:20]
top_up   <- rownames(F_lfc)[F_lfc[, k] > 0 & order(F_lfc[, k], decreasing = TRUE)]  # positive
top_down <- rownames(F_lfc)[F_lfc[, k] < 0 & order(F_lfc[, k], decreasing = FALSE)] # negative

# Significant features (lfsr < 0.05)
sig_k <- rownames(F_lfc)[F_lfsr[, k] < 0.05]
```

### GBCD vs flashier comparison

| Property | flashier (point-exponential) | GBCD |
|----------|------------------------------|------|
| L structure | Continuous, non-negative | Binary-like ∈ [0,1] |
| F structure | Non-negative weights | Signed LFC |
| F uncertainty | Not available directly | `F$lfsr` per feature |
| K selection | Auto (ELBO) | Auto (bootstrap-like) |
| Runtime | Fast (~minutes) | Slow (~5–10 min for typical size) |
| Best for | Gradients, continuous programs | Subgroup memberships, IN/OUT |

### Runtime note

GBCD is slow due to three-stage optimisation (greedy flash → signature estimation → backfit).
**COVID-19 proteomics (240×481, Kmax=15):** 471 s (7.9 min).
**COVID-19 delta matrix (120×481, Kmax=10):** ~4 min.
Plan for at least 5–10 minutes; cache results to disk.

```r
if (!file.exists("results/gbcd_fit.rds")) {
  set.seed(99)
  fit <- fit_gbcd(Y = X, Kmax = 15, ...)
  saveRDS(fit, "results/gbcd_fit.rds")
} else {
  fit <- readRDS("results/gbcd_fit.rds")
}
```

---

## 9. Known gotchas (consolidated)

### flashier

1. **`flash_get_ldf()` and `nfactors()` do not exist** in flashier v1.0.58. Use `fit$L_pm`,
   `fit$F_pm`, `fit$n_factors`, `fit$pve` directly.

2. **Structure plot requires non-negative L.** Any plot calling `verify.nonnegative.matrix()`
   will fail if L has negative values. These include: `plot_type="structure"`,
   `plot_type="heatmap"` with `pm_which="loadings"`, and `plot_type="heatmap"` with
   `pm_which="factors"` when `pm_groups` is also passed (it triggers select_loadings).
   **Fix:** use the asymmetric list prior `list(ebnm_point_exponential, ...)` to enforce
   non-negative L.

3. **`pm_colors` semantics differ by plot type.** In structure plots, `pm_colors` = K
   factor fill colours. In histogram/bar, `pm_colors` = colours for each unique group
   level in **sorted order** — not the order samples appear in the data.

4. **`pm_groups` in factor heatmap triggers the L check.** Passing `pm_groups` to
   `plot_type="heatmap", pm_which="factors"` fails with signed L. Omit `pm_groups` for
   factor heatmaps when L is signed.

5. **gbcd install updates flashier.** Installing gbcd via remotes::install_github()
   upgrades flashier from any earlier version. If flashier was already loaded in the
   current R session, restart R before running the next script.

6. **Input orientation.** `flash()` expects n × p (samples in rows). This is OPPOSITE
   to the `NMF` package which expects p × n (features in rows). Always check `dim(X)`.

7. **Rownames are preserved** from the input matrix. Set `rownames(X)` to sample IDs
   and `colnames(X)` to feature IDs before calling `flash()`, and they will appear
   in `fit$L_pm` and `fit$F_pm` respectively.

### GBCD

8. **GBCD L range with signed input.** When the input Y is signed (e.g., a delta matrix),
   GBCD L can have slight negatives (e.g., range −0.54 to 1.0) even though the GB prior
   targets [0,1]. This is expected; values near 0 and 1 are still interpretable as
   "out" and "in" memberships.

9. **GBCD result structure ≠ flashier.** `fit_gbcd()` returns a plain list, not a
   flashier fit object. Do not call `plot(fit_gbcd_result, ...)` — it will not work.
   Use `fit$L`, `fit$F$lfc`, `fit$F$lfsr` directly; build custom visualisations.

---

## 10. Quick-start templates

### Standard EBNMF (NMF with EB priors)

```r
library(flashier); library(ebnm)

X_shifted <- sweep(X, 2, apply(X, 2, min), "-")  # non-negative
rownames(X_shifted) <- sample_ids
colnames(X_shifted) <- feature_ids

set.seed(42)
fit <- flash(data = X_shifted, ebnm_fn = ebnm_point_exponential,
             greedy_Kmax = 30, backfit = TRUE, nullcheck = TRUE, verbose = 1L)
saveRDS(fit, "ebnmf_fit.rds")

K   <- fit$n_factors
L   <- fit$L_pm     # samples × K
F   <- fit$F_pm     # features × K
pve <- fit$pve
```

### Asymmetric EBMF (best for signed change matrices)

```r
set.seed(55)
fit <- flash(
  data        = delta_mat,
  ebnm_fn     = list(ebnm_point_exponential, ebnm_point_laplace),
  greedy_Kmax = 20, backfit = TRUE, nullcheck = TRUE, verbose = 1L
)
# L >= 0 (all flashier plots work); F signed sparse
```

### GBCD

```r
library(gbcd)
set.seed(99)
fit <- fit_gbcd(Y = X, Kmax = 15, prior = ebnm::ebnm_generalized_binary,
                maxiter1 = 500, maxiter2 = 200, maxiter3 = 500, verbose = 1)
# Results: fit$L (n×K), fit$F$lfc (p×K), fit$F$lfsr (p×K)
```
