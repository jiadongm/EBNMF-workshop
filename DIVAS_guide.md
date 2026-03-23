# DIVAS Practical Guide: Data Integration Via Analysis of Subspaces

> Reference for future Claude Code sessions.
> All patterns confirmed on the COVID-19 proteomics CovidCase project (2026-03).
> DIVAS v0.1.1. Script: `CovidCase/08_divas.R`.

---

## 1. Package and installation

```r
# DIVAS is not on CRAN — install from GitHub
remotes::install_github("reagan-lee-unc/DIVAS")

library(DIVAS)
```

---

## 2. When to use DIVAS vs NMF stacking

### NMF stacking limitations

When data has multiple blocks (time points, modalities), NMF can be run by stacking
the data in two ways:

- **Sample stacking** (rows): stack T1 and T2 observations → 2N × p matrix. Factors
  that load similarly at T1 and T2 are stable programs; factors that differ capture
  time-dependent change. But you cannot distinguish stable from time-specific factors
  without post-hoc analysis.

- **Feature stacking** (columns): concatenate columns across modalities → N × (p₁ + p₂).
  A factor spanning both captures coordinated programs. But interpretation requires
  knowing which columns came from which modality.

### What DIVAS adds

DIVAS finds which components are:
- **Jointly shared**: present in all blocks with the same variation across samples
- **Partially shared**: present in a subset of blocks (e.g., T1 and T2 but not metabolomics)
- **Individual**: unique to one block

Unlike stacking, DIVAS makes this structural inference **statistically formal** via
rotational bootstrap (nsim simulations per block). The number of shared/individual
components is inferred from the data, not fixed in advance.

---

## 3. Data preparation

### 3.1 Input orientation

```r
# DIVAS convention: p_features × n_samples per block
# This is OPPOSITE to flashier (which uses n × p)

# Example (COVID-19 proteomics):
# Block T1: 481 proteins × 120 patients  ← features in rows, samples in columns
# Block T2: 481 proteins × 120 patients

dim(Prot_T1)  # [481, 120]
dim(Prot_T2)  # [481, 120]
```

### 3.2 Column alignment (critical)

All blocks must have **identical column order** — same samples in the same order.
DIVAS does not match samples by name; it assumes columns are aligned.

```r
# Ensure patient_id order matches across blocks
t1_meta <- meta[meta$time_point == "T1", ] |> arrange(patient_id)
t2_meta <- meta[meta$time_point == "T2", ] |> arrange(patient_id)
stopifnot(all(t1_meta$patient_id == t2_meta$patient_id))   # must be TRUE

patient_ids <- as.character(t1_meta$patient_id)
Prot_T1 <- as.matrix(prot_raw[, t1_meta$sample])
Prot_T2 <- as.matrix(prot_raw[, t2_meta$sample])
colnames(Prot_T1) <- colnames(Prot_T2) <- patient_ids     # set to patient IDs
```

### 3.3 Feature centering (rowCent)

```r
# rowCent=TRUE: centre features (rows) within each block
# This removes the per-feature mean within each block, so DIVAS finds
# variation across samples rather than global feature level differences.
# Recommended in almost all analyses.

# colCent=FALSE: do NOT centre samples (columns)
# Centering samples would remove sample-level effects (e.g., batch).
# Keep FALSE to preserve between-sample variation.
```

### 3.4 Longitudinal (same assay, multiple time points)

```r
# T1 and T2 of the same assay
data_list <- list(
  T1 = Prot_T1,   # 481 proteins × 120 patients (features × samples)
  T2 = Prot_T2    # same orientation, same patient column order
)
# Component naming: "T1+T2-k" (shared), "T1-Individual-k", "T2-Individual-k"
```

### 3.5 Multiomics (different assays, same samples)

```r
# Different assays: proteins, metabolites, genes
data_list <- list(
  Prot  = Prot_mat,    # 481 × N
  Metab = Metab_mat,   # 763 × N
  GEX   = GEX_mat      # 20091 × N
)
# All blocks must share the same N samples in the same column order.
# rowCent=TRUE especially important here (different feature units/scales).
# Component naming generalises: "Prot+Metab+GEX-k" (fully shared),
# "Prot+Metab-k" (partially shared), "Prot-Individual-k" (block-specific).
```

---

## 4. Fitting: DIVASmain()

### Full function signature

```r
divas_res <- DIVASmain(
  datablock    = data_list,   # named list of (p_k × N) matrices
  nsim         = 400,         # bootstrap simulations for rank estimation
  colCent      = FALSE,       # centre samples (usually FALSE)
  rowCent      = TRUE,        # centre features within each block (usually TRUE)
  seed         = 123,         # reproducibility
  ReturnDetail = TRUE         # MUST be TRUE to use getTopFeatures()
)
```

### Parameters

- **`nsim`**: Number of bootstrap simulations for signal rank estimation.
  400 is reliable. Fewer simulations (e.g., 100) are faster but may give unstable
  rank estimates. Use 400 for publication-quality results.

- **`ReturnDetail=TRUE`**: Required to access loadings via `getTopFeatures()`. If
  `ReturnDetail=FALSE`, the `$Loadings` field will be absent and `getTopFeatures()`
  will fail.

- **Auto K**: DIVAS has no `Kmax` parameter. Signal rank is estimated per block via
  bootstrap, then the shared/individual decomposition is performed. K is fully
  data-driven.

### Runtime

The COVID-19 2-block proteomics analysis (481 proteins × 120 patients, nsim=400)
completed in approximately 10–20 minutes. Runtime scales with nsim and block size.
Always cache results.

```r
if (!file.exists("results/divas_result.rds")) {
  set.seed(123)
  divas_res <- DIVASmain(datablock = data_list, nsim = 400, colCent = FALSE,
                          rowCent = TRUE, seed = 123, ReturnDetail = TRUE)
  saveRDS(divas_res, "results/divas_result.rds")
} else {
  divas_res <- readRDS("results/divas_result.rds")
}
```

---

## 5. Accessing results

### 5.1 Component naming convention

The naming convention encodes which blocks share which component:

| Name pattern | Meaning |
|-------------|---------|
| `"T1+T2-k"` | Shared between T1 and T2; k-th such component |
| `"T1-Individual-k"` | Individual to T1 only |
| `"T2-Individual-k"` | Individual to T2 only |
| `"Prot+Metab-k"` | Shared between Prot and Metab but not other blocks |
| `"Prot-Individual-k"` | Individual to Prot only |

Parse component names programmatically:

```r
comp_names <- colnames(divas_res$Scores)

# Categorise
comp_df <- data.frame(component = comp_names, stringsAsFactors = FALSE)
comp_df$category <- sapply(comp_names, function(nm) {
  if (grepl("\\+", nm))         "Shared"       # has '+' separator
  else if (grepl("-Individual", nm)) "Individual"
  else                           "Other"
})
table(comp_df$category)
```

### 5.2 Scores matrix

```r
scores_mat <- divas_res$Scores
# Dimensions: N_samples × K_total
# Rows = samples (patients, cells), columns = DIVAS components
# Colnames = component names (e.g., "T1+T2-2", "T1-Individual-3")

dim(scores_mat)  # [120, 24] in the COVID-19 2-block example
```

### 5.3 Loadings

```r
# Only available if ReturnDetail=TRUE
# $Loadings is a list indexed by component name (INTERNAL naming — see §6 gotcha)
# Dimensions per component: p_k × r_i (features × rank for that component)
names(divas_res$Loadings)   # list of internal component names
```

---

## 6. Extracting top features: getTopFeatures()

### Function signature

```r
feats <- getTopFeatures(
  divasRes  = divas_res,
  compName  = "T1+T2-2",   # component name (see gotcha below)
  modName   = "T1",         # block name EXACTLY as in data_list
  n_top_pos = 10,           # number of top positive-loading features
  n_top_neg = 10            # number of top negative-loading features
)
feats$top_positive   # character vector of feature names
feats$top_negative   # character vector of feature names
```

### Critical gotcha: individual component name remapping

**The Scores matrix and `getTopFeatures()` use DIFFERENT names for individual components.**

| Scores column name | getTopFeatures compName |
|-------------------|------------------------|
| `"T1-Individual-3"` | `"T1-3"` |
| `"T2-Individual-1"` | `"T2-1"` |
| `"T1+T2-2"` | `"T1+T2-2"` (same) |

For shared components, the name is the same in both. For individual components, strip
`"-Individual"` before passing to `getTopFeatures()`.

### Confirmed working pattern

```r
# Helper function to build the right (compName, modName) queries
get_feature_queries <- function(comp_name, block_names) {
  if (grepl("\\+", comp_name)) {
    # Shared component: query each block with the same comp_name
    lapply(block_names, function(b) {
      if (grepl(b, comp_name)) list(comp = comp_name, mod = b)
      else NULL
    }) |> Filter(Negate(is.null), x = _)
  } else {
    # Individual component: strip "-Individual-" → "Block-k"
    # e.g. "T1-Individual-3" → comp="T1-3", mod="T1"
    for (b in block_names) {
      pattern <- paste0("^", b, "-Individual-")
      if (grepl(pattern, comp_name)) {
        k <- sub(pattern, "", comp_name)
        return(list(list(comp = paste0(b, "-", k), mod = b)))
      }
    }
    list()
  }
}

# Example usage for 2-block T1/T2 longitudinal data
block_names <- names(data_list)   # c("T1", "T2")

for (comp in colnames(divas_res$Scores)) {
  queries <- get_feature_queries(comp, block_names)
  for (q in queries) {
    tryCatch({
      feats <- getTopFeatures(divas_res, compName = q$comp, modName = q$mod,
                              n_top_pos = 10, n_top_neg = 10)
      cat(sprintf("%s | %s\n", comp, q$mod))
      cat("  +:", paste(head(feats$top_positive, 5), collapse = "; "), "\n")
      cat("  -:", paste(head(feats$top_negative, 5), collapse = "; "), "\n")
    }, error = function(e) {
      cat(sprintf("  [%s | %s FAILED: %s]\n", comp, q$mod, conditionMessage(e)))
    })
  }
}
```

### Simplified pattern for 2-block T1/T2 (exact confirmed code from 08_divas.R)

```r
# This pattern was tested and confirmed working
get_feature_queries <- function(comp_name) {
  if (grepl("^T1\\+T2-", comp_name)) {
    # Shared: query both modalities directly
    list(list(comp = comp_name, mod = "T1"),
         list(comp = comp_name, mod = "T2"))
  } else if (grepl("^T1-Individual-", comp_name)) {
    k <- sub("^T1-Individual-", "", comp_name)
    list(list(comp = paste0("T1-", k), mod = "T1"))
  } else if (grepl("^T2-Individual-", comp_name)) {
    k <- sub("^T2-Individual-", "", comp_name)
    list(list(comp = paste0("T2-", k), mod = "T2"))
  } else {
    list()
  }
}
```

---

## 7. Diagnostic plot: DJIVEAngleDiagnosticJP()

```r
# Produces the DIVAS angle diagnostic — shows which components are signal vs noise.
# Must write to a graphics device (pdf/png) — does not return a ggplot object.

pdf("results/figs/divas_diagnostic.pdf", width = 10, height = 8)
DJIVEAngleDiagnosticJP(
  datablock = data_list,                  # same data_list used for DIVASmain
  dataname  = names(data_list),           # block names: c("T1", "T2")
  outstruct = divas_res,                  # DIVASmain result
  randseed  = 123,
  titlestr  = "DIVAS: Proteomics T1 vs T2"
)
dev.off()
```

---

## 8. Associating components with metadata

```r
scores_mat <- divas_res$Scores   # N × K
comp_names <- colnames(scores_mat)
K_total    <- ncol(scores_mat)

# Spearman correlation of each component score with a metadata variable
severity_vec <- ...   # N-vector, one value per sample

cor_vals <- apply(scores_mat, 2, function(x)
  cor(x, severity_vec, method = "spearman", use = "complete.obs"))

pval_vals <- sapply(seq_len(K_total), function(k)
  cor.test(scores_mat[, k], severity_vec, method = "spearman")$p.value)

assoc_df <- data.frame(
  component  = comp_names,
  rho        = round(cor_vals, 3),
  pval       = round(pval_vals, 4),
  fdr        = round(p.adjust(pval_vals, method = "BH"), 4)
)
assoc_df <- assoc_df[order(abs(assoc_df$rho), decreasing = TRUE), ]
print(assoc_df)
```

### Interpretation guidance

From the COVID-19 2-block proteomics analysis:

- **Shared components (T1+T2-k)** tend to carry the strongest severity associations.
  The top shared component (T1+T2-2) had rho = 0.730 with T1 severity (FDR < 0.0001) —
  a time-stable severity program captured at both admission and follow-up.

- **Individual components** typically have weaker severity associations. This makes
  biological sense: severity-related variation is stable across time points, so it is
  captured by shared components. Individual components capture time-specific biology
  (e.g., early inflammatory response vs late resolution).

- **WHO delta (T2−T1) associations** are strongest for shared components, not individual
  ones. T1+T2-7 had rho = −0.344 with delta (FDR = 0.003) — patients with higher scores
  improved more.

---

## 9. Known gotchas (consolidated)

1. **Individual component name remapping.** `divas_res$Scores` column names use
   `"T1-Individual-k"`, but `getTopFeatures()` requires `"T1-k"`. Always parse with
   the helper function in §6 — do not pass Scores column names to `getTopFeatures()`
   directly for individual components.

2. **Input orientation.** DIVAS expects **p_features × n_samples** per block. This is
   OPPOSITE to flashier (which uses n × p). Always verify `dim()` before fitting.

3. **Column alignment is the user's responsibility.** DIVAS does not check or enforce
   that columns represent the same samples across blocks. Misaligned columns will produce
   wrong results silently. Always `stopifnot(all(block1_patients == block2_patients))`
   before building `data_list`.

4. **`ReturnDetail=TRUE` is required** for `getTopFeatures()`. If omitted, the `$Loadings`
   field will be absent and `getTopFeatures()` will throw an error.

5. **`dataname` in diagnostic must match `names(data_list)`.** Pass the block names
   exactly as they appear in the list to `DJIVEAngleDiagnosticJP()`.

6. **Component naming generalises.** For >2 blocks, DIVAS produces component names like
   `"Prot+Metab+GEX-k"` (jointly shared across all three), `"Prot+Metab-k"` (partially
   shared), and `"Prot-Individual-k"` (block-specific). Parse by splitting on `"+"`.

7. **K is determined by bootstrap.** The number of shared and individual components
   depends on the signal-to-noise in each block (estimated via `nsim` permutations per
   block). Results may differ slightly across runs if `seed` is not fixed. Fix `seed`
   for reproducibility.

8. **Score column name parsing.** To extract which blocks are in a shared component:
   ```r
   # "Prot+Metab+GEX-2" → c("Prot", "Metab", "GEX")
   nm <- "Prot+Metab+GEX-2"
   blocks_in_comp <- strsplit(sub("-[0-9]+$", "", nm), "\\+")[[1]]
   ```

---

## 10. Multi-block generalisation (>2 blocks)

```r
# Same API; data_list can have K > 2 entries
data_list <- list(
  Prot  = Prot_mat,    # p1 × N
  Metab = Metab_mat,   # p2 × N
  GEX   = GEX_mat      # p3 × N
)

divas_res <- DIVASmain(
  datablock    = data_list,
  nsim         = 400,
  colCent      = FALSE,
  rowCent      = TRUE,
  seed         = 123,
  ReturnDetail = TRUE
)

# Component naming generalises automatically:
# "Prot+Metab+GEX-k"   ← fully shared
# "Prot+Metab-k"       ← partially shared (Prot and Metab only)
# "Prot-Individual-k"  ← individual to Prot
```

### getTopFeatures for multi-block

```r
# For shared component "Prot+Metab-2", query Prot and Metab separately:
feats_prot  <- getTopFeatures(divas_res, compName = "Prot+Metab-2", modName = "Prot",
                               n_top_pos = 10, n_top_neg = 10)
feats_metab <- getTopFeatures(divas_res, compName = "Prot+Metab-2", modName = "Metab",
                               n_top_pos = 10, n_top_neg = 10)

# For individual "Prot-Individual-3", use the remapped name "Prot-3":
feats_prot_indiv <- getTopFeatures(divas_res, compName = "Prot-3", modName = "Prot",
                                    n_top_pos = 10, n_top_neg = 10)
# Pattern: strip "-Individual" from the block-individual prefix
```

---

## 11. Quick-start template

```r
library(DIVAS)

# 1. Build data_list (features × samples per block, aligned columns)
data_list <- list(T1 = Prot_T1, T2 = Prot_T2)

# 2. Fit (cache to disk)
if (!file.exists("results/divas.rds")) {
  set.seed(123)
  divas_res <- DIVASmain(datablock = data_list, nsim = 400, colCent = FALSE,
                          rowCent = TRUE, seed = 123, ReturnDetail = TRUE)
  saveRDS(divas_res, "results/divas.rds")
} else {
  divas_res <- readRDS("results/divas.rds")
}

# 3. Inspect
scores_mat <- divas_res$Scores
comp_names <- colnames(scores_mat)
cat(sprintf("Components: %d\n", ncol(scores_mat)))
print(table(
  sapply(comp_names, function(nm)
    if (grepl("\\+", nm)) "Shared" else "Individual")
))

# 4. Associate with metadata
rho_sev <- apply(scores_mat, 2, function(x) cor(x, t1_severity, method = "spearman"))
print(sort(abs(rho_sev), decreasing = TRUE))

# 5. Top features for strongest component
top_comp <- comp_names[which.max(abs(rho_sev))]
feats <- getTopFeatures(divas_res, compName = top_comp, modName = "T1",
                         n_top_pos = 20, n_top_neg = 20)
```
