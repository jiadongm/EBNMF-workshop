# CLAUDE.md — CovidCase

> Session notes for the `CovidCase/` directory.
> Goal: apply EBNMF tools to the COVID-19 multiomics cohort from Su et al. (Cell 2020).

---

## Dataset overview

**Source:** Su et al. (2020). "Multi-Omics Resolves a Sharp Disease-State Shift between Mild
and Moderate COVID-19." *Cell* **183**, 1479–1495.

**Cohort:** 120 COVID-19 patients measured at two time points (T1 ≈ days 1–4, T2 ≈ days 7–10
post admission). Severity scored on WHO ordinal scale 1–7.

**Sample naming:** `COVID_[patient_id]_T[1|2]` — 240 samples total (120 × 2 time points).

### Data files (in `data/240_input_4omics/`)

| File | Rows (features) | Columns (samples) | Description |
|------|----------------|-------------------|-------------|
| `proteomics_120patients.csv` | 481 proteins | 240 | Olink NPX, log2 scale |
| `metabolomics_120patients.csv` | 763 metabolites | 240 | Metabolon, median-scaled ratios |
| `sc_pro_120patients.csv` | 192 surface proteins | 240 | CITE-seq pseudobulk, log scale |
| `sc_pro_120patients_clr_normalized.csv` | 192 surface proteins | 240 | CLR-normalised version |
| `sc_gex_120patients_aligned.csv` | 20,091 genes | 240 | scRNA-seq pseudobulk, mean ln(CPM+1) |

All files share the same 240 sample columns in the same order (verified).

---

## Proteomics data: platform and normalisation

**Platform:** Olink Proximity Extension Assay (PEA), five panels (Cardiovascular II,
Inflammation, Metabolism, Immune Response, Organ Damage). 481 proteins after QC.

**Output unit:** NPX (Normalized Protein eXpression) — Olink's proprietary log2-scale unit.
Computed from raw Ct values by:
1. Normalisation against an extension control
2. Normalisation against an inter-plate control
3. Batch correction via pooled overlapping reference samples across plates

NPX can be negative (protein near or below the limit of detection). This is biologically
valid and not a data error.

**Normalisation decision:** Use NPX values as-is. The data file already contains
analysis-ready log2-normalised values per the paper's STAR Methods. No further
transformation is applied in `prepareData.R`.

**Note on NMF:** Point-exponential prior in flashier requires non-negative input. When
running NMF, a per-protein minimum-shift (`X[,j] - min(X[,j])`) will be applied in the
analysis script. This is kept separate from data preparation to preserve original NPX values.

---

## Scripts

### `prepareData.R` — Status: done

Loads and QC-checks the proteomics data. Run from `CovidCase/`:

```r
Rscript prepareData.R
```

**What it does:**
1. Loads `proteomics_120patients.csv` (481 proteins × 240 samples)
2. Parses sample metadata (patient_id, time_point) from column names
3. Computes per-protein QC statistics (mean, SD, min, max, n_negative)
4. Produces three QC figures
5. Transposes to samples × proteins (240 × 481) and saves as `results/X_prot.rds`

**Key QC findings:**
- No missing values
- NPX range: −3.775 to 17.570 (log2)
- 83/481 proteins have at least one negative NPX (near LOD); none are entirely negative
- Per-sample mean NPX: 4.59 ± 0.08 — very stable, no outlier samples

**Outputs:**

| File | Description |
|------|-------------|
| `results/X_prot.rds` | 240 × 481 matrix, samples × proteins, NPX log2 |
| `results/sample_meta.rds` | Data frame: sample, patient_id, time_point |
| `results/prot_qc_summary.csv` | Per-protein mean, SD, min, max, n_neg |
| `results/figs/prot_npx_distribution.pdf` | Violin plot of NPX by time point |
| `results/figs/prot_sample_means.pdf` | Per-sample mean NPX (sorted) |
| `results/figs/prot_mean_sd.pdf` | Mean vs SD scatter per protein |

---

## Scripts status

| Script | Status | Description |
|--------|--------|-------------|
| `prepareData.R` | ✅ Done | Load, QC, save X_prot.rds + sample_meta.rds |
| `01_PCA.R` | ✅ Done | PCA baseline, trajectories, save pca_scores.rds |
| `02_NMF_stacked.R` | ✅ Done | NMF on stacked 240×481 (S1 NMF K=8, S2 EBNMF auto K=14, S3 GBCD K=10) |
| `03_NMF_delta.R` | ✅ Done | NMF on delta 120×481 (D1 NMF K=6, D2 Laplace K=20, D3 Normal K=20, D4 GBCD K=18, D5 Asym EBMF K=20) |
| `04_loading_heatmaps.R` | ✅ Done | Pheatmap of L matrices for S1/S2/S3; style mirrors GBCD HNSCC vignette |
| `05_factor_plots.R` | ✅ Done | F-matrix pheatmaps + flashier vignette-style plots (scree, histogram, structure, heatmap, scatter, bar) with two patient colourings |
| `06_delta_loading_heatmaps.R` | ✅ Done | Pheatmap of L matrices for D1–D5; diverging scale for signed; non-neg for D1/D5; 6 PDFs (5 individual + all5 combined 30×9 in) |
| `07_delta_factor_plots.R` | ✅ Done | F-matrix pheatmaps for D1–D5 + flashier plots for D2 (limited: signed L) and D5 (full: non-neg L); D5 enables structure/L-heatmap/bar |
| `08_divas.R` | ✅ Done | DIVAS 2-block (proteomics T1 vs T2); 10 shared + 8 T1-individual + 6 T2-individual components |

## Next steps

- [x] Run `03_NMF_delta.R` — primary analysis (D2 point-Laplace on T2-T1 delta)
- [x] DIVAS analysis — `08_divas.R` complete (2 blocks: proteomics T1 vs T2)
- [ ] Extend `prepareData.R` to load and normalise sc and metabolomics data
- [ ] Multi-block DIVAS (add metabolomics, sc_pro, sc_gex blocks)

---

## Results: 02_NMF_stacked.R (stacked 240×481 matrix)

**Run date:** 2026-02-28. All three variants completed successfully.

**Input:** `X_shifted` (240 × 481), per-protein min-shift of NPX log2 values → range [0, 12.18].
**Association outcome:** Spearman rho(L_T2 − L_T1, who_delta) per factor (n=69 patients with
non-NA who_delta; 51 excluded due to ambiguous WHO scores at T1 or T2).

### Runtimes

| Variant | K | Runtime |
|---------|---|---------|
| S1 Standard NMF (brunet, 10 runs) | 8 (fixed) | 31.5 s |
| S2 EBNMF (point-exponential, auto K) | 14 (auto) | 39.9 s |
| S3 GBCD (generalised binary, Kmax=15) | 10 (auto) | 471 s |

### S1 — Standard NMF (K=8, 31.5 s)

No factor reaches BH FDR < 0.2 for who_delta association. Top associations:
- **F1:** rho = +0.262, p = 0.030, FDR = 0.238
- **F4:** rho = −0.203, p = 0.095, FDR = 0.378
- **F2:** rho = −0.173, p = 0.156, FDR = 0.416

### S2 — EBNMF point-exponential (auto K=14, 39.9 s)

Auto-selected **K=14** (greedy stopped at factor 15; backfit hit max iterations — fit
is exploratory but PVE sum = 0.90 shows good coverage).

**PVE per factor:** 18.5, 7.4, 8.7, 30.2, 1.4, 1.5, 6.2, 4.5, 2.1, 0.7, 0.3, 2.6, 1.2, 4.7 %

Note: F4 dominates with 30.2% PVE — likely a global severity/time-point axis.
F1 (18.5%) is the second largest; together F1+F4 account for ~49% of variance.

Strongest who_delta associations (best across all three variants):
- **F5:** rho = −0.344, p = 0.004, FDR = 0.054 — near-significant after correction
- **F1:** rho = −0.273, p = 0.023, FDR = 0.161
- **F8:** rho = +0.226, p = 0.062, FDR = 0.257
- **F2:** rho = +0.209, p = 0.085, FDR = 0.257
- **F7:** rho = +0.205, p = 0.092, FDR = 0.257

Negative rho for F5 and F1: patients with higher loading change (T2−T1 increase) in these
factors tended to **improve** (more negative who_delta). Interpretation pending protein-level
inspection of F5 and F1 gene weights.

### S3 — GBCD (auto K=10, 471 s)

Auto-selected **K=10** out of a possible maximum of 29 (2×Kmax−1). Runtime dominated by the
three-stage GBCD optimisation (greedy flash 418 s + signature estimation 2 s).

No factor reaches BH FDR < 0.2. Top associations:
- **F10:** rho = +0.256, p = 0.034, FDR = 0.240
- **F8:**  rho = −0.232, p = 0.055, FDR = 0.240
- **F6:**  rho = +0.218, p = 0.072, FDR = 0.240

GBCD gives binary-like memberships (L); the weaker association signal vs S2 may reflect
the harder thresholding of the GB prior. The LFC-scale F matrix (F$lfc) gives directly
interpretable log-fold changes per protein per program, with per-protein lfsr (F$lfsr).

### S2 — ebnm_fn single-function behaviour and sparsity

**`ebnm_fn` with a single function applies to both L and F.** The same prior family is
used for both matrices. To use different priors per matrix, pass a two-element list:
`ebnm_fn = list(fn_for_L, fn_for_F)` — first element → L, second → F.

Consequence here: `ebnm_point_exponential` enforces non-negativity on both patient
loadings (L) and protein weights (F), which is the correct NMF constraint.

**Sparsity of S2 fit** (near-zero defined as |value| < 0.01):

*L matrix (240 samples × 14 factors) — overall 9.8% sparse — very dense:*

| Factor | L sparsity | Notes |
|--------|-----------|-------|
| F06 | 35.4% | ~85 patients absent |
| F07 | 22.1% | |
| F14 | 21.7% | |
| F03, F11 | 10.4% | |
| F08, F10 | 0.0% | every patient loads on these |
| remaining | < 6% | nearly all patients present |

Most patients contribute to most factors. Point-exponential shrinks small loadings but
does not zero them aggressively when the proteomics signal is dense across samples.

*F matrix (481 proteins × 14 factors) — overall 40.7% sparse — much sparser than L:*

| Factor | F sparsity | Active proteins (~) | who_delta rho |
|--------|-----------|---------------------|---------------|
| F08 | **98.8%** | ~6 | +0.226, FDR=0.257 |
| F11 | **97.5%** | ~12 | — |
| F10 | **96.9%** | ~15 | — |
| F05 | 71.9% | ~135 | **−0.344, FDR=0.054** |
| F02 | 74.0% | ~125 | +0.209, FDR=0.257 |
| F04, F09, F12 | ~35% | ~310 | — |
| F01, F03, F13 | <5% | ~460 | F01: −0.273, FDR=0.161 |

Key asymmetry: the prior learned far more sparsity in F (protein profiles) than in L
(patient memberships). This is biologically sensible — programs should be defined by
subsets of proteins, but patients can participate in multiple programs simultaneously.

The most clinically associated factors (F5, F8) are also among the sparsest in F —
focused protein modules rather than diffuse background covariation. F1 is a notable
exception: near-dense in F (almost all 481 proteins) yet has rho=−0.273, suggesting
a global proteome-wide shift that is still correlated with improvement.

### Cross-variant comparison

| Variant | K | Best |rho| | Best FDR |
|---------|---|-----------|----------|
| S1 NMF  | 8  | 0.262 (F1) | 0.238 |
| S2 EBNMF | 14 | 0.344 (F5) | **0.054** |
| S3 GBCD | 10 | 0.256 (F10) | 0.240 |

S2 EBNMF gives the strongest who_delta signal. **S2-F5** (rho=−0.344, FDR=0.054) is the
primary candidate for a proteomics program associated with clinical improvement T1→T2.

### Loading heatmaps (04_loading_heatmaps.R)

Style mirrors the GBCD HNSCC vignette (`pheatmap`, gray96 → red, [0,1] scale):
- **Rows:** 240 samples ordered by time_point (T1 first) then WHO severity ascending
- **Columns:** programs ordered by PVE for S2; as-fit for S1 and S3
- **Row annotation:** time_point (blue/red) + WHO severity (white → dark blue)
- **Column annotation:** |Spearman rho| with who_delta (white → orange)
- **Normalisation:** per-column max-normalisation to [0,1] for S1 and S2;
  S3 GBCD L is already in [0,1] by construction

Outputs:
- `results/figs/nmf/heatmap_S1_loadings.pdf` — S1 (6×9 in)
- `results/figs/nmf/heatmap_S2_loadings.pdf` — S2 (8×9 in)
- `results/figs/nmf/heatmap_S3_loadings.pdf` — S3 (7×9 in)
- `results/figs/nmf/heatmap_all3_loadings.pdf` — all three side-by-side (22×9 in)

## Results: 03_NMF_delta.R (delta 120×481 matrix, T2-T1)

**Run date:** 2026-02-28. All four variants completed successfully.

**Input:** `delta_mat` (120 × 481), T2 − T1 NPX per patient, signed. Range: [−6.19, 7.28].
**Association:** Spearman rho(L_k, who_delta) directly (delta rows = patients; n=69 non-NA).

### Runtimes

| Variant | K | Runtime |
|---------|---|---------|
| D1 Standard NMF (brunet, 10 runs, shifted) | 6 (fixed) | 12.8 s |
| D2 EBMF (point-Laplace, auto K) | 20 (auto) | 28.6 s |
| D3 EBMF (point-normal, auto K) | 20 (auto) | 13.0 s |
| D4 GBCD (generalised binary, Kmax=10) | 18 (auto) | ~4 min |
| D5 Asymmetric EBMF (point-exp L, point-Laplace F) | 20 (auto) | 57 s |

### D1 — Standard NMF (K=6, shifted delta)

One factor reaches BH FDR < 0.05:
- **F2:** rho = −0.316, p = 0.008, **FDR = 0.049** ✓

### D2 — EBMF point-Laplace (auto K=20, PRIMARY)

PVE sum = 0.724. PVE per factor: 13.5, 10.0, 3.6, 9.6, 5.2, 8.4, 2.5, 2.8, ...
No factor reaches BH FDR < 0.2. Best associations:
- **F17:** rho = +0.279, p = 0.020, FDR = 0.213
- **F5:** rho = +0.259, p = 0.032, FDR = 0.213
- **F2:** rho = +0.250, p = 0.038, FDR = 0.213
- **F4:** rho = −0.245, p = 0.043, FDR = 0.213

### D3 — EBMF point-normal (auto K=20)

PVE sum = 0.721. No factor reaches BH FDR < 0.2. Best:
- **F16:** rho = +0.263, p = 0.029, FDR = 0.309
- **F5:** rho = +0.260, p = 0.031, FDR = 0.309

### D5 — Asymmetric EBMF (point-exponential L, point-Laplace F; auto K=20)

Non-negative L (range [0.001, 4.489]) with signed sparse F (range [−2.896, 2.221]).
PVE sum = 0.840 — higher than D2 (0.724) or D3 (0.721) despite non-negativity constraint
on L, because non-negativity forces F to absorb all directionality more efficiently.

No factor reaches BH FDR < 0.2. Best association:
- **F1:** rho = −0.256, p = 0.034, FDR = 0.380

**Key advantage over D2/D3:** Non-negative L enables all flashier plot types (structure,
L heatmap, bar) which are blocked by signed L. Patient participation scores ≥ 0 give a
cleaner interpretation: each patient "uses" each program to a non-negative degree.

**ebnm_fn semantics:** `ebnm_fn = list(ebnm_point_exponential, ebnm_point_laplace)` —
first element applies to L, second to F. A single function applies to both.

### D4 — GBCD (auto K=18, binary-like memberships)

**Strongest association signal across all delta variants:**
- **F15:** rho = −0.363, p = 0.002, **FDR = 0.039** ✓
- **F3:** rho = −0.307, p = 0.010, FDR = 0.093
- **F4:** rho = −0.274, p = 0.023, FDR = 0.136
- **F8:** rho = −0.251, p = 0.038, FDR = 0.170
- **F13:** rho = −0.234, p = 0.053, FDR = 0.187
- **F2:** rho = −0.226, p = 0.062, FDR = 0.187

6 factors achieve |rho| > 0.2 with FDR < 0.2 — best performance on the delta matrix.
All top associations are **negative rho**: higher factor loading in delta = more improvement (lower who_delta).

**Note on D4 L range:** L values range from −0.54 to 1.0 (slight negatives expected when
the input data is signed; the GB prior pushes L toward 0/1 but can't guarantee [0,1] with
signed data). Factor memberships near 0 and 1 are still biologically interpretable.

### Cross-variant comparison (delta matrix)

| Variant | K | prior_L | prior_F | Best \|rho\| | Best FDR |
|---------|---|---------|---------|------------|----------|
| D1 NMF  | 6  | non-neg (shifted) | non-neg (shifted) | 0.316 (F2) | **0.049** |
| D2 EBMF Laplace | 20 | point-Laplace | point-Laplace | 0.279 (F17) | 0.213 |
| D3 EBMF Normal | 20 | point-normal | point-normal | 0.263 (F16) | 0.309 |
| D4 GBCD | 18 | gen. binary | signed LFC | 0.363 (F15) | **0.039** |
| D5 Asym EBMF | 20 | point-exponential | point-Laplace | 0.256 (F1) | 0.380 |

**D4 GBCD gives the strongest signal on the delta matrix** (F15 FDR=0.039).
D1 NMF F2 (FDR=0.049) is second. Both reach FDR < 0.05.
D2/D3/D5 do not reach FDR < 0.2 (delta signal is diffuse in continuous models).

**Key insight:** GBCD's binary-like prior is better suited to the delta matrix than
point-Laplace or point-normal — it identifies patient subgroups whose proteomics
changed together, rather than finding a continuous loading gradient.
D5's value is interpretability (non-neg L) rather than association strength.

### Saved objects (delta)

| File | Description |
|------|-------------|
| `results/nmf/delta_D1_nmf.rds` | NMFfit (shifted delta, K=6) |
| `results/nmf/delta_D2_laplace.rds` | flashier fit (point-Laplace, K=20) |
| `results/nmf/delta_D3_normal.rds` | flashier fit (point-normal, K=20) |
| `results/nmf/delta_D4_gbcd.rds` | gbcd fit (L, F$lfc, F$lfsr, K=18) |
| `results/nmf/delta_D5_asym_ebmf.rds` | flashier fit (asym: point-exp L, point-Laplace F, K=20) |
| `results/nmf/delta_association.csv` | Spearman rho + BH FDR, all 5 variants × factors |

### Delta loading heatmaps (06_delta_loading_heatmaps.R)

**Run date:** 2026-02-28. 6 PDFs produced successfully.

- Rows = 120 patients, ordered by who_delta_cat then T1 severity ascending
- Row annotation: who_delta_cat (Improved/Stable/Worsened/Ambiguous) + T1_severity
- Column annotation: |Spearman rho| with who_delta (white → orange)
- **D1 (non-neg shifted):** col-norm [0,1], gray96→red
- **D2/D3 (signed L):** diverging blue→white→red; symmetric ±max
- **D4 GBCD (binary-like):** col-norm [0,1], gray96→red
- **D5 (non-neg L):** col-norm [0,1], gray96→red; columns in PVE order

| File | Description |
|------|-------------|
| `results/figs/nmf/heatmap_delta_D1_loadings.pdf` | D1 NMF loading heatmap (8×8 in) |
| `results/figs/nmf/heatmap_delta_D2_loadings.pdf` | D2 EBMF Laplace heatmap, diverging (10×8 in) |
| `results/figs/nmf/heatmap_delta_D3_loadings.pdf` | D3 EBMF Normal heatmap, diverging (10×8 in) |
| `results/figs/nmf/heatmap_delta_D4_loadings.pdf` | D4 GBCD heatmap (10×8 in) |
| `results/figs/nmf/heatmap_delta_D5_loadings.pdf` | D5 Asym EBMF heatmap (10×8 in) |
| `results/figs/nmf/heatmap_delta_all5_loadings.pdf` | All 5 side-by-side (30×9 in) |

### Delta factor plots (07_delta_factor_plots.R)

**Run date:** 2026-02-28. Successfully produced PDFs for D2 (limited) and D5 (full).

**Patient colourings:** who_delta_cat (Improved=steelblue, Stable=gold3, Worsened=firebrick,
Ambiguous=gray70) and T1_severity (who_score_num at T1; white→dark blue gradient).

**Section A — F-matrix pheatmaps (protein × factor):**
- D1/D2/D3: diverging scale where F signed; D1 non-neg col-norm; D4 GBCD LFC diverging
- D5: diverging scale (signed F); rows clustered ward.D2; column annotation = |rho|
- Files: `heatmap_F_D1.pdf` through `heatmap_F_D5.pdf`

**Section B — flashier plots for D2 (point-Laplace, signed L):**
- Works: scree, histogram × 2 (delta_cat, severity), scatter × 3 (F17, F5, F2), heatmap factors (no pm_groups)
- Skipped: structure plot, L heatmap, factors heatmap + pm_groups (all fail: signed L → verify.nonnegative.matrix)
- Bar plots for D2 (loadings): also work

**Section D — flashier plots for D5 (asymmetric, non-negative L):**
- ALL plot types work including structure, L heatmap, bar
- d5_D3a/b: structure × delta_cat and × severity (enabled by non-neg L)
- d5_D4a/b: L heatmap × delta_cat and × severity (enabled by non-neg L)
- d5_D5: top-protein F heatmap (kset = c(1,12,18,4) — most-associated D5 factors)
- d5_D6: scatter for F1, F12, F18
- d5_D7a/b: bar loadings × delta_cat and × severity

---

### Factor plots (05_factor_plots.R)

**Run date:** 2026-02-28. All 16 PDFs produced successfully.

Three blocks of visualisations for the S2 EBNMF fit (K=14) plus F-matrix pheatmaps for all three variants.

#### A. F-matrix pheatmaps (protein × factor)

- Rows = 481 proteins, hierarchically clustered (ward.D2), no row labels
- Columns = factors (S2: PVE order); column annotation = |rho| with who_delta
- S1/S2: per-column max-normalisation → [0,1], gray96 → red scale
- S3: signed LFC, blue → white → red diverging scale
- Files: `heatmap_F_S1.pdf`, `heatmap_F_S2.pdf`, `heatmap_F_S3.pdf`

#### B. flashier single-cell vignette style (fit_S2)

| Plot | File |
|------|------|
| Scree/PVE bar chart | `s2_B1_scree.pdf` |
| Histogram × time_point | `s2_B2a_histogram_timepoint.pdf` |
| Histogram × severity | `s2_B2b_histogram_severity.pdf` |
| Structure plot × time_point | `s2_B3a_structure_timepoint.pdf` |
| Structure plot × severity | `s2_B3b_structure_severity.pdf` |
| L heatmap × time_point | `s2_B4a_heatmap_L_timepoint.pdf` |
| L heatmap × severity | `s2_B4b_heatmap_L_severity.pdf` |
| Top-protein F heatmap (F5,F1,F8) | `s2_B5_heatmap_F_top_proteins.pdf` |
| F scatter F5 (rho=-0.344) | `s2_B6_scatter_F05.pdf` |
| F scatter F1 (rho=-0.273) | `s2_B6_scatter_F01.pdf` |
| F scatter F8 (rho=+0.226) | `s2_B6_scatter_F08.pdf` |

Scatter plots: x = protein weight (F_pm), y = mean protein NPX across all 240 samples
(computed from stored data in fit), top 10 proteins by |F_pm| labelled.

**Note on structure plot:** `pm_colors` in flashier's structure plot sets FACTOR fill
colours (K colours required), not patient group colours. Patient groups (time_point or
severity) are shown via gaps and x-axis labels from `pm_groups`. Do not pass
`pm_colors` to structure plot calls when the goal is patient grouping.

**who_score_num values:** 1, 1.5, 2, 3, 4, 5, 6, 7 — note 1.5 = ambiguous patients
(n=85). Severity gradient: light → dark blue (8 levels).

#### C. flashier_intro vignette style — bar plots

| Plot | File |
|------|------|
| Bar plot of L × time_point | `s2_C1_bar_loadings_timepoint.pdf` |
| Bar plot of L × severity | `s2_C2_bar_loadings_severity.pdf` |

Each panel = one factor; each bar = one patient-timepoint sample; bars coloured by group.

---

### Saved objects

| File | Description |
|------|-------------|
| `results/nmf/stacked_S1_nmf.rds` | NMFfit object (NMF package) |
| `results/nmf/stacked_S2_ebnmf.rds` | flashier fit (point-exponential) |
| `results/nmf/stacked_S3_gbcd.rds` | gbcd fit (L, F$lfc, F$lfsr) |
| `results/nmf/stacked_association.csv` | Spearman rho + BH FDR, all variants × factors |
| `results/figs/nmf/stacked_S1_factor_proteins.pdf` | Top proteins per factor — S1 |
| `results/figs/nmf/stacked_S1_trajectories.pdf` | L_T1 vs L_T2 scatter, coloured by who_delta |
| `results/figs/nmf/stacked_S2_pve.pdf` | PVE bar chart — S2 |
| `results/figs/nmf/stacked_S2_factor_proteins.pdf` | Top proteins per factor — S2 |
| `results/figs/nmf/stacked_S2_trajectories.pdf` | L_T1 vs L_T2 scatter, coloured by who_delta |
| `results/figs/nmf/stacked_S3_factor_proteins.pdf` | LFC bar chart — S3 (shaded = lfsr≥0.05) |
| `results/figs/nmf/stacked_S3_trajectories.pdf` | L_T1 vs L_T2 scatter, coloured by who_delta |
| `results/figs/nmf/stacked_association_summary.pdf` | Rho waterfall, all variants |

### Installation note

`gbcd` install (first run only) updated `flashier` 1.0.7 → 1.0.58 and installed
`fastTopics` (GitHub dev). The install succeeded but caused a namespace conflict within the
same R session (old flashier 1.0.7 already loaded). Solution: run the script twice — first
run installs packages and fails at `library(gbcd)`; second run loads the updated flashier
cleanly. S1 and S2 cache to disk on the first run so only S3 actually re-runs. Subsequent
runs load all three from cache instantly.

## Results: 08_divas.R (proteomics T1 vs T2, 2-block DIVAS)

**Run date:** 2026-03-01. Completed successfully.

**Input:** 2 blocks, each 481 proteins × 120 patients (Olink NPX log2, column-centred — proteins centred within each block).
Block 1 = T1 (admission), Block 2 = T2 (follow-up). Same 120 patients in both blocks.

**Parameters:** `nsim=400`, `colCent=TRUE`, `rowCent=FALSE`, `seed=123`, `ReturnDetail=TRUE`.

### Decomposition

| Category | Count | Description |
|----------|-------|-------------|
| Shared (T1+T2) | 10 | Proteomics variation present at both time points |
| T1-individual | 8 | Variation unique to admission |
| T2-individual | 6 | Variation unique to follow-up |
| **Total** | **24** | |

### Severity associations

**T1 severity (Spearman rho):**

| Component | Category | rho | FDR |
|-----------|----------|-----|-----|
| **T1+T2-1** | Shared | **−0.721** | **<0.0001** |
| **T1+T2-2** | Shared | **−0.500** | **<0.0001** |
| T2-Individual-4 | T2-indiv | 0.229 | 0.094 |

**WHO delta (T2−T1 change, Spearman rho):**

| Component | Category | rho | FDR |
|-----------|----------|-----|-----|
| **T1+T2-6** | Shared | **−0.334** | **0.005** |
| T2-Individual-1 | T2-indiv | −0.202 | 0.318 |
| T1-Individual-5 | T1-indiv | −0.184 | 0.318 |

**Key findings:**
- **T1+T2-1** is the strongest severity-associated component (rho=−0.721, highly significant).
  This is a **time-stable severity program** — the same proteomics axis separates mild
  from severe at both admission and follow-up. Negative rho: high scores = lower severity.
- **T1+T2-6** is the strongest delta-associated component (rho=−0.334, FDR=0.005).
  Patients with higher scores on this shared component **improved** more (negative delta).
- Most individual components (T1-specific, T2-specific) have weak severity associations
  (all FDR > 0.3), confirming severity-related variation is predominantly **time-stable**.

### DIVAS naming conventions

- Scores matrix: `divas_res$Scores` (120 patients × 24 components)
- Score column names: `T1+T2-k` (shared), `T1-Individual-k`, `T2-Individual-k`
- For `getTopFeatures()`: shared use `compName="T1+T2-k"`, `modName="T1"` or `"T2"`;
  individual use `compName="T1-k"` (NOT "T1-Individual-k"), `modName="T1"`

### Outputs

| File | Description |
|------|-------------|
| `results/divas_2block.rds` | DIVAS result object |
| `results/divas_associations.csv` | Spearman rho with T1 severity and WHO delta, all 24 components |
| `results/figs/divas/divas_diagnostic.pdf` | DIVAS angle diagnostic plot |
| `results/figs/divas/divas_score_heatmap.pdf` | Score heatmap (components × patients) |
| `results/figs/divas/divas_severity_association.pdf` | Association waterfall |
| `results/figs/divas/divas_top_scatter.pdf` | T1+T2-1 scatter vs severity and delta |
| `results/figs/divas/divas_all_components_vs_severity.pdf` | All 24 components vs T1 severity |

---

## Key data facts (for NMF scripts)

- `who_delta` non-NA: 69 of 120 patients (51 excluded due to ambiguous WHO scores)
- `who_delta` range: -4 to +2; 42 stable (0), 16 improved (-1), 4 worsened (+1), 1 at +2
- X_shifted range: [0, 12.18] (per-protein min-shift of 240×481 NPX matrix)
- delta_mat range: [-6.19, 7.28] (T2 - T1 per patient, 120×481, signed)
- delta_shifted range: [0, 11.14] (for D1 NMF only)
