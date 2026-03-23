---
output:
  html_document:
    toc: true
    toc_depth: 3
  pdf_document:
    toc: true
    toc_depth: 3
---

# DIVAS — Multi-block Data Decomposition: Companion Document

This document is a companion to **EBNMF.md** (the mini-workshop reference handout). It covers
multi-block data analysis strategies and the DIVAS method, illustrated with the COVID-19
multiomics cohort from Su et al. (*Cell* 2020) re-analysed in Sun, Marron, Le Cao & Mao (*bioRxiv* 2026).

**Key reference:** Sun Y, Marron JS, Le Cao KA, Mao J (2026). "DIVAS: an R package for identifying
shared and individual variations of multiomics data." *bioRxiv*.

---

## The COVID-19 multiomics cohort

This dataset (Su et al., *Cell* 2020; re-analysed in the DIVAS paper, Sun et al. 2026)
is a concrete example for discussing how to decompose data with multiple blocks and
repeated measures.

**Structure:**

- 139 COVID-19 patients, severity score 1–7 (ordinal)
- **Two time points:** T1 (early, ~days 1–4 post admission) and T2 (later, ~days 7–10)
- scRNA-seq pseudo-bulked for **4 immune cell types**: CD4 T, CD8 T, CD14 monocytes, NK cells
- Bulk proteomics (plasma proteins) and bulk metabolomics
- Total: **6 data blocks** sharing the same set of patients

For methods that require independent samples, 60 patients contributed T1 observations
and 60 different patients contributed T2 observations (120 total), breaking the
repeated-measures structure.

**Key biological finding (DIVAS analysis):** The component most correlated with COVID-19
severity was a *partially shared* component spanning CD4 T, CD14 monocytes, NK cells,
proteomics, and metabolomics — but **not** CD8 T cells. This would have been missed by
methods that only seek jointly shared (all-block) variation. The component captured
activated but metabolically stressed CD4 T cells, cytokine storm markers (CCL7, TNFRSF10B),
and negative correlation between IFNGR2 (myeloid activation) and PDP1 (oxidative metabolism).

---

## Strategy 1 — Sample stacking (time points as additional rows)

Concatenate T1 and T2 samples into a single matrix:

$$\mathbf{X} = \begin{pmatrix} \mathbf{X}_{T1} \\ \mathbf{X}_{T2} \end{pmatrix} \in \mathbb{R}^{(n_1 + n_2) \times p}$$

where rows are samples ($n_1$ at T1, $n_2$ at T2) and columns are genes. NMF or EBMF then finds
programs in this combined matrix.

**What the loadings reveal:**

- A program with similar $L_{ik}$ at T1 and T2 is **temporally stable** (a persistent biological state)
- A program with higher $L_{ik}$ at T2 than T1 is **time-increasing** (disease progression, immune activation)
- A program with loadings specific to T1 or T2 captures **time-point-specific** biology

**Limitation:** if T1 and T2 samples differ in library size or batch effects, these technical
differences will show up as the dominant programs and must be addressed (e.g., normalise
jointly before stacking).

---

## Strategy 2 — Feature stacking (blocks as additional columns)

Concatenate gene columns from multiple blocks (e.g., two cell types) into a single matrix:

$$\mathbf{X} = \begin{pmatrix} \mathbf{X}_{\text{CD4}} \;\big|\; \mathbf{X}_{\text{CD8}} \end{pmatrix} \in \mathbb{R}^{n \times (p_{\text{CD4}} + p_{\text{CD8}})}$$

where rows are samples and columns are genes concatenated across cell types.

**What the factors reveal:**

- A factor whose $\mathbf{F}$ sub-vector is non-zero in both the CD4 and CD8 partitions
  captures a **coordinated program shared across cell types**
- A factor whose $\mathbf{F}$ sub-vector is concentrated in only one partition
  represents a **cell-type-specific program**
- $\mathbf{F}$ naturally decomposes into sub-vectors: $\mathbf{F} = [\mathbf{F}_{\text{CD4}}^\top \;|\; \mathbf{F}_{\text{CD8}}^\top]^\top$

**Limitation:** features from different blocks may be on different scales and have different
sparsity. Normalise within each block before concatenation, and consider down-weighting
blocks with many more features.

---

## Strategy 3 — Block-wise NMF then cross-block comparison

Run NMF independently on each block:

$$\mathbf{X}_b = \mathbf{L}_b \mathbf{F}_b^\top + \mathbf{E}_b, \quad b = 1, \ldots, B.$$

Then compare the loading matrices $\mathbf{L}_1, \ldots, \mathbf{L}_B$ across blocks.
Programs that are shared across blocks will have correlated loading vectors (i.e.,
$\text{cor}(\ell_k^{(b)}, \ell_k^{(b')}) \approx 1$ for the same biological program
discovered in blocks $b$ and $b'$). Programs unique to one block will have no counterpart.

**Practical tool for comparison:** cosine similarity between loading vectors across blocks,
or Hungarian algorithm matching of factors by correlation.

**Advantage:** block-wise NMF avoids scale and dimensionality mismatch between blocks.
**Limitation:** no principled test for whether similarity between blocks is statistically
significant, and programs may not align across blocks even when they represent the same biology
(rotation non-uniqueness of NMF).

---

## DIVAS: principled multi-block decomposition

DIVAS (Dimension-reduction and Integration of multi-block data via Angle-based Subspace
analysis; Sun, Marron, Le Cao, Mao 2026) generalises beyond the three stacking strategies
by explicitly modelling which subset of blocks shares each component.

**Model.** Each block decomposes as:

$$\mathbf{X}_k = \mathbf{A}_k + \mathbf{E}_k, \qquad \mathbf{A}_k = \sum_{\mathbf{i} \ni k} \mathbf{L}_{\mathbf{i},k}\, \mathbf{S}_{\mathbf{i}}^\top,$$

where:

- $\mathbf{i} \subseteq \{1,\ldots,B\}$ is a *subset* of block indices
- $\mathbf{S}_{\mathbf{i}} \in \mathbb{R}^{r_{\mathbf{i}} \times N}$ is a **score matrix shared across all blocks in subset $\mathbf{i}$**
- $\mathbf{L}_{\mathbf{i},k} \in \mathbb{R}^{d_k \times r_{\mathbf{i}}}$ is a **block-specific loading matrix**, translating shared scores to block-$k$-specific feature space

The sum runs over all subsets $\mathbf{i}$ that contain block $k$: a component that appears
in all 6 blocks contributes to every block; a component shared by only CD4 and CD14 contributes
only to those two blocks.

**Note on conventions:** DIVAS uses features-as-rows ($\mathbf{X}_k \in \mathbb{R}^{d_k \times N}$,
opposite to the EBNMF.md convention). The shared scores $\mathbf{S}_{\mathbf{i}}$ are the
quantities directly comparable across blocks; the block-specific loadings $\mathbf{L}_{\mathbf{i},k}$
map them to feature space.

**Algorithm (two steps):**

1. **Denoising.** Apply PCA to each $\mathbf{X}_k$ separately. The number of PCs is
   selected by optimal singular value shrinkage (Gavish & Donoho 2017): discard singular
   values below the median Marchenko-Pastur threshold, which is the optimal hard threshold
   under a white-noise assumption.

2. **Subspace identification.** Use angle-based subspace analysis (Prothero et al. 2024):
   find shared subspace directions that are close (small principal angles) to all included
   blocks and far (large principal angles) from excluded blocks. This proceeds hierarchically:
   jointly shared components (all $B$ blocks) are found first, then components shared by
   $B-1$ blocks, then $B-2$, down to individual block-specific components.

**Inference.** Statistical significance of each component is assessed by a **rotational
bootstrap**: resample and re-run the subspace decomposition; components that are stable
across resamples are retained.

**Comparison with NMF stacking:**

+-------------------+-----------------+-----------------+-----------------+------------------------------+
|                   | Sample stacking | Feature stacking| Block-wise NMF  | DIVAS                        |
+===================+=================+=================+=================+==============================+
| K selection       | Manual          | Manual          | Manual per block| Auto (SVD + bootstrap)       |
+-------------------+-----------------+-----------------+-----------------+------------------------------+
| Partial sharing   | Implicit only   | Yes, from F     | Yes, from L     | Explicit, tested statistically|
+-------------------+-----------------+-----------------+-----------------+------------------------------+
| Scale/dim. mismatch| Moderate       | Severe          | None (separate) | None (PCA per block)         |
+-------------------+-----------------+-----------------+-----------------+------------------------------+
| Inference         | None            | None            | None            | Rotational bootstrap         |
+-------------------+-----------------+-----------------+-----------------+------------------------------+
| Interpretability  | Time-varying    | Cross-block     | Per-block       | Labelled block subsets       |
|                   | programs        | gene coord.     | programs        |                              |
+-------------------+-----------------+-----------------+-----------------+------------------------------+

---

## COVID-19 DIVAS application result

DIVAS applied to 6 multiomics blocks (4 pseudobulk GEX + proteomics + metabolomics) from
120 COVID-19 samples identified over 90 statistically significant components.

**Top component** (most correlated with severity): **CD4-CD14-NK-proteomics-metabolomics**
(partially shared; **CD8 T cells excluded**). Biological signal:

- CD4 T: activated but metabolically stressed (mTORC1 suppression, glycolytic shift)
- CD14 monocytes: high IFNGR2 (cytokine receptor), low PDP1 (oxidative metabolism enzyme)
- NK cells and proteomics: cytokine storm markers (CCL7, TNFRSF10B, quinolinate)

This partial sharing — CD8 T cells uninvolved — is a biologically specific finding that
would not emerge from methods restricted to jointly shared (all-block) variation. Sample
stacking would merge T1 and T2 variation; feature stacking would require committing to a
pre-specified block concatenation order; block-wise NMF provides no statistical test for
whether the observed correlation between blocks is real. DIVAS identifies and confirms
this partial structure via the rotational bootstrap.

**COVID-19 dataset source:** Su et al. (*Cell* 2020); scRNA-seq accession E-MTAB-9357.
