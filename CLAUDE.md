# CLAUDE.md — EBNMF Folder

> Session notes for the `EBNMF/` directory.
> This folder serves as the starting point for a mini-workshop on EBNMF, topic models, and DIVAS.

---

## Purpose of this folder

This folder contains source material and the generated document `EBNMF.md`, which will provide
a self-contained technical introduction for a **1-hour mini-workshop** covering:

1. The Empirical Bayes Normal Means (EBNM) problem and the `ebnm` R package
2. Empirical Bayes Matrix Factorization (EBMF) and the `flashier` R package
3. Topic models and their equivalence to Poisson NMF (`fastTopics`)
4. Interpreting programs via GoM DE (`fastTopics` + `ashr`)
5. Multi-block extensions: DIVAS for multiomics
6. The Generalized Binary Covariance Decomposition (GBCD) method

**Target audience:** Mixed — bioinformaticians with limited statistics background and
quantitative researchers. The document is a reference handout; the workshop uses concrete
biological examples throughout.

---

## Practical API guides (confirmed working, 2026-03)

| File | Covers |
|------|--------|
| `EBNMF_guide.md` | flashier (all prior families, result accessors, plot.flash gotchas) + GBCD (`fit_gbcd`, `F$lfc`, `F$lfsr`) |
| `DIVAS_guide.md` | DIVAS (`DIVASmain`, `getTopFeatures`, component naming gotchas, longitudinal + multiomics prep) |

These guides consolidate confirmed API patterns from the CovidCase analysis scripts
(`02_NMF_stacked.R`, `03_NMF_delta.R`, `08_divas.R`). Load them at the start of any
new session that uses flashier, GBCD, or DIVAS.

---

## Reference papers in this folder

| File | Reference | Key content |
|------|-----------|-------------|
| `ebnm.pdf` | Willwerscheid, Carbonetto & Stephens (*JSS* 2025) | EBNM problem, prior families, `ebnm` package API |
| `GBCD.pdf` | Liu, Carbonetto et al. (*Nature Genetics* 2025) | GBCD method, two-stage decomposition, HNSCC/PDAC applications |
| `GBCD_Supp.pdf` | Supplementary to GBCD paper | Nature Portfolio reporting summary only (3 pages) |
| `DIVAS_software.pdf` | Sun, Marron, Le Cao, Mao (2026, bioRxiv) | DIVAS R package, COVID-19 multiomics application |
| `CarbonettoTopicModel.pdf` | Carbonetto, Sarkar, Wang, Stephens (*arXiv* 2021) | Poisson NMF = topic model equivalence; `fastTopics` algorithms |
| `GoM_DE.pdf` | Carbonetto, Luo, Sarkar et al. (*Genome Biology* 2023) | GoM DE: interpreting NMF programs via graded membership DE |
| `DIVAS.pdf` | DIVAS theory paper | Theoretical background for DIVAS |

---

## Key dataset: COVID-19 multiomics cohort (DIVAS paper)

This is the **primary running example** for the workshop. It illustrates how the same data
can be decomposed using NMF in different ways depending on how the matrix is assembled.

**Dataset structure:**
- 139 COVID-19 patients, **two time points** (T1 and T2), severity score 1–7
- scRNA-seq pseudo-bulked for **4 immune cell types**: CD4 T, CD8 T, CD14 monocytes, NK cells
- Bulk proteomics + bulk metabolomics
- Total: **6 data blocks** (4 pseudobulk GEX + proteomics + metabolomics)
- For analysis: 120 samples (60 patients at T1, 60 different patients at T2, to satisfy independence)
- Source: Su et al. (*Cell* 2020); scRNA-seq accession E-MTAB-9357

**Why this dataset is pedagogically ideal:**
- Multiple blocks → demonstrates DIVAS multi-block decomposition
- Two time points + multiple cell types → demonstrates different matrix stacking strategies for NMF
- Clear biological gradient (severity) → programs can be validated against a known outcome
- Partially shared structure (e.g., CD8 T not in top component) → shows value of partial sharing

### Matrix stacking strategies for NMF on this data

The COVID-19 data illustrates three distinct ways to assemble X before running NMF:

| Strategy | Matrix X dimensions | What rows are | What columns are | What programs mean |
|----------|-------------------|---------------|------------------|--------------------|
| **Sample stacking** | (N·T) × p | Patients stacked across time points | Genes (one cell type) | Programs that vary across patients and/or time |
| **Feature stacking** | N × (p·B) | Patients (one time point) | Genes concatenated across blocks | Programs that are coordinated across modalities |
| **Block-wise NMF** | N × p_b separately per block | Patients | Genes per block | Block-specific programs, then compare across blocks |

Sample stacking example: stack T1 and T2 observations for CD14 monocytes → matrix (120 × p_CD14);
a factor that loads similarly at T1 and T2 is a stable program; one that differs captures
time-dependent change. Feature stacking: concatenate CD4 T and CD8 T gene columns → matrix
(N × (p_CD4 + p_CD8)); a factor spanning both captures coordinated T cell programs.
DIVAS generalises beyond NMF by identifying which blocks share which components without
fixing the stacking structure in advance.

---

## Relationship between the methods (Stephens group lineage)

```
EBNM (ebnm package)
  └─ EBMF (flashier package)            ← replaces Poisson NMF with EB priors; auto-selects K
       ├─ point-exponential prior → EBNMF (NMF with EB)
       ├─ point-normal prior → sparse MF
       ├─ generalized binary prior → GBCD
       └─ normal prior → PCA-like

Poisson NMF (fastTopics package)        ← equivalent likelihood to topic model (Carbonetto 2021)
  └─ GoM DE (fastTopics + ashr)         ← post-hoc: interprets programs via graded-membership DE

DIVAS                                   ← multi-block extension; finds jointly/partially shared
                                           components across K data blocks; uses subspace geometry
                                           (not EB); outputs components labelled by which blocks share them
```

**Key equivalences:**
- Poisson NMF ≡ multinomial topic model (same likelihood; Carbonetto et al. 2021, Lemma 1)
- GoM DE uses `ashr` (= EBNM with constrained nonparametric prior) for LFC shrinkage and lfsr
- `fastTopics` → better topic fits than EM; `flashier` → EB priors replace regularisation
- flashier's L_pm ↔ topic model's L (cell memberships); flashier's F_pm ↔ topic model's F (gene profiles)
  (key difference: flashier L is not row-normalised; topic model L sums to 1 per cell)

---

## Document improvement plan for EBNMF.md

The following changes are planned for the next editing session. Status: **pending**.

### User-requested improvements

**1. Expand EBNM → L/F bridge (currently one paragraph)**
Add a concrete derivation showing that fixing F and estimating row l_i of L is an EBNM
problem: the "observations" are the residual inner products R_i·f_k / ||f_k||² and the
standard errors are σ / ||f_k||. This is the conceptual engine of EBMF.

**2. Add an opening section: "EB in bioinformatics — you already know this"**
Before any formal math, connect EBNM to tools bioinformaticians already use:
- `limma`: EB shrinkage of gene-wise variance estimates toward a global mean (moderated t-test)
- `edgeR`/`DESeq2`: EB shrinkage of dispersion parameters (trend-based prior)
- `apeglm`/`ashr`: EB shrinkage of log-fold changes — direct precursor to `ebnm`
- One-sentence punchline: all of these borrow information across genes; EBNM formalises this idea.

**3. Add prior choice decision table (new §2.6)**
Replace scattered guidance with an actionable table:
- Decision axes: (i) sign constraint (positive-only vs signed), (ii) expected sparsity,
  (iii) tail behaviour
- Practical rule: fit two candidate priors, compare `logLik()`, pick the higher
- Bioinformatics defaults: NMF loadings → point-exponential; LFC → point-Laplace;
  PCA-like → normal; DE testing → point-normal or ash

**4. Fix matrix dimension convention throughout**
Standardise: X is n×p where n = cells/spots/samples (rows), p = genes/features (columns).
L is n×K, F is p×K. Fix all inconsistencies:
- EBMF section: already uses n×p but add explicit biological labels
- GBCD section: uses N×J (correct) but different symbols — unify to n×p
- Code examples: explain why `Matrix::t(lc_sparse)` is needed (SPE stores genes×spots)

### Additional improvements

**5. Add motivating toy example before the equations**
One paragraph: "Imagine 500 genes, 490 truly null. MLE uses all raw values. Point-normal
EBNM shrinks the 490 to near-zero, leaves the 10 signal genes near their observed values.
This is the same gain in power as increasing sample size — but from pooling information
across genes."

**6. Add new Part: Topic models and Poisson NMF**
Based on CarbonettoTopicModel.pdf and GoM_DE.pdf:
- The multinomial topic model (L row-normalised, F column-normalised, X = LF^T)
- Equivalence to Poisson NMF (Lemma 1: same likelihood up to row totals)
- EM = multiplicative updates; `fastTopics` CD with extrapolation outperforms EM
- Structure plot: the canonical visualisation of L (stacked bar, rows = cells)
- GoM DE: after fitting L, run per-gene Poisson regression with L as design; shrink LFCs
  with `ashr`; report lfsr; identify distinctive genes by least-extreme LFC
- Contrast with flashier: topic model fixes K and imposes row-sum-to-1 on L; flashier
  selects K automatically via ELBO and does not normalise L

**7. Add COVID-19 data stacking section (new Part before DIVAS)**
Use the COVID-19 multiomics cohort to illustrate:
- What "stacking along the sample dimension" means (T1 and T2 rows stacked; programs
  that are stable vs time-varying)
- What "stacking along the feature dimension" means (CD4 and CD8 columns concatenated;
  programs coordinated across cell types)
- What running NMF/topic model separately per block then comparing means
- How DIVAS generalises all three: it finds which blocks share components without
  committing to a stacking strategy

**8. Add "reading your NMF output" section (new §4.7)**
- Rank genes in factor k: `order(fit$F_pm[,k], decreasing=TRUE)[1:20]`
- Interpret fit$pve: steep drop-off = clean structure; flat tail = noisy factors
- Spot/cell loading maps: use fit$L_pm columns for spatial or UMAP visualisation
- When to merge or split factors

**9. Add clarification on log1p units and F interpretation**
When X = log1p(counts), entries of F are approximately log-fold changes in log1p space
relative to the zero-membership baseline. They are not in natural count units. Note this
when discussing gene weight interpretation.

**10. Connect lfsr to local FDR (one sentence)**
"lfsr is the Bayesian analogue of the local FDR (Efron 2004): it is a posterior probability
— no Bonferroni or BH correction is needed."

**11. Add method-choice summary table (end of document)**
One-page table: standard NMF, EBNMF (flashier), topic model (fastTopics), GoM DE,
GBCD, DIVAS — keyed by question, data type, K selection, sparsity mechanism, output.

### Revised document structure (target)

```
[New] Preamble  — EB in bioinformatics: limma, edgeR, ashr → EBNM connection
Part I          — EBNM problem (current §1.1–1.3, + toy intuition example)
Part II         — Prior families (current §2.1–2.5, + decision table §2.6)
Part III        — ebnm package (current §3.1–3.4, + lfsr↔local FDR note)
[Expanded]
Part IV         — EBMF and flashier (current, + expanded §4.2 bridge derivation,
                  + output interpretation §4.7, + log1p units note)
[New]
Part V          — Topic models and Poisson NMF (Carbonetto 2021 + GoM DE)
                  §5.1 The multinomial topic model
                  §5.2 Equivalence to Poisson NMF
                  §5.3 Algorithms: EM vs fastTopics CD
                  §5.4 Structure plot visualisation
                  §5.5 GoM DE: graded-membership differential expression
                  §5.6 flashier vs topic model: comparison
[New]
Part VI         — Multi-block data and stacking strategies (COVID-19 example)
                  §6.1 Sample-dimension stacking (time points stacked)
                  §6.2 Feature-dimension stacking (modalities concatenated)
                  §6.3 Block-wise NMF then comparison
                  §6.4 DIVAS: principled multi-block decomposition
Part VII        — GBCD (current §5 renumbered, convention fixes)
Part VIII       — Method comparison (current §6, + summary choice table)
Part IX         — Biological applications (current §7, now includes COVID-19 DIVAS result)
Part X          — Pipeline summary (current §8, convention fixes)
```

---

## Key equations from new reference papers

### CarbonettoTopicModel.pdf — Poisson NMF ≡ topic model

| Equation | Description |
|----------|-------------|
| $x_{ij} \mid \mathbf{H},\mathbf{W} \sim \text{Pois}((\mathbf{HW}^\top)_{ij})$ | Poisson NMF likelihood |
| $x_{i\cdot} \sim \text{Multin}(s_i;\, (\mathbf{LF}^\top)_{i\cdot})$ | Multinomial topic model |
| $p_\text{PNMF}(\mathbf{X}\mid\mathbf{H},\mathbf{W}) = p_\text{MTM}(\mathbf{X}\mid\mathbf{L},\mathbf{F})\prod_i\text{Pois}(t_i;s_i)$ | Lemma 1: same likelihood up to row totals |
| $\mathbf{H}^\text{ext} \leftarrow P_+[\mathbf{H}^\text{new} + \beta^{(t)}(\mathbf{H}^\text{new}-\mathbf{H}^{(t-1)})]$ | Extrapolated CD update (faster than EM) |

### GoM_DE.pdf — GoM differential expression

| Equation | Description |
|----------|-------------|
| $x_{ij} \sim \text{Pois}(s_i \sum_k l_{ik} p_{jk})$ | GoM DE model (L fixed from topic fit) |
| $\text{LFC}_{k,l}(j) = \log_2(p_{jk}/p_{jl})$ | Pairwise LFC, topic k vs l |
| $\text{LFC}_k^\text{l.e.}(j) = \text{LFC}_{k,l^*}(j),\; l^* = \arg\min_{l\ne k}\lvert\text{LFC}_{k,l}(j)\rvert$ | Least-extreme LFC: distinctive only if DE vs all others |
| $\text{lfsr}(k,j) = \min\{P(\text{LFC}_k^\text{l.e.}(j) \le 0\mid\cdot), P(\text{LFC}_k^\text{l.e.}(j) \ge 0\mid\cdot)\}$ | lfsr from MCMC posterior |

### DIVAS_software.pdf — multi-block decomposition

| Equation | Description |
|----------|-------------|
| $X_k = A_k + E_k$ | Each block = signal + noise |
| $A_k = \sum_{\mathbf{i}\ni k} L_{\mathbf{i},k} S_\mathbf{i}^\top$ | Block decomposition: loading × shared score |
| $S_\mathbf{i} \in \mathbb{R}^{r_\mathbf{i}\times N}$, $L_{\mathbf{i},k} \in \mathbb{R}^{d_k\times r_\mathbf{i}}$ | Scores shared across blocks in subset $\mathbf{i}$; loadings block-specific |
| Rotational bootstrap | Inference for which components are statistically significant |

---

## Mini-workshop structure (revised, 1 hour)

**Format:** 45 min exposition + 15 min hands-on. EBNMF.md is the handout.

### Block 1 (10 min) — EB is everywhere in bioinformatics
- limma, edgeR, ashr as familiar examples
- EBNM as the unifying framework
- Key idea: borrow information across features; the prior is estimated, not fixed

### Block 2 (15 min) — From scalar shrinkage to matrix factorization
- EBNM problem and point-exponential prior (intuition + one equation)
- The EBMF model: X = LF^T + E
- How fixing F turns the L-estimation into repeated EBNM (the bridge)
- Greedy + backfit + nullcheck: auto K selection
- Prior family → decomposition type table (normal/PCA, point-normal/sparse, point-exp/NMF, GB/GBCD)

### Block 3 (10 min) — Topic models: the bioinformatics cousin
- Topic model = Poisson NMF (same likelihood, Lemma 1)
- Structure plot: reading membership matrices
- GoM DE: how to annotate programs with genes (least-extreme LFC + lfsr)
- fastTopics vs flashier: when to use which

### Block 4 (10 min) — Multi-block data: stacking strategies and DIVAS
- COVID-19 example: 6 blocks, 2 time points
- Three stacking strategies and what each reveals
- DIVAS: finds shared/partially shared/individual components with bootstrap inference
- Key result: partially shared CD4-CD14-NK-pro-metab component most correlated with severity

### Block 5 (15 min) — Hands-on (R)
Participants run on the COVID-19 or a toy dataset:
1. `ebnm_point_exponential()` on a vector → inspect shrinkage and lfsr
2. `flash()` on one data block → inspect fit$pve and structure plot of fit$L_pm
3. (Optional) `fastTopics::fit_topic_model()` → compare to flash output
4. Discuss: what would change if you stacked T1 and T2?

---

## Notes for future sessions

- DIVAS_software.pdf read in full. COVID-19 dataset details extracted (see above).
- CarbonettoTopicModel.pdf and GoM_DE.pdf read in full. Key equations extracted above.
- The DIVAS.pdf (theory paper) has not yet been read; read it if the workshop covers
  the algorithmic details of the subspace angle decomposition.
- ebnm.pdf pages 21–32 (appendices with benchmarking) not in EBNMF.md; consider adding
  for workshop practical session on optimiser speed.
- GBCD_Supp.pdf: Nature Portfolio summary only (3 pages); full supplementary figures
  not downloaded.
- EBNMF.md edits: not yet made. The improvement plan above is the task for the next session.
  Implement sections in order: Preamble → §2.6 prior table → §4.2 bridge → new Parts V and VI.
