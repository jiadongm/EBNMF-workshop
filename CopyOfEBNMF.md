---
title: "Mini-workshop: Matrix factorisation for longitudinal omics data with EBNMF and DIVAS"
author:
  - name: Jiadong Mao
    email: chiatungmao@gmail.com
output:
  html_document:
    mathjax: "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"
    toc: true
    toc_depth: 3
  pdf_document:
    toc: true
    toc_depth: 3
---


---

### Statement about the use of AI 

Claude Opus 4.6 and Sonnet 4.6 were used when generating this document. A draft was first generated
by Claude Code using an outline and references provided by me. I then developed a case study using 
the COVID-19 data using Claude. The output was then incorporated into the draft. All writting and code 
were checked and revised.


---


## Preamble

[To Claude: Summarise this document, i.e. it's a workshop handout for introducing EBNMF and 
related tools and DIVAS. The following references section should be divided into two
parts, one for EBNMF and one for DIVAS.]

**Primary references**

- Willwerscheid J, Carbonetto P, Stephens M (2025). "ebnm: An R Package for Solving the Empirical
  Bayes Normal Means Problem Using a Variety of Prior Families." *Journal of Statistical
  Software*, **114**(3). doi:10.18637/jss.v114.i03
- Liu Y, Carbonetto P, Willwerscheid J, Oakes SA, Macleod KF, Stephens M (2025). "Dissecting
  tumor transcriptional heterogeneity from single-cell RNA-seq data by generalized binary
  covariance decomposition." *Nature Genetics*, **57**, 263–273. doi:10.1038/s41588-024-01997-z
- Carbonetto P, Sarkar A, Wang Z, Stephens M (2022). "Non-negative Matrix Factorization
  Algorithms Generally Improve Topic Model Fits." *arXiv* 2105.13440.
- Carbonetto P, Luo K, Sarkar A, Hung A, Tayeb K, Pott S, Stephens M (2023). "GoM DE:
  interpreting structure in sequence count data with differential expression analysis allowing
  for grades of membership." *Genome Biology*, **24**, 236.

[To Claude: Add a code chunk here guiding the user to install all used packages in this session.]

---

## Part I — The Empirical Bayes Normal Means Problem


### 1.0 You might have already used empirical Bayes bioinformatic tools

If you have used any of the following tools, 
you have already applied empirical Bayes (EB) to omics data.

**`limma` (Smyth 2004).** When testing differential expression, `limma` does not use the
gene-wise variance estimate $s_j^2$ directly (these are noisy with few replicates). Instead
it borrows information across all genes to estimate a global variance distribution, then
*shrinks* each $s_j^2$ toward that estimate. The result is a "moderated" t-statistic that
behaves like a t-test with more degrees of freedom than the experiment actually has. The
borrowing is done by fitting a scaled inverse-chi-squared prior to the $s_j^2$: the prior
is estimated from the data (empirical Bayes), then used to form posterior estimates.

**`edgeR` and `DESeq2`.** For RNA-seq count data, each gene has an overdispersion parameter
$\phi_j$ that must be estimated. With typical sample sizes ($n = 3$–10), the per-gene MLE
$\hat{\phi}_j$ is extremely noisy. Both packages estimate a dispersion trend across genes
(the EB prior), then shrink each $\hat{\phi}_j$ toward that trend. The shrinkage avoids
inflated test statistics for genes with unusually low estimated dispersion.

**The unifying idea.** In all three cases the logic is the same:

1. You have many noisy estimates $\hat{\theta}_1, \ldots, \hat{\theta}_n$ (variances, dispersions, LFCs).
2. These estimates share a common underlying distribution (the prior $g$).
3. Estimating $g$ from the data and using it to form posterior estimates $\mathbb{E}[\theta_i \mid \hat{\theta}_i, \hat{g}]$ *always* reduces mean squared error compared to using $\hat{\theta}_i$ directly.
4. The prior is *not* specified in advance — it is estimated from the data. This is what makes it *empirical* Bayes.

The **Empirical Bayes Normal Means (EBNM)** problem formalises this idea for the case where
each $\hat{\theta}_i$ is normally distributed with known standard error $s_i$. The `ebnm` R
package solves this problem for a wide range of prior families. When applied to the rows and
columns of a matrix, EBNM becomes **Empirical Bayes Matrix Factorization (EBMF)**, implemented
in the `flashier` package. This is the thread we follow in this document.

### 1.1 The normal means model

**Technical explanation.** The starting point is the *normal means model*: given $n$ observations
with known standard deviations, we assume

$$x_i \overset{\text{ind.}}{\sim} \mathcal{N}(\theta_i,\, s_i^2), \quad i = 1, \ldots, n,$$

where the $\theta_i \in \mathbb{R}$ are the unknown true means we wish to estimate. Because
$s_i$ is known, the MLE is trivially $\hat{\theta}_i = x_i$. However, as Stein (1956) famously
showed, this MLE is *inadmissible* under squared error loss when $n \ge 3$: there exist
estimators with uniformly lower risk. The key insight is that each $x_i$ carries information
not only about its own $\theta_i$ but also about the collective distribution of all the
$\theta_i$, and exploiting this shared information improves estimation.

**Intuition.** Imagine a DE experiment with $n = 500$ genes. The truth is that 490 genes have
zero effect ($\theta_i = 0$) and 10 have a real effect ($\theta_i \ne 0$). The 490 null genes
produce $x_i$ values clustered near zero (noise only); the 10 signal genes scatter far from
zero. When we maximise the marginal likelihood

$$L(g) = \prod_{i=1}^n \int p(x_i \mid \theta_i)\, g(\theta_i)\, d\theta_i$$

over $g$, the algorithm *discovers from the histogram of $x_i$* that the data are best
explained by a $g$ placing a large spike at zero (accounting for the cluster of null genes)
plus a slab (accounting for the outlying signal genes). No one told the algorithm that most
genes are null — it inferred this structure from the data. The MLE uses all 500 raw
observations at face value, so the 490 null genes contaminate rankings. An EB estimator
that learns "most genes are null" shrinks the 490 null genes hard toward zero while leaving
the 10 signal genes near their observed values — improving power and ranking simultaneously,
at no extra cost.

```r
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
```

### 1.2 The empirical Bayes approach

The EBNM framework assumes the true means are drawn from a common prior in some family $\mathcal{G}$:

$$\theta_i \overset{\text{ind.}}{\sim} g \in \mathcal{G}.$$

Rather than fixing $g$ in advance (fully Bayes), EBNM *estimates* $g$ from the data by
maximising the marginal likelihood:

$$\hat{g} := \underset{g \in \mathcal{G}}{\arg\max}\ L(g), \qquad
L(g) := \prod_{i=1}^n \int p(x_i \mid \theta_i, s_i)\, g(\theta_i)\, d\theta_i.$$

Given $\hat{g}$, posterior summaries of each $\theta_i$ follow from Bayes' theorem:

$$p(\theta_i \mid x_i, s_i, \hat{g}) \propto p(x_i \mid \theta_i, s_i)\, \hat{g}(\theta_i).$$

The posterior mean $\hat{\theta}_i = \mathbb{E}(\theta_i \mid x_i, s_i, \hat{g})$ is the shrinkage
estimator. Borrowing of information occurs because $\hat{g}$ is estimated by pooling all $n$
observations: regions where the prior places little mass get shrunk most strongly.

### 1.3 Why shrinkage works

**Intuition.** Shrinkage trades a small amount of bias for a large reduction in variance —
this is always worthwhile when the signal is weak relative to the noise. Consider the two
extreme cases:

- **Null gene ($\theta_i = 0$):** The MLE is $x_i$, which equals pure noise. The EB posterior
  places the estimate near zero, almost eliminating the noise contribution. The bias is
  essentially zero (we moved toward the truth) and the variance is greatly reduced.
- **Signal gene (large $|\theta_i|$):** A large $|x_i|$ is only weakly pulled toward zero;
  the posterior mean stays close to $x_i$ because the prior slab places mass at large values.
  Shrinkage is minimal for strong signals.

The EB estimator automatically calibrates the degree of shrinkage to the signal strength —
heavy shrinkage where signal is absent, minimal shrinkage where signal is strong. This is
impossible for any fixed-regularisation method.

```r
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
```

**Technical explanation.** With a normal prior the posterior is Gaussian:

$$\theta_i \mid x_i, s_i, \hat{g} \sim \mathcal{N}\!\left(\frac{\hat{\sigma}^2}{\hat{\sigma}^2 + s_i^2} x_i + \frac{s_i^2}{\hat{\sigma}^2 + s_i^2} \hat{\mu},\;\frac{\hat{\sigma}^2 s_i^2}{\hat{\sigma}^2 + s_i^2}\right),$$

where $\hat{\mu}$ and $\hat{\sigma}^2$ are estimated from data. Observations with small $s_i$
stay close to $x_i$; those with large $s_i$ are pulled toward $\hat{\mu}$. For sparse priors
(point-normal, point-exponential) the shrinkage is nonlinear: observations near zero are
collapsed to exactly zero while large observations survive largely intact.

---

## Part II — Prior Families in `ebnm`

The choice of $\mathcal{G}$ governs the shrinkage profile and what the posterior can express.
The `ebnm` package provides a unified interface for a wide range of families.

### 2.1 Normal family

$$\mathcal{G}_{\text{norm}} = \{ g : g(x) = \mathcal{N}(x;\mu,\sigma^2),\; \sigma^2 \ge 0,\; \mu \in \mathbb{R} \}$$

The posterior is Gaussian with a closed form. For homoskedastic data ($s_i = s$):

$$\hat{\mu} = \frac{1}{n}\sum_{i} x_i, \qquad
\hat{\sigma}^2 = \max\!\left\{0,\; \frac{1}{n}\sum_{i} (x_i - \hat{\mu})^2 - s^2\right\}.$$

**Limitation:** The normal prior over-shrinks extreme observations (light tails) and cannot
produce exact zeros. It is appropriate when you expect all signals to be non-zero and
roughly Gaussian.

### 2.2 Sparse (spike-and-slab) priors

Many genomics applications call for priors that capture *sparsity* — most $\theta_i$ are
zero or near-zero, a few are large.

**Point-normal** (`"point_normal"`, support $\pm$):

$$\mathcal{G}_{\text{pn}} = \{ g : g(x) = \pi_0 \delta_0(x) + (1-\pi_0)\mathcal{N}(x;\mu,\sigma^2),\; 0 \le \pi_0 \le 1,\; \sigma^2 > 0 \}$$

The spike at zero captures exact zeros; the Gaussian slab captures non-zero signals. The
posterior is a mixture of a point mass and a Gaussian (soft-thresholding). This prior is
symmetric and suited to quantities that can be positive or negative (e.g., DE LFCs,
PC loadings).

**Point-Laplace** (`"point_laplace"`, support $\pm$):

$$\mathcal{G}_{\text{pl}} = \{ g : g(x) = \pi_0 \delta_0(x) + (1-\pi_0)\text{Laplace}(x;\mu,a),\; 0 \le \pi_0 \le 1,\; a \ge 0 \}$$

The Laplace slab has heavier tails than the Gaussian slab, leading to less over-shrinkage
of large signals. Johnstone & Silverman (2005) showed this generally improves accuracy.
Good default for signed LFCs and gene signatures.

**Point-exponential** (`"point_exponential"`, support $+$ only):

$$\mathcal{G}_{\text{pe}} = \{ g : g(x) = \pi_0 \delta_0(x) + (1-\pi_0)\text{Exp}(x;a),\; 0 \le \pi_0 \le 1,\; a \ge 0 \}$$

The exponential slab has support only on the positive reals, enforcing non-negativity on
the estimated $\theta_i$. The spike allows exact zeros. This is the prior used for NMF
in `flashier`: gene weights and sample loadings are pushed to be either zero (not in
program) or positive (contributes to program), with sparsity learned from the data.

**Horseshoe** (`"horseshoe"`, support $\pm$):

$$\mathcal{G}_{\text{hs}} = \{ g : g(x) = \text{Horseshoe}(x;\tau) \}$$

Places appreciable mass near zero without an exact point mass. Useful when true sparsity
is uncertain or when you expect a continuous spectrum from small to large signals.

### 2.3 Constrained nonparametric priors

**Normal scale mixture** (`"normal_scale_mixture"`, support $\pm$):

$$\mathcal{G}_{\text{smn}} = \left\{ g : g(x) = \int_0^\infty \mathcal{N}(x;\,0,\sigma^2)\, dh(\sigma^2) \right\}$$

All distributions in this family are symmetric and unimodal at zero. This is the family
used by `ashr` for LFC shrinkage in DE analysis — it is the direct precursor to `ebnm`.

**Unimodal symmetric** (`"unimodal_symmetric"`, support $\pm$) and
**unimodal nonnegative** (`"unimodal_nonnegative"`, support $+$):

$$\mathcal{G}_{\text{uni-nn}} = \left\{ g : g(x) = \int_0^\infty \text{Unif}(x;\,0,a)\, dh(a) \right\}$$

More flexible than normal scale mixtures: the shape is learned from data subject only to
unimodality. Useful when you are unsure about tail behaviour.

### 2.4 Nonparametric maximum likelihood estimate (NPMLE)

When $\mathcal{G}$ is the unconstrained family of all distributions, $\hat{g}$ is the
nonparametric MLE. In practice, $\mathcal{G}$ is approximated by a finite mixture of
point masses on a dense grid $\{\mu_1, \ldots, \mu_K\}$:

$$\tilde{\mathcal{G}}_{\text{npmle}} = \left\{ g : g(x) = \sum_{k=1}^K \pi_k \delta_{\mu_k}(x) \;\middle|\; \pi_k \ge 0,\; \sum_{k} \pi_k = 1 \right\}$$

Maximising the log-marginal likelihood over $\boldsymbol{\pi}$ is a convex problem solved
by `mixsqp`. **Key limitation:** the resulting posteriors are discrete, making posterior
interval estimates unreliable.

### 2.5 Nesting relationships

The families are nested in terms of flexibility:

$$\mathcal{G}_{\text{norm0}} \subset \mathcal{G}_{\text{pn}} \subset \mathcal{G}_{\text{smn}} \subset \mathcal{G}_{\text{symm-u}} \subset \mathcal{G}_{\text{npmle}}$$

More flexible families achieve higher marginal log-likelihood but may overfit: the RMSE
of NPMLE posteriors is typically *worse* than point-normal because the discrete posterior
overfits the observations.

### 2.6 How to choose a prior — decision guide

This is the most practically important question for new users. The table below summarises
prior behaviour along four axes; use it to narrow down candidates, then confirm with `logLik()`.

+----------------------+-----------+----------+-------------+----------------------------------+
| Prior                | Support   | Sparsity | Tail        | Typical use                      |
+======================+===========+==========+=============+==================================+
| Normal               | $\pm$     | None     | Gaussian    | PCA-like, dense factors          |
+----------------------+-----------+----------+-------------+----------------------------------+
| Point-normal         | $\pm$     | Moderate | Gaussian    | Sparse signed loadings (PCs)     |
+----------------------+-----------+----------+-------------+----------------------------------+
| Point-Laplace        | $\pm$     | Moderate | Heavy       | DE log-fold changes, gene sigs   |
+----------------------+-----------+----------+-------------+----------------------------------+
| Point-exponential    | $+$ only  | Strong   | Light       | NMF loadings, expression progs   |
+----------------------+-----------+----------+-------------+----------------------------------+
| Horseshoe            | $\pm$     | Soft     | Very heavy  | Uncertain sparsity, robust DE    |
+----------------------+-----------+----------+-------------+----------------------------------+
| Unimodal (sym/nn)    | $\pm/+$   | Flexible | Flexible    | Unknown tail, data-adaptive      |
+----------------------+-----------+----------+-------------+----------------------------------+
| NPMLE                | $\pm$     | Maximum  | Free        | Large $n$, maximum flexibility   |
+----------------------+-----------+----------+-------------+----------------------------------+

**Practical rule for comparing priors:**

```r
library(ebnm)
fit_pe  <- ebnm_point_exponential(x, s)
fit_pl  <- ebnm_point_laplace(x, s)
fit_pn  <- ebnm_point_normal(x, s)

logLik(fit_pe)   # pick the family with the highest value
logLik(fit_pl)
logLik(fit_pn)
```

The `logLik()` method returns $\log L(\hat{g})$, the marginal log-likelihood at the
estimated prior. Because all three call the same model ($x_i \sim \mathcal{N}(\theta_i, s_i^2)$,
$\theta_i \sim g$), the values are directly comparable — a higher log-likelihood means the
family fits the data better and will typically give better-calibrated posteriors.

**Bioinformatics defaults:**

- **NMF loadings and gene weights** → `ebnm_point_exponential` (must be non-negative)
- **DE log-fold changes** → `ebnm_point_laplace` (heavy tails, symmetric around zero)
- **PCA-like global structure** → `ebnm_normal` (dense, no zeros)
- **`ashr` use case** (LFC shrinkage after DE) → `ebnm_normal_scale_mixture` or `ebnm_unimodal_symmetric` (see below)

**Prior shape visualisation.** The following code plots the four most commonly used priors
side by side. Each panel shows the continuous slab component (curve) plus, where applicable,
a vertical segment representing the point mass at zero.

```r
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
```

---

## Part III — The `ebnm` Package: Interface and Outputs *(Optional)*

### 3.1 Main function signature

```r
ebnm(x,               # vector of observations x_1, ..., x_n
     s,               # vector (or scalar) of standard errors s_1, ..., s_n
     prior_family,    # string: "normal", "point_normal", "point_exponential", etc.
     mode = 0,        # prior mode (or "estimate" to estimate from data)
     scale = "estimate",  # prior scale parameter (or estimated)
     g_init = NULL,   # initial estimate of g (speeds up MLE search)
     fix_g = FALSE,   # if TRUE, compute posteriors with g fixed to g_init
     output = ebnm_output_all())
```

Alternative interface: each prior family has its own function, e.g.,
`ebnm_point_exponential(x, s, ...)`, `ebnm_normal(x, s, ...)`, etc.

### 3.2 Two-step optimisation procedure

**Step 1 — Estimate the prior.** Compute:

$$\hat{g} := \underset{g \in \mathcal{G}}{\arg\max}\ L(g), \qquad
L(g) := \prod_{i=1}^n \int p(x_i \mid \theta_i, s_i)\, g(\theta_i)\, d\theta_i.$$

For parametric priors ($\mathcal{G}_{\text{pn}}, \mathcal{G}_{\text{pe}}, \mathcal{G}_{\text{pl}}$):
the likelihood has a closed form, optimised numerically via `nlm()` with analytical gradients.
For constrained nonparametric priors: delegates to `ashr` via mix-SQP.
For NPMLE: `mixsqp` directly.

**Step 2 — Compute posterior quantities.** Given $\hat{g}$, compute summaries of:

$$p(\theta_i \mid x_i, s_i, \hat{g}) \propto p(x_i \mid \theta_i, s_i)\, \hat{g}(\theta_i).$$

For parametric priors, posteriors have analytical forms (mixtures of Gaussians, truncated
distributions, etc.).

### 3.3 Outputs of `ebnm()`

+---------------------+--------------------------------------------+
| Output              | Description                                |
+=====================+============================================+
| `fitted_g`          | The estimated prior $\hat{g}$              |
+---------------------+--------------------------------------------+
| `log_likelihood`    | $\log L(\hat{g})$                          |
+---------------------+--------------------------------------------+
| `posterior$mean`    | Posterior means $\hat{\theta}_i$           |
+---------------------+--------------------------------------------+
| `posterior$sd`      | Posterior standard deviations              |
+---------------------+--------------------------------------------+
| `posterior$lfsr`    | Local false sign rates                     |
+---------------------+--------------------------------------------+

**S3 methods:** `coef()`, `vcov()`, `fitted()`, `residuals()`, `logLik()`, `confint()`,
`predict()`, `simulate()`.

### 3.4 Local false sign rates (lfsr)

$$\text{lfsr}(i) := \min\!\left\{ P(\theta_i \le 0 \mid x_i, s_i, \hat{g}),\; P(\theta_i \ge 0 \mid x_i, s_i, \hat{g}) \right\}$$

This is the minimum probability that we are wrong about the *sign* of $\theta_i$.
$\text{lfsr}(i) = 0.05$ means at most a 5% chance that $\theta_i$ has the wrong sign.

**Connection to multiple testing.** The lfsr is the Bayesian analogue of the *local false
discovery rate* (local FDR, Efron 2004): both are posterior probabilities that a discovery
is spurious. Unlike a p-value, lfsr does not require a multiple-testing correction — the
pooling of information via $\hat{g}$ already accounts for multiplicity. If you are
familiar with `ashr` output (`localfsr` column), the lfsr in `ebnm` is the same quantity.

---

## Part IV — Empirical Bayes Matrix Factorization (`flashier`)

### 4.1 The EBMF model

**Convention used throughout this document:** rows index observations (cells, spots, or
samples); columns index features (genes or proteins). Explicitly:

- $\mathbf{X} \in \mathbb{R}^{n \times p}$ — data matrix; $n$ cells/spots (rows), $p$ genes (columns)
- $\mathbf{L} \in \mathbb{R}^{n \times K}$ — loading matrix; row $i$ gives cell $i$'s membership in each of $K$ programs
- $\mathbf{F} \in \mathbb{R}^{p \times K}$ — factor matrix; row $j$ gives gene $j$'s weight in each program
- $\mathbf{E} \in \mathbb{R}^{n \times p}$ — noise; $e_{ij} \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0, \sigma^2)$

The EBMF model is:

$$\mathbf{X} = \mathbf{L}\mathbf{F}^\top + \mathbf{E}$$

Independent EB priors are placed on each column of $\mathbf{L}$ and $\mathbf{F}$:

$$\ell_{ik} \sim g_\ell^{(k)} \in \mathcal{G}_\ell, \qquad f_{jk} \sim g_f^{(k)} \in \mathcal{G}_f,$$

where $\mathcal{G}_\ell$ and $\mathcal{G}_f$ are specified prior families.

> **Naming note.** The term *loading* is used inconsistently across the literature. In this
> document **L** ($n \times K$) contains **cell/sample loadings** — the degree to which each
> observation participates in each program. Some authors (e.g. mixOmics, MOFA)
> call these *scores*, reserving *loadings* for the feature-weight matrix. The
> `flashier` and GBCD convention: **L = cell/sample loadings (participation scores);
> F = gene/feature weights (factors).** The DIVAS convention below will be different.


### 4.2 Algorithm: from K=0 to convergence

**Intuition.** The `flash()` algorithm builds the factorization incrementally, adding one
program at a time, then polishing all programs jointly. 

**Step 0 — Empty model.** Start with $K = 0$: no programs. The residual matrix is
$\mathbf{R} = \mathbf{X}$.

**Step 1 — Greedy addition of factor 1.** Find the best rank-1 approximation to $\mathbf{R}$:

1. Initialise $\boldsymbol{\ell}_1$ and $\mathbf{f}_1$ from the leading left and right singular vectors of $\mathbf{R}$.
2. Alternately update $\boldsymbol{\ell}_1$ (holding $\mathbf{f}_1$ fixed) via EBNM, then update $\mathbf{f}_1$ (holding $\boldsymbol{\ell}_1$ fixed) via EBNM.
3. Iterate until the ELBO (evidence lower bound) converges.
4. Accept factor 1 if it improves the ELBO; otherwise stop.

**Step 2 — Continue greedily.** Update the residual $\mathbf{R} \leftarrow \mathbf{X} - \boldsymbol{\ell}_1\mathbf{f}_1^\top$
and repeat Step 1 to find factor 2. Continue until no new factor improves the ELBO, or
`greedy_Kmax` is reached. This gives an **automatic choice of $K$** — no need to specify
it in advance.

**Step 3 — Backfit.** All $K$ factors were added greedily without looking back. Backfitting
corrects this: cycle through factors $k = 1, \ldots, K$, each time computing partial residuals
(removing all other factors' contributions), then re-estimating $\boldsymbol{\ell}_k$ and
$\mathbf{f}_k$ via EBNM. Repeat until convergence.

**Step 4 — Null check.** After backfitting, remove any factor $k$ whose deletion improves
the ELBO (it adds no real signal). This is the final automatic rank selection step.

```
  K=0 ──greedy──► K=1 ──greedy──► K=2 ──···──► K*
                                                 │
                                              backfit
                                                 │
                                           K* (refined) ──nullcheck──► K** (final)
```

[To Claude: need a small paragraph here explains why the problem of estimating l or f
vector in Step 1 was an EBNM problem: the EBNM problem we talked about about was slightly
different from the current one. The current one is a regression problem, not just mean.
I see that you explain this in more details in 4.3. But a brief introduction here will
be helpful.]

**Why this works for omics.** The greedy stage identifies dominant programs first. The
backfit refines each program in the context of all others, preventing double-counting. The
null check removes spurious factors, making the final $K$ entirely data-driven. In practice,
`flash()` on log1p-normalised data typically converges in 5–30 minutes for $n \le 10{,}000$
spots and $p \le 30{,}000$ genes.

### 4.3 From EBNM to matrix factorization — the bridge *(Optional)*

**Why EBMF reduces to repeated EBNM.** Fix all factors except the $k$-th, and consider
estimating row $i$ of $\mathbf{L}$ (i.e., $\ell_{ik}$). Define the *partial residual*:

$$r_{ij}^{(k)} := x_{ij} - \sum_{k' \ne k} \hat{\ell}_{ik'} f_{jk'},$$

which removes the contribution of all other programs from observation $(i,j)$. Under the
model, $r_{ij}^{(k)} = \ell_{ik} f_{jk} + e_{ij}$, so the partial residuals for cell $i$
look like simple linear regression on the gene weights $f_{jk}$.

Note that a sufficient statistic for $\ell_{ik}$ is
$$\tilde{x}_{ik} := \frac{\sum_{j=1}^p f_{jk}\, r_{ij}^{(k)}}{\sum_{j=1}^p f_{jk}^2},$$
since
$$
\tilde{x}_{ik} \mid \ell_{ik} \sim \mathcal{N}(\ell_{ik},\, \tilde{s}_{ik}^2), \quad \text{with} \quad
\tilde{s}_{ik} := \frac{\sigma}{\sqrt{\sum_{j=1}^p f_{jk}^2}}.
$$

This is
exactly an **EBNM problem**: observe $\tilde{x}_{ik}$ with known standard error $\tilde{s}_{ik}$,
estimate $\ell_{ik}$ under prior $g_\ell^{(k)}$.

The same argument applies symmetrically to each column $k$ of $\mathbf{F}$ (fixing
$\mathbf{L}$). **Each iteration of EBMF calls `ebnm` $K \times (n + p)$ times** — once per
factor per cell (for $\mathbf{L}$) and once per factor per gene (for $\mathbf{F}$). This
is why the speed and flexibility of `ebnm` are critical to the `flashier` engine.

**Informally:** EBMF "asks each cell $i$, for each program $k$: given what you know about the
gene weights ($\mathbf{f}_k$) and the residuals after removing all other programs, how much
of program $k$ does this cell express?" The EBNM prior $g_\ell^{(k)}$ regularises the answer,
enforcing sparsity or non-negativity as specified.

### 4.4 Inference via variational EM *(Optional)*

Exact inference in EBMF is intractable because posteriors of $\ell_{ik}$ and $f_{jk}$
interact. EBMF uses *mean-field variational inference*: the joint posterior is approximated by
a product of independent distributions,

$$p(\mathbf{L}, \mathbf{F}, g_\ell, g_f \mid \mathbf{X}) \approx \prod_{k,i} q(\ell_{ik}) \prod_{k,j} q(f_{jk}),$$

and the ELBO (Evidence Lower Bound) is maximised over $q(\cdot)$, $g_\ell^{(k)}$, and
$g_f^{(k)}$. Each update of $q(\ell_{ik})$ is an EBNM posterior computation (Section 4.3).

### 4.5 Effect of prior family on the decomposition

+------------------------------------------------------+------------------------------------------+
| Prior on $\mathcal{G}_\ell$ and $\mathcal{G}_f$     | Resulting decomposition                  |
+======================================================+==========================================+
| Normal                                               | Similar to (truncated) SVD / PCA         |
+------------------------------------------------------+------------------------------------------+
| Point-normal (both)                                  | Sparse matrix factorization              |
+------------------------------------------------------+------------------------------------------+
| Point-exponential (both)                             | Non-negative sparse MF (NMF)             |
+------------------------------------------------------+------------------------------------------+
| Point-exponential ($\mathcal{G}_\ell$) +             | Semi-NMF: non-negative loadings,         |
| point-normal ($\mathcal{G}_f$)                       | signed factors                           |
+------------------------------------------------------+------------------------------------------+
| Generalised binary ($\mathcal{G}_\ell$) +            | GBCD (see Part VI)                       |
| point-Laplace ($\mathcal{G}_f$)                      |                                          |
+------------------------------------------------------+------------------------------------------+

### 4.6 The `flash()` function (flashier v1.0.7)

```r
fit <- flash(
  data        = X,                        # n x p matrix (rows = cells, columns = genes)
  ebnm_fn     = ebnm_point_exponential,   # prior family; or list c(G_l, G_f) for asymmetric
  greedy_Kmax = 30,                       # maximum factors to add greedily
  backfit     = TRUE,
  nullcheck   = TRUE,
  verbose     = 1L
)
```

**Key outputs:**

+------------------+-------------+----------------------------------------------------------+
| Field            | Dimensions  | Description                                              |
+==================+=============+==========================================================+
| `fit$L_pm`       | $n\times K$ | Posterior means of $\mathbf{L}$: rows = cells/spots,    |
|                  |             | columns = programs                                       |
+------------------+-------------+----------------------------------------------------------+
| `fit$F_pm`       | $p\times K$ | Posterior means of $\mathbf{F}$: rows = genes,          |
|                  |             | columns = programs                                       |
+------------------+-------------+----------------------------------------------------------+
| `fit$n_factors`  | scalar      | Number of programs retained ($K$, chosen automatically)  |
+------------------+-------------+----------------------------------------------------------+
| `fit$pve`        | length $K$  | Proportion of variance explained by each program         |
+------------------+-------------+----------------------------------------------------------+

**Note on input orientation.** Data in Bioconductor (SpatialExperiment, SingleCellExperiment)
is stored as genes × cells (features in rows). `flashier` expects cells × genes (observations
in rows). Always transpose before calling `flash()`:

```r
lc_matrix <- log1p(assay(sce, "counts"))   # p x n (genes x cells)
X          <- Matrix::t(lc_matrix)          # n x p (cells x genes) — required by flash()
```

### 4.7 Interpreting NMF output

After fitting, the key objects are `fit$L_pm` (cell/spot loadings, $n \times K$) and
`fit$F_pm` (gene weights, $p \times K$).

**Identifying the top genes in a program:**

```r
k <- 1                                          # program of interest
gene_ranks <- order(fit$F_pm[, k], decreasing = TRUE)
top_genes  <- rownames(fit$F_pm)[gene_ranks[1:20]]
```

Genes with high $F_{jk}$ values define program $k$ — they are the most up-regulated
features in cells/spots that strongly express this program.

**Reading `fit$pve`:** the proportion of variance explained (PVE) per program. The interpretation
is similar to PVE in PCA, when NMF factors are *orthogonal*. When the factors are not orthogonal,
the sum of all proportions could exceed 1, since there are overlaps between factors.

**What the loadings in `fit$L_pm` represent:** each column $\ell_k \in \mathbb{R}^n$
is a sample score vector: how strongly each cell or spot expresses program $k$. For spatial
transcriptomics, plot $\ell_k$ on the tissue image (colour spots by $L_{ik}$) to see
the spatial pattern of program $k$. For single-cell data, can colour a UMAP by $L_{ik}$.

**Units of `fit$F_pm` when input is log-scaled counts.** For example, when data is
log1p-transformed, we have $X_{ij} = \log(1 + c_{ij})$
where $c_{ij}$ is the raw count. In this case the EBMF model decomposes log-space counts. The entry
$F_{jk}$ approximately represents the increment in $\log(1+\text{count})$ for gene $j$
that comes from a one-unit increase in loading $L_{ik}$. For moderate-to-high counts
$\log(1+c) \approx \log c$, so $F_{jk}$ approximates a log-fold change in expression
between a cell with $L_{ik} = 1$ and a cell with $L_{ik} = 0$ (all else equal).

**When to merge or split programs.** We can merge two programs $k$ and $k'$ if their loading
vectors are highly correlated (e.g. when $\text{cor}(\ell_k, \ell_{k'}) > 0.9$) and their gene
weights overlap substantially. We can also split a program if its loading map shows two spatially
or biologically distinct regions that differ only in intensity. Merging and splitting may also 
happen when we reduce or increase the number of factors.

---

## Part V — Topic Models and Their Relationship to NMF *(Optional)*

***Warning: This whole part is entirely AI-generated based on references, and hasn't been checked yet.***

The topic model literature — originally developed for text documents — provides a
complementary view of matrix factorization that is widely used in single-cell genomics.
The `fastTopics` R package implements this framework.

### 5.1 The multinomial topic model

**Data:** $\mathbf{X} \in \mathbb{Z}_{\ge 0}^{n \times p}$, an integer count matrix with
$n$ cells (rows) and $p$ genes (columns). Entry $x_{ij}$ is the UMI count of gene $j$
in cell $i$; row total $s_i = \sum_j x_{ij}$ is the library size.

**Model:** Each row of $\mathbf{X}$ is drawn from a multinomial:

$$x_{i1}, \ldots, x_{ip} \mid \mathbf{L}, \mathbf{F} \sim \text{Multinomial}\!\left(s_i;\; \pi_{i1}, \ldots, \pi_{ip}\right), \qquad \pi_{ij} = \sum_{k=1}^K l_{ik} f_{jk},$$

where:

- $\mathbf{L} \in \mathbb{R}_+^{n \times K}$: topic membership matrix, $l_{ik} \ge 0$, rows sum to 1 ($\sum_k l_{ik} = 1$ for all $i$)
- $\mathbf{F} \in \mathbb{R}_+^{p \times K}$: topic-gene frequency matrix, $f_{jk} \ge 0$, columns sum to 1 ($\sum_j f_{jk} = 1$ for all $k$)

The "topics" (columns of $\mathbf{F}$) are distributions over genes; $l_{ik}$ is the
proportion of cell $i$'s expression attributable to topic $k$. Unlike hard clustering
($l_{ik} \in \{0,1\}$), this is a *graded* or *mixed-membership* model.

### 5.2 Equivalence to Poisson NMF

The **Poisson NMF** model replaces the multinomial with independent Poisson observations:

$$x_{ij} \mid \mathbf{H}, \mathbf{W} \sim \text{Pois}\!\left((\mathbf{HW}^\top)_{ij}\right), \quad h_{ik} \ge 0, \; w_{jk} \ge 0.$$

Here $\mathbf{H} \in \mathbb{R}_+^{n \times K}$ corresponds to cell memberships and
$\mathbf{W} \in \mathbb{R}_+^{p \times K}$ to gene-topic weights — the same roles as
$\mathbf{L}$ and $\mathbf{F}$ in the multinomial model, but without the sum-to-one normalisation.

**Lemma 1 (Carbonetto et al. 2022).** The two likelihoods are proportional:

$$p_{\text{PNMF}}(\mathbf{X} \mid \mathbf{H}, \mathbf{W}) = p_{\text{MTM}}(\mathbf{X} \mid \mathbf{L}, \mathbf{F}) \cdot \prod_{i=1}^n \text{Pois}(s_i;\, s_i),$$

where $\mathbf{L}$, $\mathbf{F}$ are obtained from $\mathbf{H}$, $\mathbf{W}$ by
normalising rows of $\mathbf{H}$ and columns of $\mathbf{W}$ to sum to one
(absorbing scale factors into the row totals).

**Consequence:** any algorithm that maximises the Poisson NMF likelihood is simultaneously
maximising the multinomial topic model likelihood. Topic model fits can therefore be computed
using fast NMF optimisation routines.

### 5.3 Algorithms: EM vs. coordinate descent

The classical algorithm for topic models is **EM** (Blei et al. 2003), which is equivalent
to Lee & Seung's multiplicative update rule for NMF:

$$h_{ik} \leftarrow h_{ik} \frac{\sum_j w_{jk}\, x_{ij} / \lambda_{ij}}{\sum_j w_{jk}}, \qquad
w_{jk} \leftarrow w_{jk} \frac{\sum_i h_{ik}\, x_{ij} / \lambda_{ij}}{\sum_i h_{ik}},$$

where $\lambda_{ij} = \sum_k h_{ik} w_{jk}$.

EM can be slow to converge and frequently gets stuck at poor local optima, especially
for large $K$. The `fastTopics` package implements **co-ordinate descent (CD) with
Nesterov extrapolation** (Ang & Gillis 2019), which applies a momentum step after each
CD update:

$$\mathbf{H}^{\text{ext}} \leftarrow P_+\!\left[\mathbf{H}^{\text{new}} + \beta^{(t)}\left(\mathbf{H}^{\text{new}} - \mathbf{H}^{(t-1)}\right)\right],$$

where $P_+$ projects onto the non-negative orthant and $\beta^{(t)} \in [0,1]$ is an
adaptive momentum parameter. Benchmarks show CD reaches substantially better (higher)
log-likelihoods than EM in less time, and on real single-cell data the estimated topics
can be qualitatively different — making the choice of algorithm matter biologically.

### 5.4 Structure plots: visualising membership matrices

The standard visualisation for the topic model $\mathbf{L}$ matrix is the **Structure plot**
(named after the population genetics software): a stacked bar chart with one bar per cell
(row of $\mathbf{L}$), segments coloured by topic. Because rows sum to one, the bars all
reach the same height. Cells are ordered by their dominant topic membership.

```r
library(fastTopics)
fit <- fit_topic_model(X, k = K)           # X is n x p (cells x genes)
structure_plot(fit)                         # stacked bar plot of fit$L
```

This plot makes it immediately visible whether topics are:

- **Discrete cell types** — cells belong almost entirely to one topic (bars are nearly solid)
- **Continuous programs** — cells have graded membership across multiple topics (bars show mixtures)
- **Shared axes** — a topic active at a low level across all cells (thin stripe across all bars)

### 5.5 GoM DE: identifying genes distinctive to each program

After fitting the topic model and obtaining $\mathbf{L}$ (cell memberships), the biological
question becomes: *which genes are distinctive to each topic/program?* Standard DE tools
(DESeq2, MAST) require discrete group assignments and cannot handle mixed membership. **GoM DE**
(Grades-of-Membership Differential Expression; Carbonetto et al. 2023) extends this to
continuous membership.

**Model.** Fix $\mathbf{L}$ from the topic fit. For each gene $j$ separately, model the
counts as:

$$x_{ij} \sim \text{Poisson}\!\left(s_i \sum_{k=1}^K l_{ik}\, p_{jk}\right),$$

where $p_{jk}$ is the (unknown) expression level of gene $j$ in topic $k$. The parameters
$p_{jk}$ are estimated gene-by-gene with $\mathbf{L}$ held fixed — a separate per-gene
Poisson regression.

**Log-fold changes.** For each topic pair $(k, l)$, the pairwise LFC is:

$$\text{LFC}_{k,l}(j) := \log_2\!\frac{p_{jk}}{p_{jl}}.$$

A gene is truly *distinctive* for topic $k$ only if it is differentially expressed relative
to **all** other topics, not just one. The **least-extreme LFC** captures this:

$$\text{LFC}_k^{\text{l.e.}}(j) := \text{LFC}_{k,\, l^*}(j), \qquad
l^* = \underset{l \ne k}{\arg\min}\left|\text{LFC}_{k,l}(j)\right|.$$

A large positive $\text{LFC}_k^{\text{l.e.}}(j)$ means gene $j$ is higher in topic $k$ than
in its most similar other topic — making it a reliable marker.

**Shrinkage and uncertainty.** LFC estimates are noisy, especially for lowly expressed genes.
GoM DE applies `ashr` (which is EBNM with a normal scale mixture prior) to shrink the
LFC estimates toward zero and compute **lfsr** (local false sign rate) for each gene-topic pair.
The lfsr is the probability that the sign of $\text{LFC}_k^{\text{l.e.}}(j)$ is wrong — the
same quantity defined in Part III, now computed from MCMC posteriors.

```r
library(fastTopics)
de <- de_analysis(fit, X)                       # GoM DE with ashr shrinkage
volcano_plot(de, k = 1)                         # volcano: LFC vs lfsr for topic 1
top_genes <- head(order(de$lfsr[, 1]), 20)      # top distinctive genes for topic 1
```

### 5.6 `flashier` vs. topic model — comparison

Both methods decompose count-derived matrices into $\mathbf{L}$ and $\mathbf{F}$ with non-negativity,
but they differ in important ways:

+-----------------------------+---------------------------+----------------------------+
| Feature                     | Topic model (`fastTopics`)| `flashier` EBNMF           |
+=============================+===========================+============================+
| Likelihood                  | Poisson (or multinomial)  | Gaussian (log-transformed) |
+-----------------------------+---------------------------+----------------------------+
| $K$ selection               | Pre-specified             | Automatic (greedy ELBO)    |
+-----------------------------+---------------------------+----------------------------+
| Sparsity                    | Via sum-to-one            | EB prior (point-exp.)      |
+-----------------------------+---------------------------+----------------------------+
| Row normalisation of L      | Rows sum to 1             | No constraint              |
+-----------------------------+---------------------------+----------------------------+
| Uncertainty                 | MLE only                  | Posterior SDs and lfsr     |
+-----------------------------+---------------------------+----------------------------+
| Gene annotation             | GoM DE (Poisson reg.)     | F factors; lfsr available  |
+-----------------------------+---------------------------+----------------------------+
| Best for                    | Count data, interpretable | Log-count data, auto K,    |
|                             | membership proportions    | formal uncertainty         |
+-----------------------------+---------------------------+----------------------------+

**When to use which:** use `fastTopics` when you want membership proportions (how much of
each cell is each topic) and intend to run GoM DE for gene annotation. Use `flashier` when
you want automatic $K$ selection, formal uncertainty quantification on the factors, or
non-count data (e.g., spatial smoothed PCs, log-normalised proteomics).

---

## Part VI — Generalized Binary Covariance Decomposition (GBCD)

### 6.1 Overview and motivation

GBCD (Liu, Carbonetto, Willwerscheid, Oakes, Macleod & Stephens, *Nature Genetics* 2025)
is a matrix factorization method for *multi-patient or multi-sample* single-cell RNA-seq.
The core challenge: single-cell data from multiple tumors is dominated by *intertumor
heterogeneity*: cell expression profiles differ so much between patients that standard NMF
or PCA finds patient-specific components first, obscuring shared biological programs.

GBCD addresses this with two key structural assumptions on $\mathbf{Y} \approx \mathbf{L}\mathbf{F}^\top$ ($n$ cells × $p$ genes):

1. **Generalised binary memberships:** $l_{ik}$ are non-negative and nearly binary (close to 0 or 1), encouraging each cell to be "in" or "out" of each gene expression program (GEP).

2. **Orthogonal GEP signatures:** $\mathbf{F}^\top\mathbf{F} \approx D$ (diagonal), preventing shared programs from being absorbed into patient-specific components.

Together these imply a sample-sample covariance decomposition:

$$\mathbf{Y}\mathbf{Y}^\top \approx \mathbf{L}\mathbf{F}^\top\mathbf{F}\mathbf{L}^\top = \mathbf{L}D\mathbf{L}^\top.$$

**NB**: With our convention $\mathbf{Y}\mathbf{Y}^\top$ is an $n\times n$ symmetric matrix of
pairwise similarities between samples/cells.


### 6.2 Data preparation

For UMI count matrix $\mathbf{X}$ ($n$ cells × $p$ genes), with size factor $s_i = \sum_j x_{ij}$
and median $\hat{s} = \text{median}(s_i)$, apply the shifted log-count normalisation
(pseudocount $c = 0.1$):

$$y_{ij} = \log\!\left(c + \frac{\hat{s}}{s_i} x_{ij}\right) - \log c. \tag{3}$$

This ensures $y_{ij} \ge 0$ and sparsity-preserving: $y_{ij} = 0$ iff $x_{ij} = 0$.

### 6.3 Stage 1 — Estimating GEP memberships $\mathbf{L}$

**Covariance decomposition model:**

$$\mathbf{Y}\mathbf{Y}^\top = \mathbf{L}\mathbf{L}^\top + \epsilon I_n + \mathbf{E}. \tag{4}$$

**Generalised binary (GB) prior:**

$$l_{ik} \sim \left(1 - \pi_k^\ell\right)\delta_0 + \pi_k^\ell\, \mathcal{N}_+(\mu_k, \sigma_k^2). \tag{5}$$

Here $\delta_0$ is a point mass at zero. The truncated normal slab produces 
memberships near 0 or near $\mu_k$,  when the ratio $\sigma_k^2 / \mu_k$ is set small,
hence the name "generalised binary". 

**Relaxed model for tractability:** 

$$\mathbf{Y}\mathbf{Y}^\top = \mathbf{L}\tilde{\mathbf{L}}^\top + \epsilon I_n + \mathbf{E}. \tag{6}$$

$\mathbf{L}$ and $\tilde{\mathbf{L}}$ are fitted independently. After fitting, retain
component $k$ only if $\text{cor}(\mathbf{l}_k, \tilde{\mathbf{l}}_k) > 0.8$.

### 6.4 Stage 2 — Estimating GEP signatures $\mathbf{F}$

Fix $\mathbf{L}$ from Stage 1. Fit:

$$\mathbf{Y} = \mathbf{L}\mathbf{F}^\top + \tilde{\mathbf{E}}, \qquad \tilde{e}_{ij} \sim \mathcal{N}(0, \varphi_j^2). \tag{7-8}$$

with gene-specific variances $\varphi_j^2$ estimated jointly. **Point-Laplace prior on $\mathbf{F}$:**

$$f_{jk} \sim \left(1 - \pi_k^f\right)\delta_0 + \pi_k^f\, \text{Laplace}(\lambda_k). \tag{9}$$

This encourages most $f_{jk}$ to be zero (genes not differentially expressed in GEP $k$)
while allowing a subset to be positive (up-regulated) or negative (down-regulated). The
entry $f_{jk}$ approximates the log-fold change of gene $j$ in cells fully in GEP $k$
relative to cells with no membership in GEP $k$.

### 6.5 Quantifying uncertainty

GBCD additionally: (1) fits gene-level linear regression with $\mathbf{L}$ as regressors
to obtain LFC estimates and standard errors; (2) applies `ashr` to improve calibration;
(3) computes posterior z-scores and local false sign rates.

### 6.6 Annotating GEPs

[To Claude: need to explain the notation below; also elaborate on the intuition.]

**Patient-specificity score** ($\upsilon_k \in [0,1]$, high = patient-specific):

$$\upsilon_k = 1 - \frac{\max_{m \ne m(k)} p_{m,k}}{p_{m(k),k}}. \tag{14}$$


### 6.7 Choosing $K$ and computational considerations

GBCD requires an upper bound $K_{\max}$; the ELBO-based null check removes uninformative
components automatically. GEPs found with smaller $K$ are generally a subset of those
found with larger $K$. The computational bottleneck is forming $\mathbf{Y}\mathbf{Y}^\top$
($O(n^2 p)$ operations); for large $n$ (e.g., 35,000 PDAC cells), an alternative
implementation works directly with the sparse matrix $\mathbf{Y}$.

### 6.8 R code example (do not run)

```r
library(gbcd)
library(flashier)
library(ebnm)

# 1. Shifted log-count normalisation (n cells x p genes)
s_i   <- rowSums(X_counts)
s_hat <- median(s_i)
Y <- log(0.1 + (s_hat / s_i) * X_counts) - log(0.1)   # Y >= 0

# 2. Fit GBCD
res <- fit_gbcd(
  Y        = Y,
  Kmax     = 20,
  prior    = ebnm::ebnm_generalized_binary,
  maxiter1 = 500, maxiter2 = 200, maxiter3 = 500,
  verbose  = 1
)

# 3. Results
L      <- res$L        # n x K cell memberships (binary-like)
F_lfc  <- res$F$lfc    # p x K signed LFC per gene per GEP
F_lfsr <- res$F$lfsr   # p x K local false sign rates

# 4. Inspect GEP k
k      <- 1
top_up <- rownames(F_lfc)[order(F_lfc[, k], decreasing = TRUE)[1:20]]
sig    <- rownames(F_lfsr)[F_lfsr[, k] < 0.05]

# 5. Membership heatmap
pheatmap::pheatmap(
  L,
  cluster_cols   = FALSE,
  annotation_row = cell_meta,
  color          = colorRampPalette(c("gray96", "red"))(50),
  show_rownames  = FALSE
)
```

[To Claude: add a new Part here to run different versions of NMF on COVID-19
two time point proteomics data. Pick a few important graphics to plot. Introduce
the data here. Also include some interpretation of the results.]

---

## Part VII — Multi-Block Data: From NMF Stacking to DIVAS

Parts IV--VI showed how to decompose a single matrix $\mathbf{X}$ into programs.
But modern experiments often produce **multiple data blocks** measured on the same
samples. For example, proteomics and metabolomics on the same patients, or the
same assay at two time points. How should we handle this?

### 7.1 The multi-block problem

[To Claude: since we will introduce the COVID data above, here only say why
this data can be viewed as two-block and what will be the benefits.]

Consider 120 COVID-19 patients measured at two time points (T1 = admission,
T2 = follow-up) with Olink proteomics (481 proteins). We want to know:

1. Which proteomics programs are **stable** across time points (present at both T1 and T2)?
2. Which programs are **specific** to admission (acute responses that resolve)?
3. Which programs are **specific** to follow-up (recovery processes that emerge later)?

NMF on a single matrix cannot directly answer these questions. Two common workarounds
involve different ways of assembling the data matrix before running NMF.

### 7.2 Stacking strategies for NMF

[To Claude: this part should be moved to the above newly created section on applying NMF to
COVID data.]

**Sample stacking.** Stack T1 and T2 observations as rows of a single matrix:

$$\mathbf{X}_\text{stacked} = \begin{pmatrix} \mathbf{X}_{T1} \\ \mathbf{X}_{T2} \end{pmatrix} \in \mathbb{R}^{(2n) \times p}$$

Each factor $k$ has loadings $\ell_{ik}$ for every patient-timepoint pair. To ask whether
factor $k$ is time-dependent, one can compare $\ell_{ik}^{T1}$ vs. $\ell_{ik}^{T2}$ post
hoc --- but this is informal. There is no built-in mechanism to determine whether a
factor is shared across time points or specific to one.

**Feature stacking.** Concatenate different modalities as columns:

$$\mathbf{X}_\text{feat} = \begin{pmatrix} \mathbf{X}_\text{prot} & \mathbf{X}_\text{metab} \end{pmatrix} \in \mathbb{R}^{n \times (p_1 + p_2)}$$

A factor that loads on both protein and metabolite columns captures coordinated
cross-modal programs. But again, there is no formal test for whether a factor is
truly shared vs. driven by one modality.

### 7.3 Limitations of the stacking approach

In the COVID-19 case study, EBNMF on the sample-stacked proteomics matrix
(240 samples x 481 proteins) found K = 14 programs. The best clinical association
was factor F5 (Spearman $\rho = -0.344$ with WHO improvement, FDR = 0.054).
These are useful results, but they leave important questions unanswered:

1. **Shared vs. individual programs are not identified.** Every factor spans both
   time points. Post-hoc comparison of T1 vs. T2 loadings can suggest which factors
   are time-stable, but this is descriptive, not inferential.

2. **The number of shared and individual components is not estimated.** NMF (including
   EBNMF) selects a total K but does not partition it into $K_\text{shared}$,
   $K_\text{T1}$, and $K_\text{T2}$.

3. **Time-specific programs may be obscured.** A protein program that is active only at
   admission will be diluted across all 240 rows, competing for variance against
   time-stable programs. Weak time-specific signals can be absorbed into larger shared
   factors or split across multiple small factors, making them difficult to interpret.

### 7.4 DIVAS: Data Integration Via Analysis of Subspaces

DIVAS (Prothero et al. 2024; Sun et al. 2026) was developed specifically for the
multi-block setting. Given $K$ data blocks measured on the same $n$ samples, DIVAS
decomposes each block's signal subspace into components that are **jointly shared**
(present in all blocks), **partially shared** (present in a subset of blocks), or
**individual** (present in one block only).

The key idea is geometric. Each data block $\mathbf{X}_k \in \mathbb{R}^{d_k \times n}$
has a signal subspace in trait space $\mathbb{R}^n$. DIVAS measures the relationships
between these subspaces via **principal angles** --- the canonical angles between pairs
of subspaces. Small angles indicate shared structure; large angles indicate independent
variation. A rotational bootstrap provides statistical inference on which angles are
significantly small (i.e. which components are genuinely shared vs. coincidentally
aligned).

For each component, DIVAS returns:
- A **score vector** $\mathbf{s} \in \mathbb{R}^n$ (shared across the relevant blocks),
  analogous to a factor loading in NMF --- it tells you which samples participate.
- A **loading vector** $\boldsymbol{\ell}_k \in \mathbb{R}^{d_k}$ per block, analogous
  to a gene/protein weight --- it tells you which features drive the component in
  each block.

The decomposition is:

$$\mathbf{A}_k = \sum_{\mathbf{i} \ni k} \mathbf{L}_{\mathbf{i},k} \, \mathbf{S}_\mathbf{i}^\top$$

where $\mathbf{i}$ indexes subsets of blocks, $\mathbf{S}_\mathbf{i}$ are the score
vectors shared by blocks in subset $\mathbf{i}$, and $\mathbf{L}_{\mathbf{i},k}$ are
the corresponding loadings for block $k$.

**Key advantages over NMF stacking:**

| | NMF on stacked matrix | DIVAS |
|---|---|---|
| Shared vs. individual | Not distinguished | Explicitly identified and labelled |
| Number of each type | Single total K | $K_\text{shared}$, $K_\text{indiv,1}$, $K_\text{indiv,2}$ estimated separately |
| Statistical inference | None on component identity | Rotational bootstrap p-values for shared structure |
| Input structure | Single matrix (stacking choice is the user's) | Separate blocks; no stacking decision needed |
| Sparsity / prior | EB priors (flashier), regularisation (standard NMF) | SVD-based; no sparsity prior |

[To Cladue: since we will include this reference at the beginning of document,
remove the below text chunk.]
> **References:**
> - Prothero J, Jiang M, Hannig J, Tran-Dinh Q, Ackerman A, Marron JS (2024).
>   "Data integration via analysis of subspaces (DIVAS)." *TEST*.
>   doi:10.1007/s11749-024-00923-z
> - Sun Y, Marron JS, Le Cao KA, Mao J (2026). "DIVAS: an R package for
>   identifying shared and individual variations of multiomics data." *bioRxiv*.
>   doi:10.64898/2026.01.12.698985

### 7.5 COVID-19 case study: proteomics T1 vs. T2

[To Claude, include code for running DIVAS and plotting. Then interprete the results.]

We applied DIVAS to the same COVID-19 proteomics data used in the NMF analyses,
treating T1 and T2 as two blocks (each 481 proteins x 120 patients, Olink NPX log2,
row-centred). DIVAS identified **24 components**: 11 shared (T1+T2), 7 T1-individual,
and 6 T2-individual.

#### Shared components and severity

The shared component **T1+T2-2** showed the strongest severity association of any
component in the entire analysis (Spearman $\rho = 0.730$ with T1 WHO score,
FDR $< 0.0001$). This is a time-stable severity program --- the same proteomics axis
separates mild from severe patients at both admission and follow-up. Top driving
proteins include CGRP, CCL7, CCDC80, NAD kinase, and KRT19 (positive loadings) versus
IL-5 isoforms and FABP2 (negative loadings). The strength of this association
($\rho = 0.730$) substantially exceeds the best NMF result ($\rho = 0.344$ for EBNMF
F5), illustrating how separating shared from individual variation sharpens the
severity signal.

#### Individual components: what NMF cannot easily reveal

The 13 individual components (7 T1, 6 T2) are the distinctive contribution of DIVAS.
These capture proteomics variation unique to one time point --- programs that are
present at admission but not follow-up, or vice versa. None of these have strong
severity associations (all FDR > 0.4 except T2-Individual-1 at FDR = 0.08), which is
itself informative: it confirms that severity is predominantly a time-stable phenomenon.
The value of individual components lies in revealing **time-specific biology**.

**T1-individual (admission-specific) themes:**

- **Acute multi-organ injury (T1-Indiv-2).** Dominated by troponin I (+0.272),
  HAVCR1/TIM-1 (kidney injury marker, +0.200), ACE2 (+0.155), and HAO1 (liver
  peroxisomal enzyme, +0.227). This is a cardiac/hepatic/renal damage axis present at
  admission but absent at follow-up --- patients either recovered or this pattern merged
  into the shared severity axis by T2.

- **Th1 antiviral activation (T1-Indiv-4).** IFN-gamma dominates (+0.315) with
  IFN-gamma-induced chemokines CXCL9, CCL19, and CCL20. IL-10 co-loads positively,
  reflecting the well-known regulatory coupling (IL-10 is induced alongside IFN-gamma as
  a brake). This acute antiviral response is admission-specific --- it resolves or
  transforms by follow-up.

- **Cardiovascular damage (T1-Indiv-6).** SERPINA12/vaspin (-0.338), troponin I
  (-0.295), and BNP (-0.237) form the clearest cardiac injury axis. MMP12 (macrophage
  elastase) adds a tissue destruction component. This pattern of cardiac biomarker
  elevation at presentation does not persist as a separate axis at T2.

- **Hepatic metabolic stress (T1-Indiv-3).** FGF21 (-0.287) and HAO1 (-0.284) mark
  liver metabolic derangement, distinct from the multi-organ damage of T1-Indiv-2
  because it emphasises metabolic/mitochondrial stress (CA5A, FBP1) rather than
  cardiac/renal injury.

**T2-individual (follow-up-specific) themes:**

- **Persistent cardiac injury (T2-Indiv-1).** Troponin I (-0.207) and BNP (-0.103) at
  follow-up, with the only near-significant delta association ($\rho = -0.246$,
  FDR = 0.08) --- patients with worse cardiac markers at T2 improved less.

- **Th1-to-neutrophil immune transition (T2-Indiv-6).** IFN-gamma (-0.267) and IL-12B
  (-0.130) opposed by CEACAM8/CD66b (+0.151, neutrophil activation), NCF2/p67phox
  (+0.157, NADPH oxidase), and RNase 3/ECP (+0.167, eosinophil degranulation). This
  captures divergence in recovery immune trajectories: some patients shift from adaptive
  (IFN-gamma-driven) to innate (neutrophil-dominated) immunity at follow-up. This axis
  does not exist at T1.

- **Mucosal inflammation (T2-Indiv-5).** CCL20 (+0.154) and oncostatin M (+0.127) mark
  ongoing mucosal/tissue remodelling at follow-up, opposed by GH1 (-0.342) and
  FBP1 (-0.194), suggesting metabolic recovery in the opposite direction.

#### What individual components reveal that NMF does not

The NMF analysis (Part IV) found 14 programs on the stacked matrix, some of which
hinted at time-dependent behaviour (e.g. factors with higher T2 than T1 loadings). But
it could not answer: "Is the IFN-gamma activation program at admission a genuinely
distinct axis, or just a noisy variant of a time-stable program?" DIVAS answers this
directly: T1-Indiv-4 (IFN-gamma) is statistically individual --- it passed the
rotational bootstrap test for being absent from the T2 signal subspace. Similarly, the
neutrophil-dominated immune state (T2-Indiv-6) is confirmed as a follow-up-specific
phenomenon, not a weak reflection of an admission axis.

The biological picture that emerges: at admission, patients vary along axes of acute
organ damage, Th1 antiviral activation, and cardiovascular stress. By follow-up, the
organ damage and Th1 axes have either resolved or been absorbed into time-stable
programs. New axes emerge: neutrophil-dominated immunity, mucosal inflammation, and
persistent cardiac injury in a subset of patients. The time-stable severity program
(T1+T2-2) persists throughout, but the **character** of the disease changes --- and
only DIVAS, not NMF, can formally separate these dynamics.
