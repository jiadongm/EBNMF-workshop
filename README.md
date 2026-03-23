# Empirical Bayes Matrix Factorization Workshop

Materials for a 1-hour mini-workshop on Empirical Bayes Normal Means (EBNM),
Empirical Bayes Matrix Factorization (EBMF), and multi-block extensions.

## Workshop materials

| Material | Format | Description |
|----------|--------|-------------|
| [Workshop slides](https://jiadongm.github.io/EBNMF-workshop/slides_workshop.html) | xaringan slides | 45 min presentation + 15 min hands-on |
| [EBNMF handout](https://jiadongm.github.io/EBNMF-workshop/EBNMF.html) | HTML document | Technical reference covering EBNM, prior families, flashier, and GBCD |
| [DIVAS handout](https://jiadongm.github.io/EBNMF-workshop/DIVAS.html) | HTML document | Multi-block decomposition via DIVAS |

## Topics covered

1. **Empirical Bayes shrinkage** -- connecting limma, edgeR, and ashr to the EBNM framework
2. **The EBNM problem** -- prior families, the `ebnm` R package, and adaptive shrinkage
3. **EBMF and flashier** -- matrix factorization with EB priors; automatic K selection
4. **GBCD** -- Generalized Binary Covariance Decomposition for shared/context-specific programs
5. **DIVAS** -- Data Integration Via Analysis of Subspaces for multi-block data

## Case study

The workshop uses a **COVID-19 multiomics cohort** (Su et al., *Cell* 2020) as a running
example: 120 patients, Olink proteomics at two time points, severity scores 1--7.

## Source files

- `EBNMF.md` -- handout source (Markdown)
- `DIVAS.md` -- DIVAS handout source
- `slides_workshop.Rmd` -- slides source (xaringan/R Markdown)
- `workshopCode.R` -- hands-on R code for participants
- `EBNMF_guide.md` -- flashier + GBCD API reference
- `DIVAS_guide.md` -- DIVAS API reference

## Key R packages

- [`ebnm`](https://cran.r-project.org/package=ebnm) -- Empirical Bayes Normal Means
- [`flashier`](https://cran.r-project.org/package=flashier) -- Empirical Bayes Matrix Factorization
- [`fastTopics`](https://github.com/stephenslab/fastTopics) -- Topic models / Poisson NMF
- [`DIVAS`](https://github.com/jiadongm/DIVAS) -- Multi-block decomposition

## References

- Willwerscheid, Carbonetto & Stephens (2025). *ebnm: an R package for solving the empirical Bayes normal means problem using a variety of prior families.* JSS.
- Liu, Carbonetto et al. (2025). *Dissecting tumor transcriptional heterogeneity from multi-tumor single-cell RNA-seq data.* Nature Genetics.
- Carbonetto, Sarkar, Wang & Stephens (2021). *Non-negative matrix factorization algorithms greatly improve topic model fits.* arXiv.
- Carbonetto, Luo, Sarkar et al. (2023). *GoM DE: interpreting structure in sequence count data with differential expression analysis allowing for grades of membership.* Genome Biology.
- Sun, Marron, Le Cao & Mao (2026). *DIVAS: Data Integration Via Analysis of Subspaces.* bioRxiv.
