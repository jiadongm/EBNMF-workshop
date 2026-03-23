# 08_divas.R — DIVAS: proteomics T1 vs T2 (2-block analysis)
#
# Two data blocks (same 481 proteins, same 120 patients):
#   Block 1 (T1): 481 proteins × 120 patients — Olink NPX (log2) at admission
#   Block 2 (T2): 481 proteins × 120 patients — Olink NPX (log2) at follow-up
#
# DIVAS decomposes into:
#   - T1+T2 shared:  proteomics variation present at both time points (stable programs)
#   - T1-individual: variation unique to time point 1 (early-only programs)
#   - T2-individual: variation unique to time point 2 (late-only programs)
#
# Consistent with 02_NMF_stacked.R and 03_NMF_delta.R which also used proteomics only.

library(DIVAS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(patchwork)

setwd("/data/gpfs/projects/punim0613/JiaDong/NatashaSharma/ns_snail_spatial/EBNMF/CovidCase")

LOAD_CACHED <- file.exists("results/divas_2block.rds")
fig_dir     <- "results/figs/divas"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────

cat("Loading data...\n")
meta     <- readRDS("results/sample_meta.rds")
prot_raw <- read.csv("data/240_input_4omics/proteomics_120patients.csv",
                     row.names = 1, check.names = FALSE)

cat(sprintf("  Proteomics: %d proteins × %d samples\n",
            nrow(prot_raw), ncol(prot_raw)))

# ── 2. Split by time point ────────────────────────────────────────────────────

t1_meta <- meta[meta$time_point == "T1", ] |> arrange(patient_id)
t2_meta <- meta[meta$time_point == "T2", ] |> arrange(patient_id)
stopifnot(all(t1_meta$patient_id == t2_meta$patient_id))

patient_ids <- as.character(t1_meta$patient_id)

# Blocks: proteins × patients (DIVAS convention: features in rows, samples in columns)
Prot_T1 <- as.matrix(prot_raw[, t1_meta$sample])
Prot_T2 <- as.matrix(prot_raw[, t2_meta$sample])
colnames(Prot_T1) <- colnames(Prot_T2) <- patient_ids

data_list <- list(T1 = Prot_T1, T2 = Prot_T2)

cat(sprintf("  T1 block: %d × %d\n", nrow(Prot_T1), ncol(Prot_T1)))
cat(sprintf("  T2 block: %d × %d\n", nrow(Prot_T2), ncol(Prot_T2)))

# ── 3. Run DIVAS ──────────────────────────────────────────────────────────────

if (!LOAD_CACHED) {
  cat("\nRunning DIVAS (2 blocks, nsim=400)...\n")
  t0 <- proc.time()

  set.seed(123)
  divas_res <- DIVASmain(
    datablock    = data_list,
    nsim         = 400,
    colCent      = T,   # centre proteins across patients within each block
    rowCent      = F,       
    seed         = 123,
    ReturnDetail = TRUE
  )

  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("DIVAS completed in %.1f min.\n", elapsed / 60))

  saveRDS(divas_res, "results/divas_2block.rds")
  cat("Saved: results/divas_2block.rds\n")
} else {
  cat("Loading cached DIVAS results...\n")
  divas_res <- readRDS("results/divas_2block.rds")
}

# ── 4. Inspect results ────────────────────────────────────────────────────────

cat("\n=== DIVAS Results ===\n")
scores_mat <- divas_res$Scores
K_total    <- ncol(scores_mat)
comp_names <- colnames(scores_mat)
cat(sprintf("Total components found: %d\n", K_total))
cat("Component names:\n")
for (nm in comp_names) cat(sprintf("  %s\n", nm))

# Categorise: T1+T2 (shared), T1-individual, T2-individual
# DIVAS naming: "T1+T2-k" for shared, "T1-Individual-k" / "T2-Individual-k" for individual
comp_df <- data.frame(
  component = comp_names,
  stringsAsFactors = FALSE
)
comp_df$category <- sapply(comp_names, function(nm) {
  if (grepl("^T1\\+T2-", nm))       "Shared (T1+T2)"
  else if (grepl("^T1-", nm))       "T1-individual"
  else if (grepl("^T2-", nm))       "T2-individual"
  else                               "Other"
})

cat("\nComponent categories:\n")
print(table(comp_df$category))

# ── 5. Associate with clinical severity ───────────────────────────────────────

cat("\nComputing severity associations...\n")

t1_who    <- t1_meta$who_score_num
t2_who    <- t2_meta$who_score_num
who_delta <- t2_who - t1_who      # negative = improvement
names(t1_who) <- names(who_delta) <- patient_ids

cor_t1 <- apply(scores_mat, 2, function(x)
  tryCatch(cor(x, t1_who, method = "spearman", use = "complete.obs"), error = function(e) NA))
cor_delta <- apply(scores_mat, 2, function(x)
  tryCatch(cor(x, who_delta, method = "spearman", use = "complete.obs"), error = function(e) NA))

pval_t1 <- sapply(seq_len(K_total), function(k)
  tryCatch(cor.test(scores_mat[,k], t1_who, method = "spearman")$p.value, error = function(e) NA))
pval_delta <- sapply(seq_len(K_total), function(k)
  tryCatch(cor.test(scores_mat[,k], who_delta, method = "spearman")$p.value, error = function(e) NA))

assoc_df <- data.frame(
  component  = comp_names,
  category   = comp_df$category,
  rho_t1     = round(cor_t1, 3),
  pval_t1    = round(pval_t1, 4),
  fdr_t1     = round(p.adjust(pval_t1, method = "BH"), 4),
  rho_delta  = round(cor_delta, 3),
  pval_delta = round(pval_delta, 4),
  fdr_delta  = round(p.adjust(pval_delta, method = "BH"), 4),
  stringsAsFactors = FALSE
)

cat("\nAll components — severity associations:\n")
print(assoc_df[order(abs(assoc_df$rho_t1), decreasing = TRUE), ])

write.csv(assoc_df, "results/divas_associations.csv", row.names = FALSE)
cat("\nSaved: results/divas_associations.csv\n")

# ── 6. Visualisations ─────────────────────────────────────────────────────────

cat("\nGenerating plots...\n")

cat_colors <- c(
  "Shared (T1+T2)"  = "#377EB8",
  "T1-individual"   = "#FF7F00",
  "T2-individual"   = "#4DAF4A"
)

# 6a. Diagnostic angle plot
pdf(file.path(fig_dir, "divas_diagnostic.pdf"), width = 10, height = 8)
DJIVEAngleDiagnosticJP(
  datablock = data_list,
  dataname  = names(data_list),
  outstruct = divas_res,
  randseed  = 123,
  titlestr  = "DIVAS: Proteomics T1 vs T2"
)
dev.off()
cat("  Saved divas_diagnostic.pdf\n")

# 6b. Score heatmap: components × patients
# scores_mat rownames = sample names (e.g. "COVID_1_T1"); align annotation to match
samp_names <- rownames(scores_mat)
ann_row <- data.frame(
  T1_severity = t1_who[match(samp_names, t1_meta$sample)],
  WHO_delta   = who_delta[match(samp_names, t1_meta$sample)],
  row.names   = samp_names
)
ann_col <- data.frame(
  Category  = comp_df$category,
  rho_T1    = assoc_df$rho_t1[match(comp_names, assoc_df$component)],
  rho_delta = assoc_df$rho_delta[match(comp_names, assoc_df$component)],
  row.names = comp_names
)

pdf(file.path(fig_dir, "divas_score_heatmap.pdf"),
    width = max(8, K_total * 0.6 + 3), height = 8)
pheatmap(
  t(scores_mat),
  annotation_col    = ann_row,
  annotation_row    = ann_col,
  annotation_colors = list(Category = cat_colors),
  cluster_rows      = FALSE,
  cluster_cols      = TRUE,
  show_colnames     = FALSE,
  fontsize          = 8,
  main = "DIVAS scores: proteomics T1 vs T2",
  color = colorRampPalette(c("#313695", "white", "#A50026"))(100)
)
dev.off()
cat("  Saved divas_score_heatmap.pdf\n")

# 6c. Severity association waterfall
p_assoc <- assoc_df |>
  mutate(component = factor(component, levels = component[order(abs(rho_t1))])) |>
  ggplot(aes(x = component, y = rho_t1, fill = category)) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = cat_colors) +
  coord_flip() +
  labs(
    title = "DIVAS component correlation with T1 WHO severity",
    x = NULL, y = "Spearman rho",
    fill = "Component type"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "divas_severity_association.pdf"),
       p_assoc, width = 8, height = max(4, K_total * 0.35 + 2))
cat("  Saved divas_severity_association.pdf\n")

# 6d. Scatter: top component by |rho_t1| vs T1 severity and WHO delta
top_comp <- assoc_df$component[which.max(abs(assoc_df$rho_t1))]
cat(sprintf("\nTop component by |rho_t1|: %s (rho=%.3f, category=%s)\n",
            top_comp, assoc_df$rho_t1[assoc_df$component == top_comp],
            assoc_df$category[assoc_df$component == top_comp]))

scatter_df <- data.frame(
  sample    = samp_names,
  score      = scores_mat[, top_comp],
  t1_who     = t1_who[match(samp_names, t1_meta$sample)],
  who_delta  = who_delta[match(samp_names, t1_meta$sample)]
)

p1 <- ggplot(scatter_df, aes(x = score, y = t1_who)) +
  geom_point(aes(colour = t1_who), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.7) +
  scale_colour_gradient(low = "#E0F3F8", high = "#023858") +
  labs(title = paste0(top_comp, " vs T1 severity"),
       x = "DIVAS score", y = "WHO score (T1)", colour = "T1 WHO") +
  theme_bw(base_size = 10)

p2 <- ggplot(scatter_df, aes(x = score, y = who_delta)) +
  geom_point(aes(colour = who_delta), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.7) +
  scale_colour_gradient2(low = "steelblue", mid = "gray90", high = "firebrick", midpoint = 0) +
  labs(title = paste0(top_comp, " vs WHO change"),
       x = "DIVAS score", y = "WHO delta (T2−T1)", colour = "WHO\ndelta") +
  theme_bw(base_size = 10)

ggsave(file.path(fig_dir, "divas_top_scatter.pdf"),
       p1 + p2, width = 10, height = 4.5)
cat("  Saved divas_top_scatter.pdf\n")

# 6e. All components vs T1 severity — faceted scatter
all_scores <- as.data.frame(scores_mat)
all_scores$sample  <- samp_names
all_scores$t1_who  <- t1_who[match(samp_names, t1_meta$sample)]

all_long <- all_scores |>
  pivot_longer(all_of(comp_names), names_to = "component", values_to = "score") |>
  left_join(comp_df, by = "component")

p_all <- ggplot(all_long, aes(x = score, y = t1_who, colour = category)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.5) +
  facet_wrap(~component, scales = "free_x", ncol = 4) +
  scale_colour_manual(values = cat_colors) +
  labs(
    title = "All DIVAS components vs T1 WHO severity",
    x = "DIVAS score", y = "WHO score (T1)", colour = "Type"
  ) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "divas_all_components_vs_severity.pdf"),
       p_all, width = 12, height = max(4, ceiling(K_total / 4) * 3))
cat("  Saved divas_all_components_vs_severity.pdf\n")

# ── 7. Top proteins per component ─────────────────────────────────────────────

# For shared components, modality = "T1" and "T2"; for individual, modality = owning block.
# DIVAS scores name individual components "T1-Individual-k" but internally stores them
# as "T1-k" / "T2-k" in the loadings. Map accordingly.
get_feature_queries <- function(comp_name) {
  if (grepl("^T1\\+T2-", comp_name)) {
    # Shared: query both T1 and T2 modalities with the score name directly
    list(list(comp = comp_name, mod = "T1"),
         list(comp = comp_name, mod = "T2"))
  } else if (grepl("^T1-Individual-", comp_name)) {
    # T1-Individual-k → loadings store as "T1-k"
    k <- sub("^T1-Individual-", "", comp_name)
    list(list(comp = paste0("T1-", k), mod = "T1"))
  } else if (grepl("^T2-Individual-", comp_name)) {
    k <- sub("^T2-Individual-", "", comp_name)
    list(list(comp = paste0("T2-", k), mod = "T2"))
  } else {
    list()
  }
}

cat("\nTop proteins per component:\n")
for (comp in comp_names) {
  queries <- get_feature_queries(comp)
  for (q in queries) {
    tryCatch({
      feats <- getTopFeatures(divas_res, compName = q$comp, modName = q$mod,
                              n_top_pos = 10, n_top_neg = 10)
      cat(sprintf("\n  %s | %s\n", comp, q$mod))
      cat("    + :", paste(head(feats$top_positive, 5), collapse = "; "), "\n")
      cat("    − :", paste(head(feats$top_negative, 5), collapse = "; "), "\n")
    }, error = function(e) {
      cat(sprintf("  [%s | %s: %s]\n", comp, q$mod, conditionMessage(e)))
    })
  }
}

# ── 8. Summary ────────────────────────────────────────────────────────────────

cat("\n=== Summary ===\n")
cat(sprintf("DIVAS decomposed 481 proteins × 120 patients into %d components:\n", K_total))
for (cat_name in unique(comp_df$category)) {
  n <- sum(comp_df$category == cat_name)
  if (n > 0) {
    cat(sprintf("  %s: %d component(s)\n", cat_name, n))
  }
}
cat("\nKey outputs:\n")
cat("  results/divas_2block.rds         — DIVAS result object\n")
cat("  results/divas_associations.csv   — severity associations\n")
cat("  results/figs/divas/              — all figures\n")
