# 08b_divas_metadata.R — DIVAS scores vs clinical metadata
#
# Associations with: T1 severity, T2 severity, WHO delta, age, BMI, sex,
#                    location (ordinal), mechanical ventilation, resp support

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(pheatmap)

setwd("/data/gpfs/projects/punim0613/JiaDong/NatashaSharma/ns_snail_spatial/EBNMF/CovidCase")

fig_dir <- "results/figs/divas/others"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────

divas_res  <- readRDS("results/divas_2block.rds")
meta       <- readRDS("results/sample_meta.rds")
scores_mat <- divas_res$Scores
comp_names <- colnames(scores_mat)
K <- ncol(scores_mat)

t1_meta <- meta[meta$time_point == "T1", ] |> arrange(patient_id)
t2_meta <- meta[meta$time_point == "T2", ] |> arrange(patient_id)
patient_ids <- as.character(t1_meta$patient_id)

# Categorise components
comp_cat <- sapply(comp_names, function(nm) {
  if (grepl("^T1\\+T2-", nm))  "Shared (T1+T2)"
  else if (grepl("^T1-", nm))  "T1-individual"
  else                          "T2-individual"
})

# Build per-patient clinical data frame
clin <- data.frame(
  patient_id   = patient_ids,
  t1_who       = t1_meta$who_score_num,
  t2_who       = t2_meta$who_score_num,
  who_delta    = t2_meta$who_score_num - t1_meta$who_score_num,
  age          = t1_meta$age,
  bmi          = t1_meta$bmi,
  sex          = t1_meta$sex,
  location_t1  = t1_meta$location,
  location_t2  = t2_meta$location,
  mech_vent_t1 = t1_meta$mech_vent,
  mech_vent_t2 = t2_meta$mech_vent,
  resp_t1      = t1_meta$resp_support,
  resp_t2      = t2_meta$resp_support,
  stringsAsFactors = FALSE
)

# Ordinal encoding: location
loc_ord <- c("Home (mobile phlebotomy)" = 1, "Clinic" = 2, "Hospital" = 3, "ICU" = 4)
clin$loc_t1_ord <- loc_ord[clin$location_t1]
clin$loc_t2_ord <- loc_ord[clin$location_t2]
clin$loc_delta  <- clin$loc_t2_ord - clin$loc_t1_ord

# Ordinal encoding: respiratory support
resp_ord <- c("None" = 0, "Nasal cannula" = 1,
              "High flow nasal cannula (HFNC)" = 2, "Other" = 3)
clin$resp_t1_ord <- resp_ord[clin$resp_t1]
clin$resp_t2_ord <- resp_ord[clin$resp_t2]

cat("Clinical variables ready. N =", nrow(clin), "patients\n")
cat("BMI non-NA:", sum(!is.na(clin$bmi)), "/ 120\n\n")

# ── 2. Comprehensive correlation matrix ──────────────────────────────────────

# Continuous/ordinal variables to test
cont_vars <- c("t1_who", "t2_who", "who_delta", "age", "bmi",
               "loc_t1_ord", "loc_t2_ord", "loc_delta",
               "resp_t1_ord", "resp_t2_ord")
cont_labels <- c("T1 severity", "T2 severity", "WHO delta", "Age", "BMI",
                 "Location T1", "Location T2", "Location change",
                 "Resp support T1", "Resp support T2")

# Spearman rho matrix: components × clinical variables
rho_mat  <- matrix(NA, K, length(cont_vars),
                   dimnames = list(comp_names, cont_labels))
pval_mat <- rho_mat

for (j in seq_along(cont_vars)) {
  y <- clin[[cont_vars[j]]]
  for (i in seq_len(K)) {
    tt <- tryCatch(
      cor.test(scores_mat[, i], y, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    if (!is.null(tt)) {
      rho_mat[i, j]  <- tt$estimate
      pval_mat[i, j] <- tt$p.value
    }
  }
}

# BH-adjusted p-values (across all tests)
fdr_vec <- p.adjust(as.vector(pval_mat), method = "BH")
fdr_mat <- matrix(fdr_vec, K, length(cont_vars),
                  dimnames = dimnames(rho_mat))

# Print top hits
cat("=== Top associations (|rho| > 0.2 and FDR < 0.1) ===\n")
for (j in seq_len(ncol(rho_mat))) {
  for (i in seq_len(K)) {
    if (!is.na(fdr_mat[i,j]) && abs(rho_mat[i,j]) > 0.2 && fdr_mat[i,j] < 0.1) {
      cat(sprintf("  %s × %s: rho=%.3f, FDR=%.4f\n",
                  rownames(rho_mat)[i], colnames(rho_mat)[j],
                  rho_mat[i,j], fdr_mat[i,j]))
    }
  }
}

# ── 3. Correlation heatmap ────────────────────────────────────────────────────

# Annotation for component rows
ann_row <- data.frame(Category = comp_cat, row.names = comp_names)
cat_colors <- list(Category = c(
  "Shared (T1+T2)" = "#377EB8",
  "T1-individual"  = "#FF7F00",
  "T2-individual"  = "#4DAF4A"
))

# Stars for significance
display_mat <- rho_mat
star_mat <- matrix("", K, length(cont_vars), dimnames = dimnames(rho_mat))
star_mat[fdr_mat < 0.001] <- "***"
star_mat[fdr_mat >= 0.001 & fdr_mat < 0.01] <- "**"
star_mat[fdr_mat >= 0.01 & fdr_mat < 0.05] <- "*"

pdf(file.path(fig_dir, "correlation_heatmap.pdf"), width = 10, height = 10)
pheatmap(
  rho_mat,
  display_numbers    = star_mat,
  fontsize_number    = 10,
  annotation_row     = ann_row,
  annotation_colors  = cat_colors,
  cluster_rows       = FALSE,
  cluster_cols       = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-0.8, 0.8, length.out = 101),
  fontsize = 8,
  main = "DIVAS scores vs clinical variables (Spearman rho)\n* FDR<0.05  ** FDR<0.01  *** FDR<0.001",
  angle_col = 45
)
dev.off()
cat("\nSaved correlation_heatmap.pdf\n")

# Save full table
assoc_full <- expand.grid(component = comp_names, variable = cont_labels,
                          stringsAsFactors = FALSE)
assoc_full$category <- comp_cat[assoc_full$component]
assoc_full$rho <- as.vector(rho_mat)
assoc_full$pval <- as.vector(pval_mat)
assoc_full$fdr <- as.vector(fdr_mat)
write.csv(assoc_full, "results/divas_metadata_associations.csv", row.names = FALSE)
cat("Saved divas_metadata_associations.csv\n")

# ── 4. T2 severity scatter plots ─────────────────────────────────────────────

# Top 6 components by |rho| with T2 severity
t2_rho <- rho_mat[, "T2 severity"]
top6_t2 <- names(sort(abs(t2_rho), decreasing = TRUE))[1:6]

scatter_df <- as.data.frame(scores_mat[, top6_t2])
scatter_df$t2_who <- clin$t2_who

p_t2 <- scatter_df |>
  pivot_longer(all_of(top6_t2), names_to = "component", values_to = "score") |>
  mutate(component = factor(component, levels = top6_t2)) |>
  ggplot(aes(x = score, y = t2_who)) +
  geom_point(aes(colour = t2_who), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.5) +
  facet_wrap(~component, scales = "free_x", ncol = 3) +
  scale_colour_gradient(low = "#E0F3F8", high = "#023858") +
  labs(title = "Top DIVAS components vs T2 severity",
       subtitle = paste(sapply(top6_t2, function(x)
         sprintf("%s: rho=%.2f", x, t2_rho[x])), collapse = "  |  "),
       x = "DIVAS score", y = "WHO score (T2)", colour = "T2 WHO") +
  theme_bw(base_size = 9) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "t2_severity_scatter.pdf"), p_t2, width = 11, height = 7)
cat("Saved t2_severity_scatter.pdf\n")

# ── 5. Age association ────────────────────────────────────────────────────────

age_rho <- rho_mat[, "Age"]
top6_age <- names(sort(abs(age_rho), decreasing = TRUE))[1:6]

scatter_age <- as.data.frame(scores_mat[, top6_age])
scatter_age$age <- clin$age

p_age <- scatter_age |>
  pivot_longer(all_of(top6_age), names_to = "component", values_to = "score") |>
  mutate(component = factor(component, levels = top6_age)) |>
  ggplot(aes(x = score, y = age)) +
  geom_point(aes(colour = age), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.5) +
  facet_wrap(~component, scales = "free_x", ncol = 3) +
  scale_colour_viridis_c(option = "C") +
  labs(title = "Top DIVAS components vs Age",
       subtitle = paste(sapply(top6_age, function(x)
         sprintf("%s: rho=%.2f", x, age_rho[x])), collapse = "  |  "),
       x = "DIVAS score", y = "Age (years)", colour = "Age") +
  theme_bw(base_size = 9) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "age_scatter.pdf"), p_age, width = 11, height = 7)
cat("Saved age_scatter.pdf\n")

# ── 6. BMI association (n=77 non-NA) ──────────────────────────────────────────

bmi_rho <- rho_mat[, "BMI"]
top6_bmi <- names(sort(abs(bmi_rho), decreasing = TRUE))[1:6]

scatter_bmi <- as.data.frame(scores_mat[, top6_bmi])
scatter_bmi$bmi <- clin$bmi

p_bmi <- scatter_bmi |>
  filter(!is.na(bmi)) |>
  pivot_longer(all_of(top6_bmi), names_to = "component", values_to = "score") |>
  mutate(component = factor(component, levels = top6_bmi)) |>
  ggplot(aes(x = score, y = bmi)) +
  geom_point(aes(colour = bmi), size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.5) +
  facet_wrap(~component, scales = "free_x", ncol = 3) +
  scale_colour_viridis_c(option = "D") +
  labs(title = sprintf("Top DIVAS components vs BMI (n=%d with BMI data)",
                       sum(!is.na(clin$bmi))),
       subtitle = paste(sapply(top6_bmi, function(x)
         sprintf("%s: rho=%.2f", x, bmi_rho[x])), collapse = "  |  "),
       x = "DIVAS score", y = "BMI", colour = "BMI") +
  theme_bw(base_size = 9) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "bmi_scatter.pdf"), p_bmi, width = 11, height = 7)
cat("Saved bmi_scatter.pdf\n")

# ── 7. Sex differences (box plots) ───────────────────────────────────────────

sex_df <- as.data.frame(scores_mat)
sex_df$sex <- clin$sex

# Wilcoxon test per component
sex_pvals <- sapply(comp_names, function(cn) {
  tryCatch(wilcox.test(scores_mat[clin$sex == "Male", cn],
                       scores_mat[clin$sex == "Female", cn])$p.value,
           error = function(e) NA)
})
sex_fdr <- p.adjust(sex_pvals, method = "BH")

# Top 6 by significance
top6_sex <- names(sort(sex_pvals))[1:6]

p_sex <- sex_df |>
  pivot_longer(all_of(top6_sex), names_to = "component", values_to = "score") |>
  mutate(component = factor(component, levels = top6_sex)) |>
  ggplot(aes(x = sex, y = score, fill = sex)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.4) +
  facet_wrap(~component, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Female" = "#E78AC3", "Male" = "#66C2A5")) +
  labs(title = "Top DIVAS components by Sex (Wilcoxon test)",
       subtitle = paste(sapply(top6_sex, function(x)
         sprintf("%s: p=%.3f, FDR=%.3f", x, sex_pvals[x], sex_fdr[x])),
         collapse = "  |  "),
       x = NULL, y = "DIVAS score") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "sex_boxplots.pdf"), p_sex, width = 11, height = 7)
cat("Saved sex_boxplots.pdf\n")

# ── 8. Location (ordinal) box plots ──────────────────────────────────────────

loc_rho_t1 <- rho_mat[, "Location T1"]
top6_loc <- names(sort(abs(loc_rho_t1), decreasing = TRUE))[1:6]

loc_df <- as.data.frame(scores_mat[, top6_loc])
loc_df$location <- factor(clin$location_t1,
                          levels = c("Home (mobile phlebotomy)", "Clinic",
                                     "Hospital", "ICU"))

p_loc <- loc_df |>
  pivot_longer(all_of(top6_loc), names_to = "component", values_to = "score") |>
  mutate(component = factor(component, levels = top6_loc)) |>
  ggplot(aes(x = location, y = score, fill = location)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.4) +
  facet_wrap(~component, scales = "free_y", ncol = 3) +
  scale_fill_brewer(palette = "YlOrRd") +
  labs(title = "Top DIVAS components by T1 Location",
       subtitle = paste(sapply(top6_loc, function(x)
         sprintf("%s: rho=%.2f", x, loc_rho_t1[x])), collapse = "  |  "),
       x = NULL, y = "DIVAS score") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(fig_dir, "location_boxplots.pdf"), p_loc, width = 11, height = 7)
cat("Saved location_boxplots.pdf\n")

# ── 9. Mechanical ventilation (box plots) ────────────────────────────────────

vent_pvals <- sapply(comp_names, function(cn) {
  tryCatch(wilcox.test(scores_mat[clin$mech_vent_t1, cn],
                       scores_mat[!clin$mech_vent_t1, cn])$p.value,
           error = function(e) NA)
})
vent_fdr <- p.adjust(vent_pvals, method = "BH")
top6_vent <- names(sort(vent_pvals))[1:6]

vent_df <- as.data.frame(scores_mat[, top6_vent])
vent_df$mech_vent <- ifelse(clin$mech_vent_t1, "Ventilated", "Not ventilated")

p_vent <- vent_df |>
  pivot_longer(all_of(top6_vent), names_to = "component", values_to = "score") |>
  mutate(component = factor(component, levels = top6_vent)) |>
  ggplot(aes(x = mech_vent, y = score, fill = mech_vent)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.4) +
  facet_wrap(~component, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Not ventilated" = "#A6D854", "Ventilated" = "#FC8D62")) +
  labs(title = "Top DIVAS components by Mechanical Ventilation (T1)",
       subtitle = paste(sapply(top6_vent, function(x)
         sprintf("%s: p=%.3f, FDR=%.3f", x, vent_pvals[x], vent_fdr[x])),
         collapse = "  |  "),
       x = NULL, y = "DIVAS score") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "ventilation_boxplots.pdf"), p_vent, width = 11, height = 7)
cat("Saved ventilation_boxplots.pdf\n")

# ── 10. T1 vs T2 severity comparison ─────────────────────────────────────────

# Scatter: T1 severity rho vs T2 severity rho per component
sev_comp <- data.frame(
  component = comp_names,
  category  = comp_cat,
  rho_T1    = rho_mat[, "T1 severity"],
  rho_T2    = rho_mat[, "T2 severity"],
  fdr_T1    = fdr_mat[, "T1 severity"],
  fdr_T2    = fdr_mat[, "T2 severity"]
)
sev_comp$sig <- case_when(
  sev_comp$fdr_T1 < 0.05 & sev_comp$fdr_T2 < 0.05 ~ "Both",
  sev_comp$fdr_T1 < 0.05                           ~ "T1 only",
  sev_comp$fdr_T2 < 0.05                           ~ "T2 only",
  TRUE                                              ~ "Neither"
)

p_comp <- ggplot(sev_comp, aes(x = rho_T1, y = rho_T2,
                                colour = category, shape = sig)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = sev_comp |> filter(abs(rho_T1) > 0.15 | abs(rho_T2) > 0.15),
    aes(label = component), size = 2.5, max.overlaps = 15, show.legend = FALSE
  ) +
  scale_colour_manual(values = c("Shared (T1+T2)" = "#377EB8",
                                  "T1-individual" = "#FF7F00",
                                  "T2-individual" = "#4DAF4A")) +
  scale_shape_manual(values = c("Both" = 16, "T1 only" = 17,
                                 "T2 only" = 15, "Neither" = 1)) +
  coord_equal(xlim = c(-0.8, 0.8), ylim = c(-0.8, 0.8)) +
  labs(title = "DIVAS components: T1 vs T2 severity association",
       subtitle = "Components near diagonal = stable severity axis; off-diagonal = time-specific",
       x = "Spearman rho with T1 severity",
       y = "Spearman rho with T2 severity",
       colour = "Component type", shape = "FDR < 0.05") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "t1_vs_t2_severity.pdf"), p_comp, width = 8, height = 7)
cat("Saved t1_vs_t2_severity.pdf\n")

# ── 11. Multi-variable bubble plot ────────────────────────────────────────────

# Show all components × all variables as a bubble plot (size = |rho|, colour = rho)
bubble_df <- assoc_full |>
  mutate(abs_rho = abs(rho),
         sig = ifelse(fdr < 0.05, "FDR<0.05", "n.s."),
         component = factor(component, levels = comp_names))

p_bubble <- ggplot(bubble_df,
                   aes(x = variable, y = component, size = abs_rho, colour = rho)) +
  geom_point(aes(shape = sig), alpha = 0.8) +
  scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-0.8, 0.8)) +
  scale_size_continuous(range = c(0.5, 6), limits = c(0, 0.8)) +
  scale_shape_manual(values = c("FDR<0.05" = 16, "n.s." = 1)) +
  labs(title = "DIVAS scores × Clinical variables",
       x = NULL, y = NULL,
       colour = "Spearman\nrho", size = "|rho|", shape = "Significance") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(colour = "grey95"))

ggsave(file.path(fig_dir, "bubble_all_associations.pdf"),
       p_bubble, width = 10, height = 10)
cat("Saved bubble_all_associations.pdf\n")

# ── 12. Respiratory support (box plots) ──────────────────────────────────────

resp_rho_t1 <- rho_mat[, "Resp support T1"]
top6_resp <- names(sort(abs(resp_rho_t1), decreasing = TRUE))[1:6]

resp_df <- as.data.frame(scores_mat[, top6_resp])
resp_df$resp <- factor(clin$resp_t1,
                       levels = c("None", "Nasal cannula",
                                  "High flow nasal cannula (HFNC)", "Other"))

p_resp <- resp_df |>
  pivot_longer(all_of(top6_resp), names_to = "component", values_to = "score") |>
  mutate(component = factor(component, levels = top6_resp)) |>
  ggplot(aes(x = resp, y = score, fill = resp)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.4) +
  facet_wrap(~component, scales = "free_y", ncol = 3) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  labs(title = "Top DIVAS components by Respiratory Support (T1)",
       subtitle = paste(sapply(top6_resp, function(x)
         sprintf("%s: rho=%.2f", x, resp_rho_t1[x])), collapse = "  |  "),
       x = NULL, y = "DIVAS score") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(fig_dir, "resp_support_boxplots.pdf"), p_resp, width = 11, height = 7)
cat("Saved resp_support_boxplots.pdf\n")

# ── Done ──────────────────────────────────────────────────────────────────────

cat("\n=== Done ===\n")
cat("All plots saved to", fig_dir, "\n")
cat("Full association table: results/divas_metadata_associations.csv\n")
