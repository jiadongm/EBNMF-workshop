# prepareData.R — COVID-19 Proteomics Data Loading and Preparation
#
# Data: Su et al. (Cell 2020) COVID-19 multiomics cohort.
#   File: data/240_input_4omics/proteomics_120patients.csv
#   Dimensions: 481 proteins × 240 samples (120 patients × 2 time points T1/T2)
#   Format: rows = proteins, columns = samples (COVID_[id]_T[1|2])
#
# Normalisation:
#   The file contains Olink NPX (Normalized Protein eXpression) values, which are
#   already in log2 scale and batch-corrected. Olink's processing pipeline computes
#   NPX from raw Ct values by normalising against an extension control, an inter-plate
#   control, and a batch correction factor derived from overlapping pooled reference
#   samples (Su et al. STAR Methods, "Plasma proteomics"). This is the primary
#   normalised form used throughout the paper for all downstream analyses.
#
#   We therefore load the NPX values as-is: no further transformation is applied.
#
#   Note on negatives: NPX can be negative (protein detected below the LOD or at low
#   abundance). For signed EBMF (point-normal or point-Laplace prior), NPX values are
#   used directly. For NMF (point-exponential prior, requires non-negative input),
#   a per-protein minimum-shift will be applied in the analysis scripts — that step
#   is kept separate to preserve the original NPX values here.

library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)

# ── Paths ─────────────────────────────────────────────────────────────────────

data_dir    <- file.path("data", "240_input_4omics")
prot_file   <- file.path(data_dir, "proteomics_120patients.csv")
meta_file   <- file.path("data", "patientInfo.xlsx")
results_dir <- "results"
fig_dir     <- file.path(results_dir, "figs")

dir.create(results_dir, showWarnings = FALSE)
dir.create(fig_dir,     showWarnings = FALSE)

# ── 1. Load raw proteomics matrix ─────────────────────────────────────────────

cat("Loading proteomics data...\n")

prot_raw <- read.csv(prot_file, row.names = 1, check.names = FALSE)

cat(sprintf("  Dimensions: %d proteins × %d samples\n",
            nrow(prot_raw), ncol(prot_raw)))
cat(sprintf("  Value range: [%.3f, %.3f]\n", min(prot_raw), max(prot_raw)))
cat(sprintf("  Missing values: %d\n", sum(is.na(prot_raw))))

# ── 2. Parse sample metadata from column names ─────────────────────────────────
# Column names follow the pattern: COVID_[patient_id]_T[1|2]

sample_names <- colnames(prot_raw)
parts        <- strsplit(sample_names, "_")

meta <- data.frame(
  sample     = sample_names,
  patient_id = as.integer(sapply(parts, `[`, 2)),
  time_point = sapply(parts, `[`, 3),
  stringsAsFactors = FALSE
)

cat(sprintf("\nSample IDs parsed: %d samples, %d patients, time points: %s\n",
            nrow(meta), length(unique(meta$patient_id)),
            paste(unique(meta$time_point), collapse = ", ")))

# ── 3. Load and join patient-level clinical metadata ───────────────────────────
# Source: Su et al. Table S1.1 (patientInfo.xlsx, sheet 'S1.1 Patient Clinical Data')
#
# ID mapping: INCOV036-1 → COVID_36_T1
#   Study Subject ID: INCOV[zero-padded patient id]  → extract numeric patient id
#   Sample ID suffix: -1 = T1, -2 = T2

cat("Loading patient clinical metadata from patientInfo.xlsx...\n")

clinical_raw <- read_excel(meta_file, sheet = "S1.1 Patient Clinical Data")

# Parse patient_id and time_point from the original sample IDs
clinical_raw$patient_id <- as.integer(
  sub("INCOV0*([0-9]+)", "\\1", clinical_raw$`Study Subject ID`)
)
clinical_raw$time_point <- ifelse(
  sub(".*-([12])$", "\\1", clinical_raw$`Sample ID`) == "1", "T1", "T2"
)

# Select and rename columns of interest
clinical <- clinical_raw |>
  select(
    patient_id,
    time_point,
    who_score       = `Who Ordinal Scale`,
    location        = `Patient Location`,
    sex             = Sex,
    age             = Age,
    bmi             = BMI,
    mech_vent       = `Mechanical Ventilation`,
    resp_support    = `Respiratory Support`,
    has_proteomics  = `Proteomics?`
  ) |>
  mutate(
    # WHO score: coerce to numeric; "1 or 2" → 1.5 (midpoint, flagged separately)
    who_score_char    = who_score,
    who_score_num     = suppressWarnings(as.numeric(who_score)),
    who_ambiguous     = who_score == "1 or 2",
    who_score_num     = ifelse(who_ambiguous, 1.5, who_score_num),
    age               = as.numeric(age),
    bmi               = as.numeric(bmi)
  )

# Join to sample-level metadata (left join: keep all 240 of our samples)
meta <- meta |>
  left_join(clinical, by = c("patient_id", "time_point"))

n_matched <- sum(!is.na(meta$who_score))
cat(sprintf("  Clinical metadata joined: %d / %d samples matched\n",
            n_matched, nrow(meta)))
cat(sprintf("  WHO score distribution (numeric):\n"))
print(table(meta$who_score_char, useNA = "always"))
cat(sprintf("  Samples with 'Proteomics?' = yes: %d\n",
            sum(meta$has_proteomics == "yes", na.rm = TRUE)))

cat(sprintf("\nSample metadata:\n"))
cat(sprintf("  Total samples:   %d\n", nrow(meta)))
cat(sprintf("  Unique patients: %d\n", length(unique(meta$patient_id))))
cat(sprintf("  T1 samples: %d,  T2 samples: %d\n",
            sum(meta$time_point == "T1"), sum(meta$time_point == "T2")))

# ── 4. Basic QC ────────────────────────────────────────────────────────────────

# Per-protein summary
prot_summary <- data.frame(
  protein  = rownames(prot_raw),
  mean     = rowMeans(prot_raw),
  sd       = apply(prot_raw, 1, sd),
  min      = apply(prot_raw, 1, min),
  max      = apply(prot_raw, 1, max),
  n_neg    = rowSums(prot_raw < 0)   # number of samples with NPX < 0 (below LOD)
)

cat(sprintf("\nPer-protein QC:\n"))
cat(sprintf("  Proteins with any NPX < 0:  %d / %d\n",
            sum(prot_summary$n_neg > 0), nrow(prot_summary)))
cat(sprintf("  Proteins with all NPX < 0:  %d\n",
            sum(prot_summary$min < 0 & prot_summary$max < 0)))
cat(sprintf("  Proteins with SD == 0:      %d\n",
            sum(prot_summary$sd == 0)))

# Per-sample library: check for outlier samples
sample_means <- colMeans(prot_raw)
cat(sprintf("\nPer-sample mean NPX: %.2f ± %.2f  [%.2f, %.2f]\n",
            mean(sample_means), sd(sample_means),
            min(sample_means), max(sample_means)))

# ── 5. QC plots ───────────────────────────────────────────────────────────────

# 4a. Distribution of NPX values (violin per time point)
prot_long <- prot_raw |>
  as.data.frame() |>
  tibble::rownames_to_column("protein") |>
  pivot_longer(-protein, names_to = "sample", values_to = "NPX") |>
  left_join(meta, by = "sample")

p_dist <- ggplot(prot_long, aes(x = time_point, y = NPX, fill = time_point)) +
  geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(title = "Proteomics: NPX distribution by time point",
       subtitle = "Olink NPX (log2 scale, already normalised)",
       x = "Time point", y = "NPX (log2)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "prot_npx_distribution.pdf"),
       p_dist, width = 5, height = 4)

# 4b. Per-sample mean NPX (check for outlier samples)
meta$mean_npx <- sample_means[meta$sample]

p_sample <- ggplot(meta, aes(x = reorder(sample, mean_npx), y = mean_npx,
                              colour = time_point)) +
  geom_point(size = 1) +
  labs(title = "Per-sample mean NPX",
       subtitle = "Sorted by mean; colour = time point",
       x = "Sample (sorted)", y = "Mean NPX (log2)",
       colour = "Time point") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(fig_dir, "prot_sample_means.pdf"),
       p_sample, width = 7, height = 4)

# 4c. Per-protein mean and SD
p_prot <- ggplot(prot_summary, aes(x = mean, y = sd)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "firebrick") +
  labs(title = "Proteomics: mean vs SD per protein",
       subtitle = "Dashed line: SD = 0.5 (low-variability threshold)",
       x = "Mean NPX (log2)", y = "SD across 240 samples") +
  theme_bw()

ggsave(file.path(fig_dir, "prot_mean_sd.pdf"),
       p_prot, width = 5, height = 4)

cat("\nQC figures saved to", fig_dir, "\n")

# ── 6. Prepare analysis matrix ─────────────────────────────────────────────────
# Convention for EBMF/flashier: X is samples × proteins (n × p).
# The input CSV is proteins × samples; we transpose here.

# Filter: remove proteins with zero variance (uninformative)
low_var <- prot_summary$protein[prot_summary$sd == 0]
if (length(low_var) > 0) {
  cat(sprintf("\nRemoving %d zero-variance proteins\n", length(low_var)))
  prot_raw <- prot_raw[!(rownames(prot_raw) %in% low_var), ]
}

# Final matrix: samples × proteins (n × p), NPX log2 values
X_prot <- t(as.matrix(prot_raw))   # 240 × 481 (samples × proteins)

cat(sprintf("\nAnalysis matrix X_prot: %d samples × %d proteins\n",
            nrow(X_prot), ncol(X_prot)))
cat(sprintf("  Values in log2 NPX scale: [%.3f, %.3f]\n",
            min(X_prot), max(X_prot)))
cat("  (Negative values reflect proteins near or below the Olink LOD.)\n")

# ── 7. Save ───────────────────────────────────────────────────────────────────

saveRDS(X_prot, file.path(results_dir, "X_prot.rds"))
saveRDS(meta,   file.path(results_dir, "sample_meta.rds"))
write.csv(prot_summary, file.path(results_dir, "prot_qc_summary.csv"),
          row.names = FALSE)

cat("\nSaved:\n")
cat("  results/X_prot.rds          — samples × proteins matrix (NPX log2)\n")
cat("  results/sample_meta.rds     — sample metadata with clinical variables\n")
cat("  results/prot_qc_summary.csv — per-protein QC statistics\n")
cat("\nClinical columns in sample_meta:\n")
cat("  patient_id, time_point, who_score_char, who_score_num, who_ambiguous,\n")
cat("  location, sex, age, bmi, mech_vent, resp_support\n")
