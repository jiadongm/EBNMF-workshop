# 08c_divas_individual_loadings.R — Bar plots of top 20 proteins per individual DIVAS component
#
# DIVAS identified 13 individual (time-specific) components:
#   7 T1-individual and 6 T2-individual.
# These are interpreted via their protein loadings — which proteins drive each axis.

library(ggplot2)
library(patchwork)

setwd("/data/gpfs/projects/punim0613/JiaDong/NatashaSharma/ns_snail_spatial/EBNMF/CovidCase")

fig_dir <- "results/figs/divas/others"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load DIVAS results ──────────────────────────────────────────────────────

divas_res <- readRDS("results/divas_2block.rds")

# Individual component loadings are stored internally as "T1-k" / "T2-k"
# in divas_res$Loadings[["T1"]] and divas_res$Loadings[["T2"]]
L_T1 <- divas_res$Loadings[["T1"]]  # 481 × 18 (11 shared + 7 individual)
L_T2 <- divas_res$Loadings[["T2"]]  # 481 × 17 (11 shared + 6 individual)

cat(sprintf("T1 loadings: %d proteins × %d columns\n", nrow(L_T1), ncol(L_T1)))
cat(sprintf("T2 loadings: %d proteins × %d columns\n", nrow(L_T2), ncol(L_T2)))

# Identify individual component columns (not shared "T1+T2-k")
t1_indiv_cols <- grep("^T1-\\d+$", colnames(L_T1), value = TRUE)
t2_indiv_cols <- grep("^T2-\\d+$", colnames(L_T2), value = TRUE)

cat(sprintf("T1 individual components: %d (%s)\n",
            length(t1_indiv_cols), paste(t1_indiv_cols, collapse = ", ")))
cat(sprintf("T2 individual components: %d (%s)\n",
            length(t2_indiv_cols), paste(t2_indiv_cols, collapse = ", ")))

# ── 2. Helper: bar plot of top 20 proteins ──────────────────────────────────────

truncate_label <- function(x, max_len = 50) {
  out <- ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len - 3), "..."), x)
  make.unique(out, sep = " ")
}

plot_top20 <- function(loading_vec, comp_name) {
  # Rank by absolute value, take top 20
  ord <- order(abs(loading_vec), decreasing = TRUE)
  top_idx <- head(ord, 20)

  df <- data.frame(
    protein  = truncate_label(names(loading_vec)[top_idx]),
    loading  = loading_vec[top_idx],
    sign     = ifelse(loading_vec[top_idx] >= 0, "Positive", "Negative"),
    stringsAsFactors = FALSE
  )
  # Order by loading value for visual clarity
  df$protein <- factor(df$protein, levels = df$protein[order(df$loading)])

  ggplot(df, aes(x = loading, y = protein, fill = sign)) +
    geom_col() +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("Positive" = "#D73027", "Negative" = "#4575B4")) +
    labs(
      title = paste0(comp_name, " - Top 20 proteins by |loading|"),
      x = "Loading", y = NULL, fill = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
}

# ── 3. Generate individual PDFs ─────────────────────────────────────────────────

all_plots <- list()

# T1-individual components
for (col in t1_indiv_cols) {
  vec <- L_T1[, col]
  names(vec) <- rownames(L_T1)
  display_name <- sub("^T1-", "T1-Individual-", col)

  p <- plot_top20(vec, display_name)
  all_plots[[display_name]] <- p

  fname <- file.path(fig_dir, sprintf("indiv_%s_top20.pdf", col))
  ggsave(fname, p, width = 7, height = 6)
  cat(sprintf("  Saved %s\n", basename(fname)))
}

# T2-individual components
for (col in t2_indiv_cols) {
  vec <- L_T2[, col]
  names(vec) <- rownames(L_T2)
  display_name <- sub("^T2-", "T2-Individual-", col)

  p <- plot_top20(vec, display_name)
  all_plots[[display_name]] <- p

  fname <- file.path(fig_dir, sprintf("indiv_%s_top20.pdf", col))
  ggsave(fname, p, width = 7, height = 6)
  cat(sprintf("  Saved %s\n", basename(fname)))
}

# ── 4. Combined PDF with all 13 panels ──────────────────────────────────────────

n_panels <- length(all_plots)
ncol_layout <- 3
nrow_layout <- ceiling(n_panels / ncol_layout)

combined <- wrap_plots(all_plots, ncol = ncol_layout) +
  plot_annotation(title = "DIVAS Individual Components — Top 20 Proteins by |Loading|")

ggsave(
  file.path(fig_dir, "indiv_all_top20.pdf"),
  combined,
  width = 21, height = nrow_layout * 6
)
cat(sprintf("\n  Saved indiv_all_top20.pdf (%d panels, %d × %d in)\n",
            n_panels, 21, nrow_layout * 6))

cat("\nDone.\n")
