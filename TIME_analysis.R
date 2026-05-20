###############################################################################
# CIBERSORTx Immune Deconvolution Workflow
#
# Purpose:
#   Clean, aggregate, visualize, and statistically compare CIBERSORTx immune
#   cell deconvolution results across biological groups.
#
# Main outputs:
#   - Relative immune composition stacked barplot
#   - Absolute immune score + relative composition combined plot
#   - Global PERMANOVA
#   - Pairwise PERMANOVA
#   - Pairwise betadisper dispersion control
#   - PCoA plot
#   - Pairwise PERMANOVA bubble plot
#
# Expected inputs:
#   - CIBERSORTx result file with LM22-like cell columns and P.value
#   - Optional second CIBERSORTx result file containing Absolute score
#   - sample metadata table with SampleID and MetaGroup columns
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(vegan)
})

# =============================================================================
# 0. USER CONFIGURATION
# =============================================================================

config <- list(
  cibersort_fraction_file = "CIBERSORTx_Job100_Results.csv",
  cibersort_absolute_file = "CIBERSORTx_Job102_Results.csv",
  sample_metadata_file    = "sample_meta_canonical.csv",
  outdir                  = "cibersort_outputs",
  pvalue_cutoff           = 0.05,
  permutations            = 9999
)

dir.create(config$outdir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. GROUP AND COLOR DEFINITIONS
# =============================================================================

group_mapping <- c(
  "Decidua"             = "Decidua",
  "Healthy_Endometrium" = "Healthy",
  "Healthy_EC"          = "HealthyEC",
  "EC"                  = "EC"
)

group_order <- c("Decidua", "Healthy", "HealthyEC", "EC")

publication_labels <- c(
  "Healthy_Endometrium" = "Healthy - Endometrium Healthy Tissue",
  "Healthy_EC"          = "Healthy - EC Adjacent Healthy Tissue",
  "Decidua"             = "Decidua - MFI",
  "EC"                  = "EC - TIME"
)

publication_group_order <- c(
  "Healthy - Endometrium Healthy Tissue",
  "Healthy - EC Adjacent Healthy Tissue",
  "Decidua - MFI",
  "EC - TIME"
)

group_colors_immune <- c(
  "Healthy - Endometrium Healthy Tissue" = "#1B6294",
  "Healthy - EC Adjacent Healthy Tissue" = "#72B9ED",
  "Decidua - MFI"                        = "#7B9E87",
  "EC - TIME"                            = "#7A0B12"
)

immune_palette <- colorRampPalette(
  c("#373D70", "#A04389", "#81B29A", "#3C6451", "#E9B157", "#E07A5F", "#EAB69F")
)(7)

# =============================================================================
# 2. INPUT HELPERS
# =============================================================================

check_required_file <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }
}

load_sample_metadata <- function(path) {
  check_required_file(path)

  meta <- readr::read_csv(path, show_col_types = FALSE)

  required_cols <- c("SampleID", "MetaGroup")
  missing_cols <- setdiff(required_cols, colnames(meta))

  if (length(missing_cols) > 0) {
    stop(
      "Sample metadata is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  meta %>%
    dplyr::mutate(
      SampleID = as.character(SampleID),
      MetaGroup = as.character(MetaGroup)
    ) %>%
    dplyr::distinct(SampleID, .keep_all = TRUE)
}

# =============================================================================
# 3. CIBERSORTx FRACTION PROCESSING
# =============================================================================

load_cibersort_fractions <- function(path, pvalue_cutoff = 0.05) {
  check_required_file(path)

  ciber <- read.csv(path, check.names = TRUE)

  required_cols <- c("Mixture", "P.value")
  missing_cols <- setdiff(required_cols, colnames(ciber))

  if (length(missing_cols) > 0) {
    stop(
      "CIBERSORTx fraction file is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  ciber <- ciber %>%
    dplyr::filter(P.value < pvalue_cutoff)

  rownames(ciber) <- ciber$Mixture

  ciber <- ciber %>%
    dplyr::select(
      -dplyr::any_of(c("Mixture", "P.value", "RMSE", "Correlation"))
    )

  ciber[] <- lapply(ciber, function(x) suppressWarnings(as.numeric(as.character(x))))

  ciber
}

aggregate_immune_cell_classes <- function(ciber) {
  cell_sets <- list(
    B.cells = c("B.cells.naive", "B.cells.memory", "Plasma.cells"),
    T.cells.CD4 = c(
      "T.cells.CD4.naive",
      "T.cells.CD4.memory.resting",
      "T.cells.CD4.memory.activated",
      "T.cells.regulatory..Tregs."
    ),
    Macrophages = c("Macrophages.M0", "Macrophages.M1", "Macrophages.M2"),
    Dendritic = c("Dendritic.cells.resting", "Dendritic.cells.activated"),
    Mast.cells = c("Mast.cells.resting", "Mast.cells.activated"),
    NK = c("NK.cells.resting", "NK.cells.activated"),
    T.cells.CD8 = "T.cells.CD8"
  )

  missing_cols <- unique(unlist(cell_sets))[!unique(unlist(cell_sets)) %in% colnames(ciber)]

  if (length(missing_cols) > 0) {
    warning(
      "Some expected CIBERSORTx columns are missing and will be treated as zero: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  out <- data.frame(row.names = rownames(ciber))

  for (cell_type in names(cell_sets)) {
    cols <- intersect(cell_sets[[cell_type]], colnames(ciber))

    if (length(cols) == 0) {
      out[[cell_type]] <- 0
    } else {
      out[[cell_type]] <- rowSums(ciber[, cols, drop = FALSE], na.rm = TRUE)
    }
  }

  out
}

# =============================================================================
# 4. GROUP-LEVEL RELATIVE COMPOSITION
# =============================================================================

summarize_group_composition <- function(ciber_agg, sample_meta) {
  meta <- sample_meta %>%
    dplyr::filter(MetaGroup %in% names(group_mapping)) %>%
    dplyr::mutate(GroupPlot = unname(group_mapping[MetaGroup]))

  common_samples <- intersect(rownames(ciber_agg), meta$SampleID)

  if (length(common_samples) == 0) {
    stop("No overlap between CIBERSORTx samples and sample metadata.", call. = FALSE)
  }

  ciber_agg <- ciber_agg[common_samples, , drop = FALSE]
  meta <- meta[match(common_samples, meta$SampleID), , drop = FALSE]

  group_df <- lapply(group_order, function(g) {
    ids <- meta$SampleID[meta$GroupPlot == g]
    ids <- intersect(ids, rownames(ciber_agg))

    if (length(ids) == 0) {
      vals <- rep(NA_real_, ncol(ciber_agg))
    } else {
      vals <- colMeans(ciber_agg[ids, , drop = FALSE], na.rm = TRUE)
    }

    data.frame(
      Group = g,
      CellType = colnames(ciber_agg),
      MeanAbundance = as.numeric(vals),
      stringsAsFactors = FALSE
    )
  }) %>%
    dplyr::bind_rows()

  group_df$Group <- factor(group_df$Group, levels = group_order)

  group_df
}

plot_relative_composition <- function(group_composition, outfile) {
  p <- ggplot2::ggplot(
    group_composition,
    ggplot2::aes(x = Group, y = MeanAbundance, fill = CellType)
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      color = "black",
      width = 0.7,
      position = "fill"
    ) +
    ggplot2::scale_fill_manual(values = immune_palette) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold")
    ) +
    ggplot2::labs(
      x = "Groups",
      y = "Relative abundance",
      fill = "Cell type",
      title = "CIBERSORTx TIME deconvolution"
    )

  ggplot2::ggsave(outfile, p, width = 1600 / 300, height = 1200 / 300, dpi = 300)

  p
}

# =============================================================================
# 5. ABSOLUTE SCORE PROCESSING
# =============================================================================

load_absolute_scores <- function(path, sample_meta, pvalue_cutoff = 0.05) {
  check_required_file(path)

  ciber_abs <- read.csv(path, check.names = TRUE)

  abs_col <- grep("Absolute", colnames(ciber_abs), value = TRUE)[1]

  if (is.na(abs_col)) {
    stop("No Absolute score column found in file: ", path, call. = FALSE)
  }

  ciber_abs <- ciber_abs %>%
    dplyr::filter(P.value < pvalue_cutoff)

  abs_df <- data.frame(
    SampleID = as.character(ciber_abs$Mixture),
    AbsoluteScore = suppressWarnings(as.numeric(as.character(ciber_abs[[abs_col]]))),
    stringsAsFactors = FALSE
  )

  abs_df <- abs_df %>%
    dplyr::inner_join(sample_meta[, c("SampleID", "MetaGroup")], by = "SampleID") %>%
    dplyr::filter(MetaGroup %in% names(group_mapping)) %>%
    dplyr::mutate(
      Group = unname(group_mapping[MetaGroup]),
      Group = factor(Group, levels = group_order)
    )

  abs_df
}

summarize_absolute_scores <- function(abs_df) {
  abs_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      mean_abs = mean(AbsoluteScore, na.rm = TRUE),
      sd_abs = sd(AbsoluteScore, na.rm = TRUE),
      n = sum(!is.na(AbsoluteScore)),
      sem_abs = sd_abs / sqrt(n),
      .groups = "drop"
    )
}

plot_absolute_plus_relative <- function(abs_df, abs_summary, group_composition, outfile) {
  p_abs <- ggplot2::ggplot(abs_summary, ggplot2::aes(x = Group, y = mean_abs)) +
    ggplot2::geom_col(fill = "grey80", color = "black", width = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_abs - sem_abs, ymax = mean_abs + sem_abs),
      width = 0.2,
      linewidth = 0.5
    ) +
    ggplot2::geom_jitter(
      data = abs_df,
      ggplot2::aes(x = Group, y = AbsoluteScore),
      inherit.aes = FALSE,
      width = 0.12,
      size = 1,
      alpha = 0.6
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      y = "Absolute score",
      title = "CIBERSORTx TIME deconvolution"
    )

  p_rel <- ggplot2::ggplot(
    group_composition,
    ggplot2::aes(x = Group, y = MeanAbundance, fill = CellType)
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      color = "black",
      width = 0.7,
      position = "fill"
    ) +
    ggplot2::scale_fill_manual(values = immune_palette) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold")
    ) +
    ggplot2::labs(
      x = "Groups",
      y = "Relative abundance",
      fill = "Cell type"
    )

  p_final <- p_abs / p_rel + patchwork::plot_layout(heights = c(1, 2))

  ggplot2::ggsave(outfile, p_final, width = 2000 / 300, height = 2000 / 300, dpi = 300)

  p_final
}

# =============================================================================
# 6. PERMANOVA AND DISPERSION ANALYSIS
# =============================================================================

prepare_permanova_matrix <- function(ciber_agg, sample_meta) {
  mat <- as.data.frame(lapply(ciber_agg, function(x) as.numeric(as.character(x))))
  rownames(mat) <- rownames(ciber_agg)

  meta <- sample_meta[match(rownames(mat), sample_meta$SampleID), , drop = FALSE]

  keep <- !is.na(meta$MetaGroup) & complete.cases(mat)
  mat <- mat[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  meta <- meta %>%
    dplyr::filter(MetaGroup %in% names(publication_labels)) %>%
    dplyr::mutate(
      GroupPub = publication_labels[MetaGroup],
      GroupPub = factor(GroupPub, levels = publication_group_order)
    )

  mat <- mat[meta$SampleID, , drop = FALSE]

  list(mat = mat, meta = meta)
}

run_global_permanova <- function(mat, meta, permutations = 9999) {
  vegan::adonis2(
    mat ~ GroupPub,
    data = meta,
    method = "bray",
    permutations = permutations
  )
}

run_pairwise_permanova <- function(mat, meta, permutations = 9999) {
  groups <- levels(droplevels(meta$GroupPub))
  pairs <- combn(groups, 2, simplify = FALSE)

  res <- lapply(pairs, function(pair) {
    keep <- meta$GroupPub %in% pair

    mat_sub <- mat[keep, , drop = FALSE]
    meta_sub <- meta[keep, , drop = FALSE]
    meta_sub$GroupPub <- droplevels(meta_sub$GroupPub)

    ad <- vegan::adonis2(
      mat_sub ~ GroupPub,
      data = meta_sub,
      method = "bray",
      permutations = permutations
    )

    data.frame(
      group1 = pair[1],
      group2 = pair[2],
      n1 = sum(meta_sub$GroupPub == pair[1]),
      n2 = sum(meta_sub$GroupPub == pair[2]),
      F = ad$F[1],
      R2 = ad$R2[1],
      p = ad$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  })

  res <- dplyr::bind_rows(res)

  res %>%
    dplyr::mutate(
      p.adj = p.adjust(p, method = "BH"),
      significance = cut(
        p.adj,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("***", "**", "*", "ns")
      )
    )
}

run_pairwise_betadisper <- function(mat, meta, permutations = 9999) {
  groups <- levels(droplevels(meta$GroupPub))
  pairs <- combn(groups, 2, simplify = FALSE)

  res <- lapply(pairs, function(pair) {
    keep <- meta$GroupPub %in% pair

    mat_sub <- mat[keep, , drop = FALSE]
    meta_sub <- meta[keep, , drop = FALSE]
    meta_sub$GroupPub <- droplevels(meta_sub$GroupPub)

    dist_sub <- vegan::vegdist(mat_sub, method = "bray")

    bd <- vegan::betadisper(dist_sub, meta_sub$GroupPub)
    bd_test <- vegan::permutest(bd, permutations = permutations)

    data.frame(
      group1 = pair[1],
      group2 = pair[2],
      F_dispersion = bd_test$tab$F[1],
      p_dispersion = bd_test$tab$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(res) %>%
    dplyr::mutate(
      p_dispersion_adj = p.adjust(p_dispersion, method = "BH")
    )
}

# =============================================================================
# 7. ORDINATION AND STATISTICAL PLOTS
# =============================================================================

plot_pcoa <- function(mat, meta, global_adonis, outfile) {
  dist_mat <- vegan::vegdist(mat, method = "bray")
  pcoa <- cmdscale(dist_mat, eig = TRUE, k = 2)

  eig <- pcoa$eig
  var_explained <- round(100 * eig[1:2] / sum(eig[eig > 0]), 1)

  global_r2 <- global_adonis$R2[1]
  global_p <- global_adonis$`Pr(>F)`[1]

  pcoa_df <- data.frame(
    SampleID = rownames(mat),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2],
    MetaGroup = meta$GroupPub,
    stringsAsFactors = FALSE
  )

  pcoa_df$MetaGroup <- factor(pcoa_df$MetaGroup, levels = publication_group_order)

  p <- ggplot2::ggplot(
    pcoa_df,
    ggplot2::aes(x = PCoA1, y = PCoA2, color = MetaGroup)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = MetaGroup),
      type = "t",
      linetype = 2,
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = group_colors_immune, drop = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(
      title = "Global immune infiltrate composition",
      subtitle = paste0(
        "PERMANOVA: R² = ",
        signif(global_r2, 3),
        ", p = ",
        signif(global_p, 3)
      ),
      x = paste0("PCoA1 (", var_explained[1], "%)"),
      y = paste0("PCoA2 (", var_explained[2], "%)"),
      color = NULL
    )

  ggplot2::ggsave(outfile, p, width = 8, height = 5, dpi = 300)

  p
}

plot_pairwise_permanova_bubbles <- function(combined_res, outfile) {
  plot_df <- combined_res %>%
    dplyr::mutate(
      group1 = factor(group1, levels = publication_group_order),
      group2 = factor(group2, levels = rev(publication_group_order)),
      neglog10_FDR = -log10(p.adj),
      label = sprintf("%.2f", R2)
    )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = group1,
      y = group2,
      size = R2,
      color = neglog10_FDR
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      color = "black",
      size = 7
    ) +
    ggplot2::scale_size_continuous(
      name = expression(R^2),
      range = c(4, 18)
    ) +
    ggplot2::scale_color_gradient(
      low = "#eef4ed",
      high = "#8da9c4",
      name = expression(-log[10](FDR))
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = ggplot2::element_text(hjust = 1, size = 12),
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      legend.key.size = grid::unit(0.3, "cm"),
      legend.spacing.y = grid::unit(0.05, "cm")
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        title.position = "top",
        override.aes = list(alpha = 0.9)
      ),
      color = ggplot2::guide_colorbar(
        title.position = "top",
        barheight = grid::unit(2, "cm"),
        barwidth = grid::unit(0.3, "cm")
      )
    ) +
    ggplot2::labs(
      title = "Pairwise PERMANOVA of immune infiltrate composition",
      subtitle = "Bubble size = R²; color = -log10(FDR)",
      x = NULL,
      y = NULL
    )

  ggplot2::ggsave(outfile, p, width = 6, height = 5, dpi = 250)

  p
}

# =============================================================================
# 8. MAIN WORKFLOW
# =============================================================================

run_cibersort_workflow <- function(config) {
  sample_meta <- load_sample_metadata(config$sample_metadata_file)

  ciber_raw <- load_cibersort_fractions(
    path = config$cibersort_fraction_file,
    pvalue_cutoff = config$pvalue_cutoff
  )

  ciber_agg <- aggregate_immune_cell_classes(ciber_raw)

  readr::write_csv(
    tibble::rownames_to_column(ciber_agg, "SampleID"),
    file.path(config$outdir, "cibersort_aggregated_cell_classes.csv")
  )

  group_composition <- summarize_group_composition(ciber_agg, sample_meta)

  readr::write_csv(
    group_composition,
    file.path(config$outdir, "cibersort_group_relative_composition.csv")
  )

  plot_relative_composition(
    group_composition,
    file.path(config$outdir, "cibersort_overview_relative_composition.png")
  )

  if (file.exists(config$cibersort_absolute_file)) {
    abs_df <- load_absolute_scores(
      path = config$cibersort_absolute_file,
      sample_meta = sample_meta,
      pvalue_cutoff = config$pvalue_cutoff
    )

    abs_summary <- summarize_absolute_scores(abs_df)

    readr::write_csv(abs_df, file.path(config$outdir, "cibersort_absolute_scores.csv"))
    readr::write_csv(abs_summary, file.path(config$outdir, "cibersort_absolute_score_summary.csv"))

    plot_absolute_plus_relative(
      abs_df,
      abs_summary,
      group_composition,
      file.path(config$outdir, "cibersort_overview_with_absolute_score.png")
    )
  } else {
    warning("Absolute score file not found; skipping absolute score panel.")
  }

  perm_obj <- prepare_permanova_matrix(ciber_agg, sample_meta)

  global_adonis <- run_global_permanova(
    mat = perm_obj$mat,
    meta = perm_obj$meta,
    permutations = config$permutations
  )

  write.csv(
    as.data.frame(global_adonis),
    file.path(config$outdir, "cibersort_global_permanova.csv")
  )

  pairwise_res <- run_pairwise_permanova(
    mat = perm_obj$mat,
    meta = perm_obj$meta,
    permutations = config$permutations
  )

  dispersion_res <- run_pairwise_betadisper(
    mat = perm_obj$mat,
    meta = perm_obj$meta,
    permutations = config$permutations
  )

  combined_res <- pairwise_res %>%
    dplyr::left_join(dispersion_res, by = c("group1", "group2"))

  readr::write_csv(
    pairwise_res,
    file.path(config$outdir, "cibersort_pairwise_permanova.csv")
  )

  readr::write_csv(
    dispersion_res,
    file.path(config$outdir, "cibersort_pairwise_betadisper.csv")
  )

  readr::write_csv(
    combined_res,
    file.path(config$outdir, "cibersort_pairwise_permanova_with_dispersion.csv")
  )

  plot_pcoa(
    perm_obj$mat,
    perm_obj$meta,
    global_adonis,
    file.path(config$outdir, "cibersort_pcoa_permanova.png")
  )

  plot_pairwise_permanova_bubbles(
    combined_res,
    file.path(config$outdir, "cibersort_pairwise_permanova_bubbleplot.png")
  )

  message("CIBERSORTx workflow completed. Outputs written to: ", config$outdir)

  invisible(list(
    ciber_aggregated = ciber_agg,
    group_composition = group_composition,
    permanova_global = global_adonis,
    permanova_pairwise = pairwise_res,
    dispersion_pairwise = dispersion_res,
    permanova_combined = combined_res
  ))
}

# =============================================================================
# 9. RUN
# =============================================================================

results <- run_cibersort_workflow(config)
