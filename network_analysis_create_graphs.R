#!/usr/bin/env Rscript

############################################################
## NETWORK ANALYSIS PIPELINE (single chunk, reviewer-hardened)
## CLI:
##   --groups <groups.xlsx> --tpm <tpm.tsv> --counts <counts.tsv>
##   [--out_prefix DCN_graphs] [--overwrite]
##   [--alpha_group 0.01] [--t_percentile 0.90] [--rho NULL]
##   [--min_overlap 8]
##   [--diff_alpha 0.01] [--diff_t 0.2] [--diff_fdr_method BH]
##   [--permute_B 0] [--bootstrap_B 0]
##
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readxl)
  library(matrixStats)
})

## =========================================================
## CLI parsing (no optparse dependency)
## =========================================================
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}
has_flag <- function(flag) flag %in% args

as_num <- function(x, default=NULL) {
  if (is.null(x)) return(default)
  if (identical(x, "NULL")) return(NULL)
  suppressWarnings({
    v <- as.numeric(x)
    if (is.na(v)) default else v
  })
}
as_int <- function(x, default=0L) {
  if (is.null(x)) return(default)
  suppressWarnings({
    v <- as.integer(x)
    if (is.na(v)) default else v
  })
}

groups_xlsx <- get_arg("--groups", NULL)
tpm_file    <- get_arg("--tpm", NULL)
counts_file <- get_arg("--counts", NULL)
out_prefix  <- get_arg("--out_prefix", "DCN_graphs")
overwrite   <- has_flag("--overwrite")

alpha_group   <- as_num(get_arg("--alpha_group", "0.01"), 0.01)
t_percentile  <- as_num(get_arg("--t_percentile", "0.90"), 0.90)
rho_fixed     <- get_arg("--rho", "NULL")
rho_fixed     <- if (identical(rho_fixed, "NULL")) NULL else as_num(rho_fixed, NULL)

min_overlap   <- as_int(get_arg("--min_overlap", "8"), 8L)

diff_alpha    <- as_num(get_arg("--diff_alpha", "0.01"), 0.01)
diff_t        <- as_num(get_arg("--diff_t", "0.2"), 0.2)   # secondary effect size only
diff_fdr_m    <- get_arg("--diff_fdr_method", "BH")

permute_B     <- as_int(get_arg("--permute_B", "0"), 0L)   # 0 disables
bootstrap_B   <- as_int(get_arg("--bootstrap_B", "0"), 0L) # 0 disables

if (is.null(groups_xlsx) || is.null(tpm_file) || is.null(counts_file)) {
  stop("Usage: Rscript network_analysis.R --groups <groups.xlsx> --tpm <tpm.tsv> --counts <counts.tsv> [--out_prefix DCN_graphs] [--overwrite] [--alpha_group 0.01] [--t_percentile 0.90] [--rho NULL] [--min_overlap 8] [--diff_alpha 0.01] [--diff_t 0.2] [--diff_fdr_method BH] [--permute_B 0] [--bootstrap_B 0]")
}

exclude_groups <- c("iperA","iperT")  # hyperplasia exclusions

stopifnot(file.exists(groups_xlsx), file.exists(tpm_file), file.exists(counts_file))

`%||%` <- function(a,b) if (!is.null(a)) a else b

## =========================================================
## Helpers: IO + sample parsing
## =========================================================
clean_sample_id <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x
}

read_salmon_merged <- function(path){
  df <- read.csv2(path, row.names = 1, check.names = FALSE, sep = "\t")
  if ("gene_id" %in% colnames(df)) df$gene_id <- NULL
  stopifnot("gene_name" %in% colnames(df))

  num_cols <- setdiff(colnames(df), "gene_name")
  df[num_cols] <- lapply(df[num_cols], function(v) as.numeric(as.character(v)))

  df_agg <- aggregate(
    df[, num_cols, drop = FALSE],
    by = list(gene_name = df$gene_name),
    FUN = sum,
    na.rm = TRUE
  )
  rownames(df_agg) <- df_agg$gene_name
  df_agg$gene_name <- NULL

  df_agg <- df_agg[rowSums(df_agg, na.rm = TRUE) >= 1, , drop = FALSE]
  colnames(df_agg) <- clean_sample_id(colnames(df_agg))
  df_agg
}

infer_metagroup_from_id <- function(sample_id){
  sid <- as.character(sample_id)
  suffix <- sub("^.*\\.", "", sid)

  if (grepl("^Dec$",   suffix, ignore.case = TRUE)) return("decidua")
  if (grepl("^Troph$", suffix, ignore.case = TRUE)) return("trophoblast")

  if (grepl("\\.sADK$", sid, ignore.case = TRUE)) return("Healthy")
  if (grepl("(^|\\.|_)san($|\\.|_)", sid, ignore.case = TRUE)) return("Healthy")

  if (grepl("\\.tADK$", sid, ignore.case = TRUE)) return("EC")

  NA_character_
}

drop_na_cols <- function(mat){
  mat[, colSums(is.na(mat)) == 0, drop = FALSE]
}

subset_group_matrix <- function(mat, sample_meta, group_name){
  ids <- sample_meta %>%
    filter(Group %in% group_name) %>%
    pull(SampleID) %>%
    intersect(colnames(mat))

  if (length(ids) < 3) {
    warning("Group '", group_name, "' has <3 samples; correlations may be unstable.")
  }

  m <- mat[, ids, drop = FALSE]
  m <- drop_na_cols(m)
  if (ncol(m) == 0) return(matrix(numeric(0), nrow = 0, ncol = 0))
  m <- m[rowMeans(m, na.rm = TRUE) > 1, , drop = FALSE]
  m
}

## =========================================================
## Load matrices
## =========================================================
tpm_raw    <- read_salmon_merged(tpm_file)
counts_raw <- read_salmon_merged(counts_file)

common_keep <- intersect(colnames(tpm_raw), colnames(counts_raw))
if (length(common_keep) < 2) stop("TPM/Counts overlap < 2 samples")

tpm_raw    <- tpm_raw[, common_keep, drop = FALSE]
counts_raw <- counts_raw[, common_keep, drop = FALSE]

## =========================================================
## Load groups.xlsx (audit only)
## =========================================================
groups_xlsx_df <- readxl::read_excel(groups_xlsx) %>% as.data.frame()
stopifnot("ID" %in% colnames(groups_xlsx_df))

groups_mfp <- groups_xlsx_df %>%
  transmute(SampleID_raw = as.character(ID),
            SampleID     = clean_sample_id(as.character(ID))) %>%
  distinct(SampleID, .keep_all = TRUE)

## =========================================================
## Build sample_meta
## =========================================================
sample_meta <- tibble(SampleID = common_keep) %>%
  mutate(
    SampleID = as.character(SampleID),
    PairID   = sub("\\..*$", "", SampleID),
    Suffix   = sub("^.*\\.", "", SampleID),
    SuffixToken = vapply(strsplit(SampleID, ".", fixed = TRUE), function(v) tail(v, 1), character(1)),
    IsHyperplasia = SuffixToken %in% exclude_groups
  ) %>%
  filter(!IsHyperplasia) %>%
  mutate(
    MetaGroup = vapply(SampleID, infer_metagroup_from_id, character(1)),
    Group = MetaGroup
  ) %>%
  filter(!is.na(Group)) %>%
  distinct(SampleID, .keep_all = TRUE)

if (nrow(sample_meta) < 2) stop("Too few samples after filtering by exclude_groups/MetaGroup mapping.")

keep_ids <- intersect(sample_meta$SampleID, colnames(tpm_raw))
tpm_raw    <- tpm_raw[, keep_ids, drop = FALSE]
counts_raw <- counts_raw[, keep_ids, drop = FALSE]

sample_meta <- sample_meta %>%
  filter(SampleID %in% keep_ids) %>%
  mutate(SampleID = factor(SampleID, levels = keep_ids)) %>%
  arrange(SampleID) %>%
  mutate(SampleID = as.character(SampleID))

message("sample_meta groups:\n")
print(table(sample_meta$Group, useNA = "ifany"))

## =========================================================
## Gene expression matrix (TPM)
## =========================================================
gene_expression_matrix <- as.matrix(tpm_raw)
mode(gene_expression_matrix) <- "numeric"

if (any(duplicated(rownames(gene_expression_matrix)))) {
  gene_expression_matrix <- rowsum(gene_expression_matrix,
                                  group = rownames(gene_expression_matrix),
                                  reorder = FALSE)
}

gene_expression_matrix_decidua <- subset_group_matrix(gene_expression_matrix, sample_meta, "decidua")
gene_expression_matrix_healthy <- subset_group_matrix(gene_expression_matrix, sample_meta, "Healthy")
gene_expression_matrix_EC      <- subset_group_matrix(gene_expression_matrix, sample_meta, "EC")

message("Matrix dims (genes x samples):")
message("  decidua: ", nrow(gene_expression_matrix_decidua), " x ", ncol(gene_expression_matrix_decidua))
message("  Healthy: ", nrow(gene_expression_matrix_healthy), " x ", ncol(gene_expression_matrix_healthy))
message("  EC:      ", nrow(gene_expression_matrix_EC),      " x ", ncol(gene_expression_matrix_EC))

## =========================================================
## DCN builder (scran correlations; Fisher-Z differential testing)
## =========================================================
make_diff_graph <- function(df1,
                            df2,
                            alpha = 0.01,
                            rho = NULL,
                            t_percentile = 0.90,
                            min_overlap = 8,
                            require_both_above_t = FALSE,
                            diff_alpha = 0.01,
                            diff_t = 0.2,
                            fdr_method = c("BH","bonferroni","holm"),
                            diff_fdr_method = c("BH","bonferroni","holm"),
                            keep_all_pairs_for_diff = TRUE) {

  fdr_method      <- match.arg(fdr_method)
  diff_fdr_method <- match.arg(diff_fdr_method)

  if (!requireNamespace("scran", quietly = TRUE))
    stop("Install 'scran' (BiocManager::install('scran')).")
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Install 'igraph' (install.packages('igraph')).")

  keep_numeric_cols <- function(x) {
    x <- as.data.frame(x, check.names = FALSE)
    num <- vapply(x, function(v) is.numeric(v) || is.integer(v), logical(1))
    x[, num, drop = FALSE]
  }
  sd0 <- function(v) isTRUE(all.equal(stats::sd(v, na.rm = TRUE), 0))

  if (is.null(rownames(df1)) || is.null(rownames(df2)))
    stop("Both inputs must have rownames = gene IDs.")

  df1 <- keep_numeric_cols(df1)
  df2 <- keep_numeric_cols(df2)

  common <- sort(intersect(rownames(df1), rownames(df2)))
  if (length(common) < 2) stop("Need >=2 shared genes.")

  M1 <- as.matrix(df1[common, , drop = FALSE])
  M2 <- as.matrix(df2[common, , drop = FALSE])

  ok1 <- rowSums(!is.na(M1)) > 1 & !apply(M1, 1, sd0)
  ok2 <- rowSums(!is.na(M2)) > 1 & !apply(M2, 1, sd0)
  keep <- ok1 & ok2
  M1 <- M1[keep, , drop = FALSE]
  M2 <- M2[keep, , drop = FALSE]

  genes <- rownames(M1)
  if (length(genes) < 2) stop("Less than 2 genes after filtering.")

  # scran computes correlations + overlaps (+ p-values in recent versions)
  tab1 <- scran::correlatePairs(M1)
  tab2 <- scran::correlatePairs(M2)

  extract_pairs <- function(tab, gnames) {
    nm <- names(tab)
    asDF <- as.data.frame(tab)

    if (all(c("first","second") %in% nm)) {
      i <- as.integer(asDF$first); j <- as.integer(asDF$second)
      g1 <- gnames[i]; g2 <- gnames[j]
    } else if (all(c("row","column") %in% nm)) {
      i <- as.integer(asDF$row); j <- as.integer(asDF$column)
      g1 <- gnames[i]; g2 <- gnames[j]
    } else if (all(c("gene1","gene2") %in% nm)) {
      g1 <- as.character(asDF$gene1); g2 <- as.character(asDF$gene2)
    } else {
      stop("Unexpected correlatePairs() columns: ", paste(nm, collapse=", "))
    }

    if (!"rho" %in% nm) stop("Missing 'rho' column from correlatePairs().")
    rho <- as.numeric(asDF$rho)
    rho <- pmin(pmax(rho, -0.999999), 0.999999)

    n_over <- if ("n" %in% nm) as.integer(asDF$n) else NA_integer_
    pval   <- if ("p.value" %in% nm) as.numeric(asDF$p.value) else NA_real_
    fdr    <- if ("FDR" %in% nm) as.numeric(asDF$FDR) else NA_real_

    out <- data.frame(gene1=g1, gene2=g2, rho=rho, n=n_over, pval=pval, FDR=fdr,
                      stringsAsFactors = FALSE)

    if (all(is.na(out$FDR)) && !all(is.na(out$pval)))
      out$FDR <- p.adjust(out$pval, method = fdr_method)

    # normalize ordering
    swap <- out$gene1 > out$gene2
    if (any(swap)) {
      tmp <- out$gene1[swap]; out$gene1[swap] <- out$gene2[swap]; out$gene2[swap] <- tmp
    }
    out
  }

  e1_raw <- extract_pairs(tab1, genes)
  e2_raw <- extract_pairs(tab2, genes)

  # Threshold for condition graphs (still heuristic but explicit)
  t_thr <- if (is.null(rho)) {
    as.numeric(stats::quantile(c(abs(e1_raw$rho), abs(e2_raw$rho)),
                               probs = t_percentile, na.rm = TRUE))
  } else rho

  filt_group <- function(E) {
    subset(E,
           is.finite(rho) & abs(rho) >= t_thr &
             (is.na(FDR) | FDR < alpha) &
             (is.na(n)   | n >= min_overlap))
  }

  e1 <- filt_group(e1_raw); names(e1)[names(e1)=="rho"] <- "rho1"
  e2 <- filt_group(e2_raw); names(e2)[names(e2)=="rho"] <- "rho2"

  if (require_both_above_t) {
    e12 <- merge(e1, e2, by=c("gene1","gene2"), all=FALSE)
    e1 <- e12[, c("gene1","gene2","rho1","FDR.x","n.x")]; names(e1)[4:5] <- c("FDR1","n1")
    e2 <- e12[, c("gene1","gene2","rho2","FDR.y","n.y")]; names(e2)[4:5] <- c("FDR2","n2")
  } else {
    names(e1)[names(e1)=="FDR"] <- "FDR1"; names(e1)[names(e1)=="n"] <- "n1"
    names(e2)[names(e2)=="FDR"] <- "FDR2"; names(e2)[names(e2)=="n"] <- "n2"
  }

  # Differential testing:
  # Prefer testing across *all* comparable pairs (not only intersection of significant edges),
  # then apply FDR + optional effect-size filter for interpretability.
  # This removes dependence on arbitrary |Δρ| alone.
  if (isTRUE(keep_all_pairs_for_diff)) {
    dtab <- merge(e1_raw, e2_raw, by=c("gene1","gene2"), suffixes=c(".1",".2"), all=FALSE)
    if (nrow(dtab)) {
      rho1_all <- pmin(pmax(as.numeric(dtab$rho.1), -0.999999), 0.999999)
      rho2_all <- pmin(pmax(as.numeric(dtab$rho.2), -0.999999), 0.999999)
      n1_all   <- if ("n.1" %in% names(dtab)) as.integer(dtab$n.1) else NA_integer_
      n2_all   <- if ("n.2" %in% names(dtab)) as.integer(dtab$n.2) else NA_integer_

      keep_n <- (is.na(n1_all) | n1_all >= min_overlap) & (is.na(n2_all) | n2_all >= min_overlap)
      dtab <- dtab[keep_n, , drop=FALSE]
      rho1_all <- rho1_all[keep_n]; rho2_all <- rho2_all[keep_n]
      n1_all <- n1_all[keep_n]; n2_all <- n2_all[keep_n]

      z1 <- atanh(rho1_all); z2 <- atanh(rho2_all)
      n1_eff <- if (!all(is.na(n1_all))) n1_all else ncol(M1)
      n2_eff <- if (!all(is.na(n2_all))) n2_all else ncol(M2)
      se_diff <- sqrt(1 / pmax(1, n1_eff - 3) + 1 / pmax(1, n2_eff - 3))
      zstat <- (z1 - z2) / se_diff
      pval  <- 2 * stats::pnorm(-abs(zstat))
      padj  <- stats::p.adjust(pval, method = diff_fdr_method)
      delta <- rho1_all - rho2_all

      eD_all <- data.frame(
        gene1 = dtab$gene1, gene2 = dtab$gene2,
        rho1 = rho1_all, rho2 = rho2_all,
        delta = delta,
        z = zstat, pval = pval, padj = padj,
        n1 = n1_eff, n2 = n2_eff,
        stringsAsFactors = FALSE
      )

      eD <- subset(eD_all,
                   is.finite(padj) & padj < diff_alpha &
                     (is.null(diff_t) | abs(delta) >= diff_t))
    } else {
      eD_all <- data.frame()
      eD <- data.frame()
    }
  } else {
    # legacy: differential on intersection of retained group edges only
    if (nrow(e1) && nrow(e2)) {
      md <- merge(e1, e2, by=c("gene1","gene2"), all=FALSE)
      if (nrow(md)) {
        z1 <- atanh(md$rho1); z2 <- atanh(md$rho2)
        n1_eff <- if (!all(is.na(md$n1))) md$n1 else ncol(M1)
        n2_eff <- if (!all(is.na(md$n2))) md$n2 else ncol(M2)
        se_diff <- sqrt(1 / pmax(1, n1_eff - 3) + 1 / pmax(1, n2_eff - 3))
        zstat   <- (z1 - z2) / se_diff
        pval    <- 2 * stats::pnorm(-abs(zstat))
        padj    <- stats::p.adjust(pval, method = diff_fdr_method)
        md$delta <- md$rho1 - md$rho2
        md$z <- zstat; md$pval <- pval; md$padj <- padj
        eD <- subset(md, padj < diff_alpha & (is.null(diff_t) | abs(delta) >= diff_t),
                     select=c("gene1","gene2","rho1","rho2","delta","z","pval","padj"))
      } else eD <- data.frame()
    } else eD <- data.frame()
    eD_all <- NULL
  }

  if (!requireNamespace("igraph", quietly = TRUE)) stop("igraph required")

  g1 <- igraph::graph_from_data_frame(
    if (nrow(e1)) e1[,c("gene1","gene2","rho1")] else data.frame(gene1=character(),gene2=character(),rho1=numeric()),
    directed=FALSE, vertices=genes
  )
  if (nrow(e1)) igraph::E(g1)$weight <- e1$rho1

  g2 <- igraph::graph_from_data_frame(
    if (nrow(e2)) e2[,c("gene1","gene2","rho2")] else data.frame(gene1=character(),gene2=character(),rho2=numeric()),
    directed=FALSE, vertices=genes
  )
  if (nrow(e2)) igraph::E(g2)$weight <- e2$rho2

  gd <- igraph::graph_from_data_frame(
    if (nrow(eD)) eD[,c("gene1","gene2","delta")] else data.frame(gene1=character(),gene2=character(),delta=numeric()),
    directed=FALSE, vertices=genes
  )
  if (nrow(eD)) igraph::E(gd)$weight <- eD$delta

  list(
    threshold_t      = t_thr,
    alpha_group      = alpha,
    diff_alpha       = diff_alpha,
    params = list(
      t_percentile = t_percentile,
      min_overlap = min_overlap,
      require_both_above_t = require_both_above_t,
      diff_t = diff_t,
      fdr_method = fdr_method,
      diff_fdr_method = diff_fdr_method,
      keep_all_pairs_for_diff = keep_all_pairs_for_diff
    ),
    edges_group1 = e1,
    edges_group2 = e2,
    edges_difference = eD,
    edges_difference_all = eD_all,
    graph_group1     = g1,
    graph_group2     = g2,
    graph_difference = gd
  )
}

## =========================================================
## Cache / skip (digest-based; no overflow)
## =========================================================
if (!requireNamespace("digest", quietly = TRUE))
  stop("Install 'digest' for caching signatures: install.packages('digest')")

graph_cache_signature <- function(label1, label2, M1, M2, params) {
  rn1 <- rownames(M1); rn2 <- rownames(M2)
  cn1 <- colnames(M1); cn2 <- colnames(M2)
  sig_list <- list(
    pair = paste(label1, "vs", label2),
    dims = list(M1 = c(nrow(M1), ncol(M1)), M2 = c(nrow(M2), ncol(M2))),
    rn1  = if (!is.null(rn1) && length(rn1)) c(length(rn1), rn1[1], rn1[length(rn1)]) else NA,
    rn2  = if (!is.null(rn2) && length(rn2)) c(length(rn2), rn2[1], rn2[length(rn2)]) else NA,
    cn1  = if (!is.null(cn1) && length(cn1)) c(length(cn1), cn1[1], cn1[length(cn1)]) else NA,
    cn2  = if (!is.null(cn2) && length(cn2)) c(length(cn2), cn2[1], cn2[length(cn2)]) else NA,
    params = params
  )
  digest::digest(sig_list, algo = "xxhash64")
}

## ---- optional stability diagnostics (slow; opt-in) ----
bootstrap_edge_sd <- function(g, M, B = 200L) {
  # g: igraph with V names = genes; M: genes x samples matrix
  if (B <= 0) return(NULL)
  if (!inherits(g, "igraph") || igraph::ecount(g) == 0) return(NULL)
  genes <- unique(as.vector(igraph::as_edgelist(g, names = TRUE)))
  genes <- intersect(genes, rownames(M))
  if (length(genes) < 2) return(NULL)
  Mb <- M[genes, , drop=FALSE]
  el <- igraph::as_edgelist(g, names = TRUE)
  i1 <- match(el[,1], rownames(Mb))
  i2 <- match(el[,2], rownames(Mb))
  n <- ncol(Mb)
  vals <- matrix(NA_real_, nrow = B, ncol = nrow(el))
  for (b in seq_len(B)) {
    cols <- sample.int(n, size=n, replace=TRUE)
    C <- stats::cor(t(Mb[, cols, drop=FALSE]), use="pairwise.complete.obs")
    vals[b,] <- C[cbind(i1, i2)]
  }
  data.frame(
    gene1 = el[,1],
    gene2 = el[,2],
    rho_mean = colMeans(vals, na.rm=TRUE),
    rho_sd = apply(vals, 2, stats::sd, na.rm=TRUE),
    stringsAsFactors = FALSE
  )
}

permute_global_rewiring <- function(M_all, idxA, idxB, params, B = 200L) {
  # Very expensive because it recomputes correlations each time.
  if (B <= 0) return(NULL)
  nA <- length(idxA); nB <- length(idxB)
  all_idx <- c(idxA, idxB)

  obs_obj <- make_diff_graph(M_all[, idxA, drop=FALSE], M_all[, idxB, drop=FALSE],
                             alpha = params$alpha, rho = params$rho, t_percentile = params$t_percentile,
                             min_overlap = params$min_overlap,
                             require_both_above_t = params$require_both_above_t,
                             diff_alpha = params$diff_alpha, diff_t = params$diff_t,
                             fdr_method = params$fdr_method, diff_fdr_method = params$diff_fdr_method,
                             keep_all_pairs_for_diff = TRUE)
  obs <- sum(abs(obs_obj$edges_difference$delta), na.rm=TRUE)

  perm <- numeric(B)
  for (b in seq_len(B)) {
    perm_ids <- sample(all_idx)
    pA <- perm_ids[seq_len(nA)]
    pB <- perm_ids[(nA+1):(nA+nB)]
    obj <- make_diff_graph(M_all[, pA, drop=FALSE], M_all[, pB, drop=FALSE],
                           alpha = params$alpha, rho = params$rho, t_percentile = params$t_percentile,
                           min_overlap = params$min_overlap,
                           require_both_above_t = params$require_both_above_t,
                           diff_alpha = params$diff_alpha, diff_t = params$diff_t,
                           fdr_method = params$fdr_method, diff_fdr_method = params$diff_fdr_method,
                           keep_all_pairs_for_diff = TRUE)
    perm[b] <- sum(abs(obj$edges_difference$delta), na.rm=TRUE)
  }
  p_emp <- (1 + sum(perm >= obs)) / (B + 1)
  list(observed = obs, perm = perm, p_emp = p_emp)
}

run_or_load_diff_graph <- function(label1, label2, M1, M2, out_prefix,
                                   params = list(),
                                   overwrite = FALSE,
                                   bootstrap_B = 0L,
                                   permute_B = 0L) {

  sig <- graph_cache_signature(label1, label2, M1, M2, params)

  f_g1   <- sprintf("%s_%s_%s_group1.rds",     out_prefix, sig, label1)
  f_g2   <- sprintf("%s_%s_%s_group2.rds",     out_prefix, sig, label2)
  f_gd   <- sprintf("%s_%s_%s_vs_%s_diff.rds", out_prefix, sig, label1, label2)
  f_meta <- sprintf("%s_%s_%s_vs_%s_meta.rds", out_prefix, sig, label1, label2)

  f_edges_all <- sprintf("%s_%s_%s_vs_%s_diff_edges_all.tsv.gz", out_prefix, sig, label1, label2)
  f_edges_sig <- sprintf("%s_%s_%s_vs_%s_diff_edges_sig.tsv.gz", out_prefix, sig, label1, label2)

  f_boot1 <- sprintf("%s_%s_%s_bootstrap_edges.tsv.gz", out_prefix, sig, label1)
  f_boot2 <- sprintf("%s_%s_%s_bootstrap_edges.tsv.gz", out_prefix, sig, label2)

  exists_all <- file.exists(f_g1) && file.exists(f_g2) && file.exists(f_gd) && file.exists(f_meta)

  if (exists_all && !isTRUE(overwrite)) {
    message(sprintf("[cache] Using existing graphs for %s vs %s (sig=%s)", label1, label2, sig))
    return(list(
      signature = sig,
      files = list(g1=f_g1, g2=f_g2, gd=f_gd, meta=f_meta,
                   edges_all=f_edges_all, edges_sig=f_edges_sig,
                   boot1=f_boot1, boot2=f_boot2),
      meta = readRDS(f_meta),
      graph_group1 = readRDS(f_g1),
      graph_group2 = readRDS(f_g2),
      graph_difference = readRDS(f_gd),
      skipped = TRUE
    ))
  }

  message(sprintf("[run] Building graphs for %s vs %s (sig=%s)", label1, label2, sig))

  diff_obj <- do.call(make_diff_graph, c(list(df1 = M1, df2 = M2), params))

  saveRDS(diff_obj$graph_group1,     f_g1)
  saveRDS(diff_obj$graph_group2,     f_g2)
  saveRDS(diff_obj$graph_difference, f_gd)

  # Persist edge tables (useful for downstream + transparency)
  if (!is.null(diff_obj$edges_difference_all) && nrow(diff_obj$edges_difference_all)) {
    tmp <- diff_obj$edges_difference_all
    tmp <- tmp[order(tmp$padj, -abs(tmp$delta)), , drop=FALSE]
    gz <- gzfile(f_edges_all, "wt"); write.table(tmp, gz, sep="\t", quote=FALSE, row.names=FALSE); close(gz)
  }
  if (nrow(diff_obj$edges_difference)) {
    tmp <- diff_obj$edges_difference
    tmp <- tmp[order(tmp$padj, -abs(tmp$delta)), , drop=FALSE]
    gz <- gzfile(f_edges_sig, "wt"); write.table(tmp, gz, sep="\t", quote=FALSE, row.names=FALSE); close(gz)
  }

  # Optional bootstrap stability for edges in the retained condition graphs
  boot1 <- boot2 <- NULL
  if (bootstrap_B > 0L) {
    message(sprintf("[stability] Bootstrapping edge SD (%s, B=%d)", label1, bootstrap_B))
    boot1 <- bootstrap_edge_sd(diff_obj$graph_group1, M1, B=bootstrap_B)
    if (!is.null(boot1) && nrow(boot1)) {
      gz <- gzfile(f_boot1, "wt"); write.table(boot1, gz, sep="\t", quote=FALSE, row.names=FALSE); close(gz)
    }
    message(sprintf("[stability] Bootstrapping edge SD (%s, B=%d)", label2, bootstrap_B))
    boot2 <- bootstrap_edge_sd(diff_obj$graph_group2, M2, B=bootstrap_B)
    if (!is.null(boot2) && nrow(boot2)) {
      gz <- gzfile(f_boot2, "wt"); write.table(boot2, gz, sep="\t", quote=FALSE, row.names=FALSE); close(gz)
    }
  }

  # Optional permutation test for global rewiring burden (very slow)
  perm_res <- NULL
  if (permute_B > 0L) {
    message(sprintf("[permute] Permutation test for global rewiring (B=%d) — may be slow", permute_B))
    # Build combined matrix for permutation (genes x samples)
    common_genes <- intersect(rownames(M1), rownames(M2))
    M_all <- cbind(M1[common_genes,,drop=FALSE], M2[common_genes,,drop=FALSE])
    idxA <- seq_len(ncol(M1))
    idxB <- (ncol(M1)+1):(ncol(M1)+ncol(M2))
    perm_res <- permute_global_rewiring(M_all, idxA, idxB, params=params, B=permute_B)
  }

  meta <- list(
    signature = sig,
    label1 = label1,
    label2 = label2,
    params = params,
    threshold_t = diff_obj$threshold_t,
    alpha_group = diff_obj$alpha_group,
    diff_alpha  = diff_obj$diff_alpha,
    timestamp = as.character(Sys.time()),
    dims = list(M1 = c(genes=nrow(M1), samples=ncol(M1)),
                M2 = c(genes=nrow(M2), samples=ncol(M2))),
    outputs = list(
      edges_all = if (file.exists(f_edges_all)) f_edges_all else NA,
      edges_sig = if (file.exists(f_edges_sig)) f_edges_sig else NA,
      bootstrap1 = if (file.exists(f_boot1)) f_boot1 else NA,
      bootstrap2 = if (file.exists(f_boot2)) f_boot2 else NA
    ),
    permutation = perm_res
  )
  saveRDS(meta, f_meta)

  list(
    signature = sig,
    files = list(g1=f_g1, g2=f_g2, gd=f_gd, meta=f_meta,
                 edges_all=f_edges_all, edges_sig=f_edges_sig,
                 boot1=f_boot1, boot2=f_boot2),
    meta = meta,
    graph_group1 = diff_obj$graph_group1,
    graph_group2 = diff_obj$graph_group2,
    graph_difference = diff_obj$graph_difference,
    skipped = FALSE
  )
}

## =========================================================
## X2) Run DCN (cached)
## =========================================================
dcn_params <- list(
  alpha = alpha_group,
  rho = 0.7,
  t_percentile = t_percentile,
  min_overlap = min_overlap,
  require_both_above_t = FALSE,
  diff_alpha = diff_alpha,
  diff_t = diff_t,
  fdr_method = "BH",
  diff_fdr_method = diff_fdr_m,
  keep_all_pairs_for_diff = TRUE
)

res_EC_vs_Healthy <- run_or_load_diff_graph(
  label1 = "EC", label2 = "Healthy",
  M1 = gene_expression_matrix_EC,
  M2 = gene_expression_matrix_healthy,
  out_prefix = out_prefix,
  params = dcn_params,
  overwrite = overwrite,
  bootstrap_B = bootstrap_B,
  permute_B = permute_B
)

res_Dec_vs_Healthy <- run_or_load_diff_graph(
  label1 = "Decidua", label2 = "Healthy",
  M1 = gene_expression_matrix_decidua,
  M2 = gene_expression_matrix_healthy,
  out_prefix = out_prefix,
  params = dcn_params,
  overwrite = overwrite,
  bootstrap_B = bootstrap_B,
  permute_B = permute_B
)

res_EC_vs_Dec <- run_or_load_diff_graph(
  label1 = "EC", label2 = "Decidua",
  M1 = gene_expression_matrix_EC,
  M2 = gene_expression_matrix_decidua,
  out_prefix = out_prefix,
  params = dcn_params,
  overwrite = overwrite,
  bootstrap_B = bootstrap_B,
  permute_B = permute_B
)

## Legacy fixed filenames (optional convenience)
saveRDS(res_EC_vs_Healthy$graph_group1,     file = "graph_groupEC.rds")
saveRDS(res_EC_vs_Healthy$graph_group2,     file = "graph_groupHealthy.rds")
saveRDS(res_EC_vs_Healthy$graph_difference, file = "graph_difference_ECvsHealthy.rds")

saveRDS(res_EC_vs_Dec$graph_group1,         file = "graph_groupEC_ECvsDecidua.rds")
saveRDS(res_EC_vs_Dec$graph_group2,         file = "graph_groupDecidua.rds")
saveRDS(res_EC_vs_Dec$graph_difference,     file = "graph_difference_ECvsDecidua.rds")

message("\n[done] Upstream graphs written. Downstream analyses should load the cached RDS + *_meta.rds + *_edges_*.tsv.gz\n")
