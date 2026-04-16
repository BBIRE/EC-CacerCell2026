## =========================================================
## RNA-seq unified pipeline
## + CIBERSORT + ICI complexes (MFI paired) + plots
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readxl)
  library(ggplot2)
  library(reshape2)
  library(DESeq2)
  library(msigdbr)
  library(fgsea)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(purrr)
  library(ggpubr)
  library(hrbrthemes)
  library(matrixStats)
  library(SummarizedExperiment)
})

set.seed(1234)
options(stringsAsFactors = FALSE)

## ---------------------------------------------------------
## Resolve masking: AnnotationDbi::select() vs dplyr::select()
## ---------------------------------------------------------
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
rename    <- dplyr::rename
arrange   <- dplyr::arrange
distinct  <- dplyr::distinct
summarise <- dplyr::summarise
pull      <- dplyr::pull
left_join <- dplyr::left_join
right_join <- dplyr::right_join
inner_join <- dplyr::inner_join
full_join  <- dplyr::full_join

## =========================================================
## 0) CONFIG
## =========================================================
cfg <- list(
  groups_xlsx = "groups.xlsx",
  tpm_file    = "final.salmon.merged.gene_tpm.tsv",
  counts_file = "final.salmon.merged.gene_counts.tsv",
  ciber_file  = "CIBERSORTx_Job100_Results.csv",

  out_dir    = "results",
  out_tables = "results/tables",
  out_figs   = "results/figures",
  out_logs   = "results/logs",

  ciber_p_cutoff = 0.05,

  padj_cutoff = 0.05,
  lfc_cutoff  = 1,

  enrichment = list(
    excel_dir      = "results/tables",
    file_pattern   = "DESeq2_",
    species        = "Homo sapiens",
    go_subcategory = "BP",
    top_n_per_comp = 25,
    seed           = 1234
  )
)

dir.create(cfg$out_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$out_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$out_figs,   showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$out_logs,   showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()), file.path(cfg$out_logs, "sessionInfo.txt"))

## =========================================================
## 1) Canonical group dictionaries (ONLY EDIT HERE)
## =========================================================

## IMMUNE/plot view (single EC group)
group_levels_immune <- c(
  "Healthy - Endometrium Healthy Tissue" = 1,
  "Healthy - EC Adjacent Healthy Tissue" = 2,
  "Decidua - MFI"                        = 3,
  "EC - TIME"                            = 4
)

group_colors_immune <- c(
  "Healthy - Endometrium Healthy Tissue" = "#1B6294",
  "Healthy - EC Adjacent Healthy Tissue" = "#72B9ED",
  "Decidua - MFI"                        = "#7B9E87",
  "EC - TIME"                            = "#7A0B12"
)

## Tissue view (biology; keep trophoblast here if you want in tissue-only plots)
group_levels_tissue <- c(
  "Healthy_Endometrium" = 1,
  "Healthy_EC"          = 2,
  "Decidua"             = 3,
  "Trophoblast"         = 4,
  "EC"                  = 5
)

group_colors_tissue <- c(
  "Healthy_Endometrium" = "#1B6294",
  "Healthy_EC"          = "#72B9ED",
  "Decidua"             = "#7B9E87",
  "Trophoblast"         = "#2B3F33",
  "EC"                  = "#7A0B12"
)

plot_group_order <- c("Healthy_Endometrium","Healthy_EC","Decidua","Trophoblast","EC")
plot_group_cols  <- group_colors_tissue

## =========================================================
## 2) Helpers
## =========================================================

clean_sample_id <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x
}

diagnose_overlap <- function(a, b, label = "overlap"){
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  common <- intersect(a, b)
  message(label, ": common n = ", length(common))
  if (length(common) < 2) {
    message("In A not B (first 30): ", paste(head(setdiff(a, b), 30), collapse = ", "))
    message("In B not A (first 30): ", paste(head(setdiff(b, a), 30), collapse = ", "))
  }
  invisible(common)
}

safe_sum <- function(mat, cols){
  cols2 <- intersect(cols, colnames(mat))
  if (length(cols2) == 0) return(rep(0, nrow(mat)))
  rowSums(mat[, cols2, drop = FALSE], na.rm = TRUE)
}

infer_metagroup_from_id <- function(sample_id, immune_class = NA_character_){
  sid <- as.character(sample_id)
  suffix <- sub("^.*\\.", "", sid)

  if (grepl("^Dec$",   suffix, ignore.case = TRUE)) return("Decidua")
  if (grepl("^Troph$", suffix, ignore.case = TRUE)) return("Trophoblast")

  if (grepl("\\.sADK$", sid, ignore.case = TRUE)) return("Healthy_EC")
  if (grepl("(^|\\.|_)san($|\\.|_)", sid, ignore.case = TRUE)) return("Healthy_Endometrium")

  if (grepl("\\.tADK$", sid, ignore.case = TRUE)) return("EC")
  NA_character_
}

make_mfi_pairs <- function(meta_df){
  meta_df %>%
    filter(MetaGroup %in% c("Decidua","Trophoblast")) %>%
    mutate(
      isDec   = grepl("^Dec$",   Suffix, ignore.case = TRUE),
      isTroph = grepl("^Troph$", Suffix, ignore.case = TRUE)
    ) %>%
    group_by(PairID) %>%
    summarise(
      DecSample   = SampleID[which(isDec)][1],
      TrophSample = SampleID[which(isTroph)][1],
      .groups = "drop"
    ) %>%
    filter(!is.na(DecSample), !is.na(TrophSample))
}

calc_gene_set_median <- function(expr_mat, genes){
  expr_mat <- as.matrix(expr_mat)
  genes_ok <- intersect(as.character(genes), rownames(expr_mat))
  if (length(genes_ok) < 2) warning("Few genes found in TPM for set: n=", length(genes_ok))
  x <- t(expr_mat[genes_ok, , drop = FALSE])  # samples x genes
  x <- as.data.frame(x)
  x$SampleID <- rownames(x)
  long <- reshape2::melt(x, id.vars = "SampleID", variable.name = "gene", value.name = "value")
  long %>%
    group_by(SampleID) %>%
    summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")
}

get_gene_expr <- function(expr_mat, gene){
  expr_mat <- as.matrix(expr_mat)
  gene_ok <- intersect(gene, rownames(expr_mat))
  if (length(gene_ok) != 1) stop("Gene not found or not unique in TPM rownames: ", gene)
  v <- as.numeric(expr_mat[gene_ok, , drop = TRUE])
  names(v) <- colnames(expr_mat)
  v
}

get_samples_by_metagroup <- function(meta_df, groups){
  meta_df %>% filter(MetaGroup %in% as.character(groups)) %>% pull(SampleID) %>% unique()
}

pick_signal <- function(expr_mat, genes, samples, inter_type){
  genes_ok   <- intersect(as.character(genes), rownames(expr_mat))
  samples_ok <- intersect(as.character(samples), colnames(expr_mat))
  if (length(samples_ok) == 0) return(setNames(numeric(0), character(0)))
  if (length(genes_ok) == 0)   return(setNames(rep(0, length(samples_ok)), samples_ok))
  sub <- expr_mat[genes_ok, samples_ok, drop = FALSE]
  rs  <- rowSums(sub, na.rm = TRUE)
  idx <- if (inter_type == "co-receptors") which.min(rs) else which.max(rs)
  out <- as.numeric(sub[idx, , drop = TRUE])
  names(out) <- colnames(sub)
  out
}

sum_antagonist <- function(expr_mat, genes, samples){
  samples_ok <- intersect(as.character(samples), colnames(expr_mat))
  if (length(samples_ok) == 0) return(setNames(numeric(0), character(0)))
  genes_ok <- intersect(as.character(genes), rownames(expr_mat))
  if (length(genes_ok) == 0) return(setNames(rep(1, length(samples_ok)), samples_ok))
  sub <- expr_mat[genes_ok, samples_ok, drop = FALSE]
  out <- colSums(sub, na.rm = TRUE)
  names(out) <- colnames(sub)
  out
}

## =========================================================
## 3) Load inputs (TPM, counts, CIBERSORT, groups.xlsx)
## =========================================================

## TPM
tpm_raw <- read.csv2(cfg$tpm_file, row.names = 1, check.names = FALSE, sep = "\t")
if ("gene_id" %in% colnames(tpm_raw)) tpm_raw$gene_id <- NULL
tpm_raw[-1] <- lapply(tpm_raw[-1], as.double)
stopifnot("gene_name" %in% colnames(tpm_raw))

tpm_raw <- aggregate(
  tpm_raw[, setdiff(colnames(tpm_raw), "gene_name"), drop = FALSE],
  tpm_raw["gene_name"], sum
)
rownames(tpm_raw) <- tpm_raw$gene_name
tpm_raw$gene_name <- NULL
tpm_raw <- tpm_raw[rowSums(tpm_raw) >= 1, , drop = FALSE]

## COUNTS
counts_raw <- read.csv(cfg$counts_file, row.names = 1, check.names = FALSE, sep = "\t")
if ("gene_id" %in% colnames(counts_raw)) counts_raw$gene_id <- NULL
counts_raw[-1] <- lapply(counts_raw[-1], as.double)
stopifnot("gene_name" %in% colnames(counts_raw))

counts_raw <- aggregate(
  counts_raw[, setdiff(colnames(counts_raw), "gene_name"), drop = FALSE],
  counts_raw["gene_name"], sum
)
rownames(counts_raw) <- counts_raw$gene_name
counts_raw$gene_name <- NULL
counts_raw <- counts_raw[rowSums(counts_raw) >= 1, , drop = FALSE]

## Harmonize matrix column names
colnames(tpm_raw)    <- clean_sample_id(colnames(tpm_raw))
colnames(counts_raw) <- clean_sample_id(colnames(counts_raw))

## CIBERSORT
ciber_raw <- read.csv(cfg$ciber_file, check.names = FALSE)
if ("Mixture" %in% colnames(ciber_raw)) ciber_raw <- ciber_raw %>% dplyr::rename(SampleID = Mixture)
if (!("SampleID" %in% colnames(ciber_raw))) colnames(ciber_raw)[1] <- "SampleID"
ciber_raw$SampleID <- clean_sample_id(ciber_raw$SampleID)

## groups.xlsx (audit)
groups_xlsx <- readxl::read_excel(cfg$groups_xlsx) %>% as.data.frame()
stopifnot(all(c("ID") %in% colnames(groups_xlsx)))

groups_mfp <- groups_xlsx %>%
  transmute(
    SampleID_raw = as.character(ID),
    ImmuneClass  = as.character(dplyr::coalesce(.data$MFP, .data$Cluster, NA_character_)),
    SampleID     = clean_sample_id(as.character(ID))
  ) %>%
  distinct(SampleID, .keep_all = TRUE)

diagnose_overlap(colnames(counts_raw), groups_mfp$SampleID, "counts vs groups.xlsx")

## =========================================================
## 4) Filter to overlap (TPM/Counts)
## =========================================================
common_keep <- Reduce(intersect, list(colnames(counts_raw), colnames(tpm_raw)))
if (length(common_keep) < 2) stop("counts_raw and tpm_raw overlap < 2 samples")

counts_filt <- counts_raw[, common_keep, drop = FALSE]
tpm_filt    <- tpm_raw[,    common_keep, drop = FALSE]
ciber_filt  <- ciber_raw %>% filter(SampleID %in% common_keep)

write.csv(counts_filt, file.path(cfg$out_tables, "counts_filtered.csv"), quote = FALSE)
write.csv(tpm_filt,    file.path(cfg$out_tables, "tpm_filtered.csv"),    quote = FALSE)
write.csv(ciber_filt,  file.path(cfg$out_tables, "cibersort_filtered.csv"), row.names = FALSE, quote = FALSE)

## =========================================================
## 5) Canonical sample metadata (single source of truth)
## =========================================================
sample_meta <- tibble(SampleID = colnames(tpm_filt)) %>%
  mutate(SampleID = as.character(SampleID)) %>%
  left_join(groups_mfp %>% dplyr::select(SampleID, ImmuneClass), by = "SampleID") %>%
  mutate(
    PairID = sub("\\..*$", "", SampleID),
    Suffix = sub("^.*\\.", "", SampleID),
    SuffixToken = vapply(strsplit(SampleID, ".", fixed = TRUE), function(v) tail(v, 1), character(1)),
    IsHyperplasia = SuffixToken %in% cfg$exclude_groups
  ) %>%
  filter(!IsHyperplasia) %>%
  mutate(
    MetaGroup = vapply(seq_len(n()), function(i){
      infer_metagroup_from_id(SampleID[i], ImmuneClass[i])
    }, character(1)),

    ## IMMUNE groups (NOTE: trophoblast excluded from immune labels on purpose)
    Group_immune = case_when(
      MetaGroup == "Healthy_Endometrium" ~ "Healthy - Endometrium Healthy Tissue",
      MetaGroup == "Healthy_EC"          ~ "Healthy - EC Adjacent Healthy Tissue",
      MetaGroup == "Decidua"             ~ "Decidua - MFI",
      MetaGroup == "EC"                  ~ "EC - TIME",
      TRUE ~ NA_character_
    ),

    ## Tissue groups (keep trophoblast for tissue-only plots/DE)
    Group_tissue = case_when(
      MetaGroup == "Healthy_Endometrium" ~ "Healthy_Endometrium",
      MetaGroup == "Healthy_EC"          ~ "Healthy_EC",
      MetaGroup == "Decidua"             ~ "Decidua",
      MetaGroup == "Trophoblast"         ~ "Trophoblast",
      MetaGroup == "EC"                  ~ "EC",
      TRUE ~ NA_character_
    ),

    level_immune = unname(group_levels_immune[Group_immune]),
    level_tissue = unname(group_levels_tissue[Group_tissue])
  ) %>%
  distinct(SampleID, .keep_all = TRUE)

write.csv(sample_meta, file.path(cfg$out_tables, "sample_meta_canonical.csv"),
          row.names = FALSE, quote = FALSE)

message("(MetaGroup distribution)")
print(table(sample_meta$MetaGroup, useNA = "ifany"))

mfi_pairs <- make_mfi_pairs(sample_meta)

## =========================================================
## 6) DESeq2 utilities 
## =========================================================
make_dds <- function(count_mat, meta_df, group_col){
  meta2 <- meta_df %>%
    transmute(SampleID = as.character(SampleID),
              Group = .data[[group_col]]) %>%
    filter(!is.na(Group)) %>%
    distinct(SampleID, .keep_all = TRUE)

  common <- intersect(colnames(count_mat), meta2$SampleID)
  if (length(common) < 2) stop("Too few overlapping samples for model: ", group_col)

  count_mat <- count_mat[, common, drop = FALSE]
  meta2 <- meta2[match(common, meta2$SampleID), , drop = FALSE]
  rownames(meta2) <- meta2$SampleID

  count_mat <- as.matrix(count_mat)
  mode(count_mat) <- "numeric"
  storage.mode(count_mat) <- "integer"

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = data.frame(Group = droplevels(factor(meta2$Group)), row.names = rownames(meta2)),
    design    = ~ Group
  )

  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]
  DESeq(dds)
}

run_contrast_export <- function(dds, contrast_vec, prefix,
                                out_dir = cfg$out_tables,
                                padj_cutoff = cfg$padj_cutoff,
                                lfc_cutoff  = cfg$lfc_cutoff) {
  res <- results(dds, contrast = contrast_vec) |> as.data.frame()
  res <- res[complete.cases(res), , drop = FALSE]
  res_sig <- res %>% filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  out_path <- file.path(out_dir, paste0("DESeq2_", prefix, ".csv"))
  write.csv(res_sig, out_path, quote = FALSE)
  invisible(list(all = res, sig = res_sig, path = out_path))
}

safe_contrast <- function(dds, case, control, prefix){
  lev <- levels(colData(dds)$Group)
  tab <- table(colData(dds)$Group)
  if (!(case %in% lev && control %in% lev)) {
    message("SKIP ", prefix, ": missing level(s). Available: ", paste(lev, collapse = ", "))
    return(NULL)
  }
  if (tab[[case]] < 2 || tab[[control]] < 2) {
    message("WARN ", prefix, ": low n (", case, "=", tab[[case]], ", ", control, "=", tab[[control]], "). Continuing.")
  }
  run_contrast_export(dds, c("Group", case, control), prefix = prefix)
}

metaA <- sample_meta %>%
  transmute(SampleID, MetaGroup) %>%
  mutate(
    Group_DE = case_when(
      MetaGroup == "Healthy_Endometrium" ~ "Healthy_Endometrium",
      MetaGroup == "Healthy_EC"          ~ "Healthy_EC",
      MetaGroup == "Decidua"             ~ "Decidua",
      MetaGroup == "Trophoblast"         ~ "Trophoblast",
      MetaGroup == "EC"                  ~ "EC",
      TRUE ~ NA_character_
    ),
    Group_DE = factor(Group_DE, levels = c("Healthy_Endometrium","Healthy_EC","Decidua","Trophoblast","EC"))
  ) %>%
  filter(!is.na(Group_DE))

ddsA <- make_dds(counts_filt, metaA, "Group_DE")

deA <- list(
  EC_vs_HealthyEC            = safe_contrast(ddsA, "EC", "Healthy_EC",              "EC_vs_HealthyEC"),
  EC_vs_Decidua              = safe_contrast(ddsA, "EC", "Decidua",                 "EC_vs_Decidua"),
  EC_vs_Trophoblast          = safe_contrast(ddsA, "EC", "Trophoblast",             "EC_vs_Trophoblast"),
  Decidua_vs_HealthyEndo     = safe_contrast(ddsA, "Decidua", "Healthy_Endometrium","Decidua_vs_HealthyEndometrium"),
  HealthyEC_vs_HealthyEndo   = safe_contrast(ddsA, "Healthy_EC", "Healthy_Endometrium","HealthyEC_vs_HealthyEndometrium")
)

## =========================================================
## 7) Enrichment heatmaps (reads DE files written above)
## =========================================================
guess_symbol_col <- function(df){
  cand <- c("symbol","SYMBOL","gene","Gene","gene_symbol","GeneSymbol","hgnc_symbol",
            "row.names","Row.names","rownames")
  hit <- cand[cand %in% names(df)][1]
  if (length(hit) == 0 || is.na(hit)) hit <- names(df)[1]
  hit
}

read_de <- function(path){
  df <- read.csv(path)
  symcol <- guess_symbol_col(df)
  df <- df %>%
    dplyr::rename(GENE = dplyr::all_of(symcol)) %>%
    dplyr::mutate(GENE = as.character(GENE),
                  stat = dplyr::coalesce(.data$stat, .data$log2FoldChange)) %>%
    dplyr::filter(!is.na(GENE), nzchar(GENE), !is.na(stat)) %>%
    dplyr::arrange(dplyr::desc(abs(stat))) %>%
    dplyr::distinct(GENE, .keep_all = TRUE)
  df
}

parse_comp <- function(fname){
  base <- tools::file_path_sans_ext(basename(fname))
  m <- stringr::str_match(base, "(.+?)_vs_(.+)")
  if(!is.na(m[1,2])) list(label = base, A = m[1,2], B = m[1,3])
  else list(label = base, A = NA_character_, B = NA_character_)
}

get_msig_list <- function(species, collection, subcollection = NULL){
  mdf <- msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection)
  split(mdf$gene_symbol, mdf$gs_name)
}

run_fgsea_list <- function(ranks_named, gsets, seed){
  set.seed(seed)
  fgsea::fgseaMultilevel(pathways = gsets, stats = ranks_named, scoreType = "std") %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(padj = p.adjust(pval, method = "BH"))
}

sym2entrez <- function(genes){
  suppressMessages({
    map <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  })
  map <- dplyr::filter(map, !is.na(ENTREZID)) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)
  tibble::deframe(map[,c("SYMBOL","ENTREZID")])
}

run_gsekegg <- function(ranks_named, organism = "hsa", seed){
  m <- sym2entrez(names(ranks_named))
  r_entrez <- ranks_named[names(ranks_named) %in% names(m)]
  names(r_entrez) <- unname(m[names(r_entrez)])
  r_entrez <- sort(r_entrez, decreasing = TRUE)
  if (length(r_entrez) < 10) return(tibble(pathway=character(), NES=numeric(), pval=numeric(), padj=numeric()))
  set.seed(seed)
  ek <- clusterProfiler::gseKEGG(geneList = r_entrez, organism = organism,
                                pAdjustMethod = "BH", verbose = FALSE)
  if (is.null(ek) || nrow(as.data.frame(ek)) == 0) {
    tibble(pathway=character(), NES=numeric(), pval=numeric(), padj=numeric())
  } else {
    as.data.frame(ek) %>%
      transmute(pathway = Description, NES = NES, pval = pvalue, padj = p.adjust(pvalue, "BH"))
  }
}

make_mat_and_heatmap <- function(res_df, title, comp_info,
                                 padj_cutoff = 0.05,
                                 top_n_per_comp = 25,
                                 cell_w_mm = 5,
                                 cell_h_mm = 7){

  if (nrow(res_df) == 0) return(NULL)
  sig_df <- res_df %>% dplyr::filter(padj < padj_cutoff)
  if (nrow(sig_df) == 0) return(NULL)

  top_paths <- sig_df %>%
    dplyr::group_by(comparison) %>%
    dplyr::slice_max(order_by = abs(NES), n = top_n_per_comp, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(abs(NES))) %>%
    dplyr::pull(pathway) %>%
    unique()

  if (length(top_paths) == 0) return(NULL)

  nes_wide <- res_df %>%
    dplyr::filter(pathway %in% top_paths) %>%
    dplyr::select(pathway, comparison, NES) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = NES) %>%
    as.data.frame()

  padj_wide <- res_df %>%
    dplyr::filter(pathway %in% top_paths) %>%
    dplyr::select(pathway, comparison, padj) %>%
    tidyr::pivot_wider(names_from = comparison, values_from = padj) %>%
    as.data.frame()

  rownames(nes_wide)  <- nes_wide$pathway
  rownames(padj_wide) <- padj_wide$pathway
  nes_wide$pathway  <- NULL
  padj_wide$pathway <- NULL

  cols_present <- intersect(comp_info$comparison, colnames(nes_wide))
  if (length(cols_present) == 0) return(NULL)

  nes_wide  <- nes_wide[, cols_present, drop = FALSE]
  padj_wide <- padj_wide[, cols_present, drop = FALSE]

  ci_sub <- comp_info %>%
    dplyr::filter(comparison %in% cols_present) %>%
    dplyr::mutate(comparison = factor(comparison, levels = cols_present)) %>%
    dplyr::arrange(comparison)

  colnames(nes_wide)  <- paste0(ci_sub$groupB, " vs ", ci_sub$groupA)
  colnames(padj_wide) <- colnames(nes_wide)

  mat_nes  <- as.matrix(nes_wide);  mode(mat_nes)  <- "numeric"
  mat_padj <- as.matrix(padj_wide); mode(mat_padj) <- "numeric"

  rng <- max(2, stats::quantile(abs(mat_nes), 0.98, na.rm = TRUE))
  col_fun <- circlize::colorRamp2(c(-rng, 0, rng), c("#3B73B9", "white", "#B83C3D"))

  star_for_p <- function(p){
    if (is.na(p)) ""
    else if (p < 0.001) "***"
    else if (p < 0.01)  "**"
    else if (p < padj_cutoff) "*"
    else ""
  }

  cell_star <- function(j, i, x, y, w, h, fill){
    s <- star_for_p(mat_padj[i, j])
    if (nzchar(s)) grid::grid.text(s, x = x, y = y, gp = grid::gpar(fontsize = 9, fontface = "bold"))
  }

  ht <- ComplexHeatmap::Heatmap(
    mat_nes,
    name = "NES",
    col = col_fun,
    na_col = "grey90",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_title = title,
    width  = grid::unit(cell_w_mm, "mm") * ncol(mat_nes),
    height = grid::unit(cell_h_mm, "mm") * nrow(mat_nes),
    row_names_gp = grid::gpar(fontsize = 9),
    heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topcenter"),
    cell_fun = cell_star
  )

  list(ht = ht, mat_nes = nes_wide, mat_padj = padj_wide, comp_info = ci_sub)
}

files <- list.files(cfg$enrichment$excel_dir, pattern = cfg$enrichment$file_pattern, full.names = TRUE)
res_h <- res_g <- res_k <- tibble()
comp_info <- tibble(comparison=character(), groupA=character(), groupB=character())

if (length(files) > 0) {
  hallmark_sets <- get_msig_list(cfg$enrichment$species, "H")
  go_sets       <- get_msig_list(cfg$enrichment$species, "C5", cfg$enrichment$go_subcategory)

  do_one <- function(path){
    de    <- read_de(path)
    ranks <- tibble::deframe(de %>% dplyr::select(GENE, stat))
    ranks <- sort(ranks, decreasing = TRUE)

    comp  <- parse_comp(path)
    label <- comp$label

    hallmark <- run_fgsea_list(ranks, hallmark_sets, cfg$enrichment$seed) %>%
      transmute(pathway, NES, pval, padj, comparison = label, ont = "Hallmark")

    gores <- run_fgsea_list(ranks, go_sets, cfg$enrichment$seed) %>%
      transmute(pathway, NES, pval, padj, comparison = label, ont = paste0("GO_", cfg$enrichment$go_subcategory))

    kegg <- run_gsekegg(ranks, organism = "hsa", seed = cfg$enrichment$seed) %>%
      mutate(comparison = label, ont = "KEGG")

    list(
      comp_info = tibble(comparison = label, groupA = comp$A, groupB = comp$B),
      hallmark  = hallmark,
      go        = gores,
      kegg      = kegg
    )
  }

  all_results <- purrr::map(files, do_one)
  comp_info <- bind_rows(purrr::map(all_results, "comp_info")) %>% arrange(comparison)
  res_h <- bind_rows(purrr::map(all_results, "hallmark"))
  res_g <- bind_rows(purrr::map(all_results, "go"))
  res_k <- bind_rows(purrr::map(all_results, "kegg"))

  out_pdf <- file.path(cfg$out_figs, "enrichment_heatmaps.pdf")
  pdf(out_pdf, width = 12, height = 20)

  out_h <- make_mat_and_heatmap(res_h, "Hallmark", comp_info, cfg$padj_cutoff, cfg$enrichment$top_n_per_comp, 5, 7)
  if (!is.null(out_h)) ComplexHeatmap::draw(out_h$ht, heatmap_legend_side = "bottom")

  out_k <- make_mat_and_heatmap(res_k, "KEGG", comp_info, cfg$padj_cutoff, cfg$enrichment$top_n_per_comp, 5, 7)
  if (!is.null(out_k)) ComplexHeatmap::draw(out_k$ht, heatmap_legend_side = "bottom")

  out_g <- make_mat_and_heatmap(res_g, paste0("GO ", cfg$enrichment$go_subcategory), comp_info,
                                cfg$padj_cutoff, cfg$enrichment$top_n_per_comp, 5, 7)
  if (!is.null(out_g)) ComplexHeatmap::draw(out_g$ht, heatmap_legend_side = "bottom")

  dev.off()
  message("Saved enrichment heatmaps: ", out_pdf)
} else {
  message("No DE files found for enrichment under: ", cfg$enrichment$excel_dir)
}

## =========================================================
## 8) CIBERSORT processing 
## =========================================================
drop_cols <- c("RMSE","Correlation")
ciber2 <- ciber_filt
ciber2 <- ciber2[, setdiff(colnames(ciber2), drop_cols), drop = FALSE]

if ("P.value" %in% colnames(ciber2)) {
  ciber2 <- ciber2 %>% filter(P.value < cfg$ciber_p_cutoff)
}

rownames(ciber2) <- ciber2$SampleID
ciber2$SampleID <- NULL
ciber2[,] <- apply(ciber2, 2, function(x) as.numeric(as.character(x)))

Bcells      <- c("B cells naive","B cells memory","Plasma cells")
TcellsCD4   <- c("T cells CD4 naive","T cells CD4 memory resting","T cells CD4 memory activated","T cells regulatory (Tregs)")
Macrophages <- c("Macrophages M0","Macrophages M1","Macrophages M2")
NKs         <- c("NK cells resting","NK cells activated")
DCs         <- c("Dendritic cells resting","Dendritic cells activated")

ciber2$B.cells     <- safe_sum(ciber2, Bcells)
ciber2$T.cells.CD4 <- safe_sum(ciber2, TcellsCD4)
ciber2$Macrophages <- safe_sum(ciber2, Macrophages)
ciber2$NK          <- safe_sum(ciber2, NKs)
ciber2$Dendritic   <- safe_sum(ciber2, DCs)

write.csv(ciber2, file.path(cfg$out_tables, "ciber2_processed.csv"), quote = FALSE)

## =========================================================
## 10) ICI complex scoring (MFI paired: Decidua immune; Troph APC only)
## =========================================================
tpm_mat <- as.matrix(tpm_filt)
mode(tpm_mat) <- "numeric"

complex_defs <- list(
  list(cname="CD112 - CD112R",              module="co-receptors",
       APC=c("PVRL2"), APC_ant=c(), ICI=c("PVRIG"), ICI_ant=c()),
  list(cname="HLA-G/HLA-F - LILRB1/LILRB2", module="alt-receptors",
       APC=c("HLA-F","HLA-G"), APC_ant=c(), ICI=c("LILRB1","LILRB2"), ICI_ant=c()),
  list(cname="HLA-G - KIR2DL4",             module="co-receptors",
       APC=c("HLA-G"), APC_ant=c("KIR2DL4"), ICI=c(), ICI_ant=c("KIR2DL4")),
  list(cname="CD47 - SIRPG/SIRPA",          module="alt-receptors",
       APC=c("CD47"), APC_ant=c(), ICI=c("SIRPA"), ICI_ant=c()),
  list(cname="HLA-E - CD94",                module="co-receptors",
       APC=c("HLA-E"), APC_ant=c(), ICI=c("KLRD1","KLRC1"), ICI_ant=c("KLRK1")),
  list(cname="LGALS9.CEACAM1 - TIM3",       module="co-receptors",
       APC=c("LGALS9","CEACAM1"), APC_ant=c(), ICI=c("HAVCR2"), ICI_ant=c()),
  list(cname="CD112/CD155 - TIGIT",         module="alt-receptors",
       APC=c("PVR","PVRL2"), APC_ant=c(), ICI=c("TIGIT"), ICI_ant=c("CD226")),
  list(cname="B7 - CTLA4",                  module="alt-receptors",
       APC=c("CD80","CD86"), APC_ant=c(), ICI=c("CTLA4"), ICI_ant=c("CD28")),
  list(cname="HLAII.GAL3 - LAG3",           module="alt-receptors",
       APC=c("HLA-DP","HLA-DM","HLA-DOA","HLA-DOB","HLA-DQ","LGALS3"),
       APC_ant=c(), ICI=c("LAG3"), ICI_ant=c("CD28"))
)

cohorts <- list(
  MFI_pair = list(group_name = "Decidua - MFI (paired APC from trophoblast)", groupA = "Decidua", groupB = "Trophoblast"),
  EC_TIME  = list(group_name = "EC - TIME",  groupA = "EC",      groupB = "EC"),
  HealthyEndometrium = list(group_name = "Healthy - Endometrium Healthy Tissue", groupA = "Healthy_Endometrium", groupB = "Healthy_Endometrium"),
  HealthyEC          = list(group_name = "Healthy - EC Adjacent Healthy Tissue", groupA = "Healthy_EC",          groupB = "Healthy_EC"),
  DeciduaMFI         = list(group_name = "Decidua - MFI",        groupA = "Decidua", groupB = "Decidua")
)

score_complex <- function(expr_mat, meta_df, cohort_key, cohort, complex_def, mfi_pairs){

  if (cohort_key == "MFI_pair") {
    if (nrow(mfi_pairs) == 0) return(NULL)

    dec_ids   <- intersect(mfi_pairs$DecSample, colnames(expr_mat))
    troph_ids <- intersect(mfi_pairs$TrophSample, colnames(expr_mat))
    if (length(dec_ids) < 1 || length(troph_ids) < 1) return(NULL)

    ## order by PairID as in mfi_pairs
    dec_ids   <- mfi_pairs$DecSample[mfi_pairs$DecSample %in% dec_ids]
    troph_ids <- mfi_pairs$TrophSample[mfi_pairs$TrophSample %in% troph_ids]

    ## Decidua supplies immune-side genes
    ICI_sig <- pick_signal(expr_mat, complex_def$ICI, dec_ids, complex_def$module)
    ICI_ant <- sum_antagonist(expr_mat, complex_def$ICI_ant, dec_ids)

    ## Trophoblast supplies APC-side genes (mapped onto decidua sample IDs)
    APC_sig <- pick_signal(expr_mat, complex_def$APC, troph_ids, complex_def$module)
    APC_ant <- sum_antagonist(expr_mat, complex_def$APC_ant, troph_ids)

    troph_to_dec <- setNames(mfi_pairs$DecSample, mfi_pairs$TrophSample)

    APC_sig2 <- APC_sig[names(APC_sig) %in% names(troph_to_dec)]
    APC_ant2 <- APC_ant[names(APC_ant) %in% names(troph_to_dec)]
    names(APC_sig2) <- troph_to_dec[names(APC_sig2)]
    names(APC_ant2) <- troph_to_dec[names(APC_ant2)]

    dec_ids2 <- dec_ids
    if (length(ICI_sig) == 0) ICI_sig <- setNames(rep(0, length(dec_ids2)), dec_ids2) else ICI_sig <- ICI_sig[dec_ids2]
    if (length(ICI_ant) == 0) ICI_ant <- setNames(rep(1, length(dec_ids2)), dec_ids2) else ICI_ant <- ICI_ant[dec_ids2]

    APC_sig_aligned <- if (all(dec_ids2 %in% names(APC_sig2))) APC_sig2[dec_ids2] else {
      mu <- if (length(APC_sig2) > 0) mean(APC_sig2, na.rm=TRUE) else 0
      setNames(rep(mu, length(dec_ids2)), dec_ids2)
    }

    APC_den_aligned <- if (all(dec_ids2 %in% names(APC_ant2))) APC_ant2[dec_ids2] else {
      mu <- if (length(APC_ant2) > 0) mean(APC_ant2, na.rm=TRUE) else 1
      setNames(rep(mu, length(dec_ids2)), dec_ids2)
    }

    sig_value <- (APC_sig_aligned / (APC_den_aligned + 0.001)) *
      (pmax(1, ICI_sig) / (ICI_ant + 0.001))

    return(tibble(
      SampleID  = dec_ids2,
      sig_value = as.numeric(sig_value),
      GroupName = "Decidua - MFI",   ## force integration into decidua immune group
      Complex   = complex_def$cname,
      InterType = complex_def$module,
      CohortKey = cohort_key
    ))
  }

  ## non-MFI cohorts: within-sample A-side, B-side for APC antagonists if needed
  sA <- get_samples_by_metagroup(meta_df, cohort$groupA)
  sB <- get_samples_by_metagroup(meta_df, cohort$groupB)
  if (length(sA) < 1) return(NULL)

  ICI_sig <- pick_signal(expr_mat, complex_def$ICI, sA, complex_def$module)
  APC_sig <- pick_signal(expr_mat, complex_def$APC, sA, complex_def$module)

  ICI_ant <- sum_antagonist(expr_mat, complex_def$ICI_ant, sA)
  APC_ant <- sum_antagonist(expr_mat, complex_def$APC_ant, sB)

  if (length(ICI_sig) == 0) ICI_sig <- setNames(rep(0, length(sA)), sA) else ICI_sig <- ICI_sig[sA]
  if (length(APC_sig) == 0) APC_sig <- setNames(rep(0, length(sA)), sA) else APC_sig <- APC_sig[sA]
  if (length(ICI_ant) == 0) ICI_ant <- setNames(rep(1, length(sA)), sA) else ICI_ant <- ICI_ant[sA]

  APC_den <- if (all(sA %in% names(APC_ant))) APC_ant[sA]
  else if (length(APC_ant) > 0) setNames(rep(mean(APC_ant, na.rm=TRUE), length(sA)), sA)
  else setNames(rep(1, length(sA)), sA)

  sig_value <- (APC_sig / (APC_den + 0.001)) * (pmax(1, ICI_sig) / (ICI_ant + 0.001))

  tibble(
    SampleID  = sA,
    sig_value = as.numeric(sig_value),
    GroupName = cohort$group_name,
    Complex   = complex_def$cname,
    InterType = complex_def$module,
    CohortKey = cohort_key
  )
}

ici_scores <- purrr::map_dfr(names(cohorts), function(ck){
  cohort <- cohorts[[ck]]
  purrr::map_dfr(complex_defs, function(cx){
    score_complex(tpm_mat, sample_meta, ck, cohort, cx, mfi_pairs)
  })
})

ici_scores <- ici_scores %>%
  filter(CohortKey %in% c("MFI_pair","EC_TIME","HealthyEndometrium","HealthyEC","DeciduaMFI"))

write.csv(ici_scores, file.path(cfg$out_tables, "ICI_complex_scores.csv"),
          row.names = FALSE, quote = FALSE)

## =========================================================
## 11) df_all (ICI + sample_meta + CIBERSORT)
## =========================================================
ciber_df <- as.data.frame(ciber2)
ciber_df$SampleID <- rownames(ciber_df)
ciber_df$SampleID <- as.character(ciber_df$SampleID)
rownames(ciber_df) <- NULL

df_all <- ici_scores %>%
  mutate(SampleID = as.character(SampleID)) %>%
  left_join(sample_meta, by = "SampleID") %>%
  left_join(distinct(ciber_df, SampleID, .keep_all = TRUE), by = "SampleID") %>%
  mutate(sig_value_log2 = log2(sig_value + 1e-6)) %>%
  ## immune plots exclude trophoblast entirely
  filter(MetaGroup != "Trophoblast")

write.csv(df_all, file.path(cfg$out_tables, "df_all_merged.csv"),
          row.names = FALSE, quote = FALSE)

## =========================================================
## 12) Plot ICI panel (immune view; troph excluded already)
## =========================================================
plot_ici_panel <- function(df, complexes, outfile, title){

  df_sub <- df %>%
    filter(!is.na(Group_immune), Complex %in% complexes) %>%
    mutate(
      Group_immune = factor(Group_immune, levels = names(group_levels_immune)),
      Complex = factor(Complex, levels = complexes)
    )

  if (nrow(df_sub) == 0) {
    message("SKIP: ", outfile, " (no data)")
    return(invisible(NULL))
  }

  candidate_comparisons <- list(
    c("EC - TIME", "Healthy - EC Adjacent Healthy Tissue"),
    c("EC - TIME", "Decidua - MFI"),
    c("Decidua - MFI", "Healthy - Endometrium Healthy Tissue")
  )
  valid_groups <- levels(df_sub$Group_immune)
  comparisons <- Filter(function(x) all(x %in% valid_groups), candidate_comparisons)

  p <- ggplot(df_sub, aes(x = Group_immune, y = sig_value_log2, color = Group_immune)) +
    geom_boxplot(aes(fill = Group_immune), color="black", outlier.shape=NA, alpha=0.7, linewidth=0.35) +
    geom_jitter(width=0.18, size=1.2, alpha=0.55) +
    facet_wrap(~ Complex, ncol = 5, scales = "free_y") +
    scale_color_manual(values = group_colors_immune, drop = FALSE) +
    scale_fill_manual(values = group_colors_immune, drop = FALSE) +
    labs(x = NULL, y = "ICI complex interaction score (log2)", color = "Group", fill = "Group") +
    ggtitle(title) +
    {
      if (length(comparisons) > 0) {
        ggpubr::stat_compare_means(
          method = "wilcox",
          comparisons = comparisons,
          label = "p.signif",
          paired = FALSE,
          hide.ns = TRUE,
          vjust = 0.6
        )
      } else NULL
    } +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.text.position ='right',
      legend.title.position = 'top',
      axis.text.x = element_blank(),
      strip.text = element_text(size=10, face="bold")
    )

  png(outfile, width=4000, height=2000, res=300)
  print(p)
  dev.off()
  message("Saved: ", outfile)
  invisible(p)
}

panel <- c("LGALS9.CEACAM1 - TIM3","HLA-E - CD94","CD112 - CD112R",
           "CD112/CD155 - TIGIT","CD47 - SIRPG/SIRPA","B7 - CTLA4",
           "HLAII.GAL3 - LAG3","HLA-G - KIR2DL4","HLA-G/HLA-F - LILRB1/LILRB2")

plot_ici_panel(df_all, panel, file.path(cfg$out_figs, "ICI_panel.png"), "ICI complexes")

## =========================================================
## 13) CIBERSORT “immune state” barplot (IQR error bars)
## Macrophages: polarized (M1+M2) vs not polarized (M0)
## =========================================================
abund_parts <- df_all %>%
  distinct(
    SampleID, Group_immune,
    `T cells CD4 memory resting`, `T cells CD4 memory activated`,
    `NK cells resting`, `NK cells activated`,
    `Dendritic cells resting`, `Dendritic cells activated`,
    `Macrophages M0`, `Macrophages M1`, `Macrophages M2`
  ) %>%
  filter(!is.na(Group_immune)) %>%
  mutate(
    CD4mem_inactive = `T cells CD4 memory resting`,
    CD4mem_active   = `T cells CD4 memory activated`,
    NK_inactive     = `NK cells resting`,
    NK_active       = `NK cells activated`,
    DC_inactive     = `Dendritic cells resting`,
    DC_active       = `Dendritic cells activated`,
    Mac_np          = `Macrophages M0`,
    Mac_pol         = `Macrophages M1` + `Macrophages M2`
  )

df_long_bar <- abund_parts %>%
  select(
    SampleID, Group_immune,
    CD4mem_inactive, CD4mem_active,
    NK_inactive, NK_active,
    DC_inactive, DC_active,
    Mac_np, Mac_pol
  ) %>%
  pivot_longer(cols = -c(SampleID, Group_immune),
               names_to = "feature", values_to = "Value") %>%
  separate(feature, into = c("Cell_Type", "state"), sep = "_") %>%
  mutate(
    Cell_Type = recode(Cell_Type,
      CD4mem = "CD4 memory",
      NK     = "NK",
      DC     = "Dendritic",
      Mac    = "Macrophages (polarization)"
    ),
    state = case_when(
      Cell_Type == "Macrophages (polarization)" & state == "np"  ~ "Not polarized (M0)",
      Cell_Type == "Macrophages (polarization)" & state == "pol" ~ "Polarized (M1+M2)",
      state == "inactive" ~ "Inactive",
      state == "active"   ~ "Active",
      TRUE ~ state
    ),
    state = factor(state, levels = c("Inactive","Active","Not polarized (M0)","Polarized (M1+M2)")),
    Group_immune = factor(Group_immune, levels = names(group_levels_immune))
  )

bar_sum <- df_long_bar %>%
  group_by(Group_immune, Cell_Type, state) %>%
  summarise(
    median_value = median(Value, na.rm=TRUE),
    q25 = quantile(Value,0.25,na.rm=TRUE),
    q75 = quantile(Value,0.75,na.rm=TRUE),
    .groups="drop"
  )

p_bar_err <- ggplot(bar_sum,
                    aes(x=Group_immune, y=median_value, fill=Group_immune, alpha=state)) +
  geom_col(position="stack", color="black", linewidth=0.2, width=0.65) +
  geom_errorbar(aes(ymin=q25, ymax=q75),
                position=position_stack(vjust=0.5),
                width=0.18, linewidth=0.3) +
  facet_wrap(~Cell_Type, ncol=2, scales="free_y") +
  scale_fill_manual(values=group_colors_immune, drop=FALSE) +
  scale_alpha_manual(values=c("Inactive"=0.3,"Active"=1,"Not polarized (M0)"=0.3,"Polarized (M1+M2)"=1)) +
  labs(x=NULL, y="Median relative abundance (IQR bars)", fill="Group", alpha="State") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90,hjust=1),
        legend.position="bottom",
        legend.direction="vertical")

ggsave(file.path(cfg$out_figs,"stacked_median_active_inactive_IQR.png"),
       p_bar_err, dpi=300, width=5, height=8, units="in")

## =========================================================
## 14) Inactivation ratios (computed per sample)
## =========================================================
df_corr <- df_all %>%
  distinct(SampleID, Group_immune, Complex, sig_value_log2,
           `T cells CD4 memory resting`, `T cells CD4 memory activated`,
           `NK cells resting`, `NK cells activated`,
           `Dendritic cells resting`, `Dendritic cells activated`,
           `Macrophages M1`, `Macrophages M2`) %>%
  filter(!is.na(Group_immune), !is.na(sig_value_log2)) %>%
  mutate(
    CD4_memory_inactivation =
      `T cells CD4 memory resting` / (`T cells CD4 memory activated` + `T cells CD4 memory resting` + 0.001),
    NK_inactivation =
      `NK cells resting` / (`NK cells activated` + `NK cells resting` + 0.001),
    DC_inactivation =
      `Dendritic cells resting` / (`Dendritic cells activated` + `Dendritic cells resting` + 0.001),
    Macrophages_inactivation =
      `Macrophages M2` / (`Macrophages M1` + `Macrophages M2` + 0.001),
    Group_immune = factor(Group_immune, levels = names(group_levels_immune))
  )

## =========================================================
## 15) TPM signatures: Phagocytosis + CD8 exhaustion
## =========================================================
genes_negative_reg_fagocitosi  <- c(
  "RACK1","CD300A","CNN2","CSK","CD300LF","FCGR2B",
  "SYT11","APPL1","HMGB1","MIR181B1","PLSCR1","PIP4P2","PRTN3",
  "ATG3","TGFB1","TLR2","DYSF","SNX3","ADIPOQ","ATG5")

genes_positive_reg_fagocitosi <- c(
  "C2","C3","CALR","CAMK1D","CCL2","CCR7","CD36","CD47","CD209B","CD300LF",
  "CFP","CLEC7A","COLEC10","COLEC11","CYBA","DNM2","DOCK2","F2RL1","FCER1G",
  "FCGR1","FCGR2B","FCGR3","FCNB","FPR2","GAS6","GATA2","HSPA8","IFNG",
  "IL15","IL15RA","IL2RB","IL2RG","ITGA2","LBP","LRP1","MBL1","MBL2","MERTK",
  "MFGE8","MYH9","NCKAP1L","NOD2","PLA2G5","PLCG2","PPARG","PROS1","PTPRC",
  "PTPRJ","PTK2","PTX3","PYCARD","RAB27A","RAP1A","RAPGEF1","SFTPA1","SFTPD",
  "SIRPA","SLC11A1","SOD1","TREM2"
)

neg_reg_fagocitosi <- calc_gene_set_median(tpm_mat, genes_negative_reg_fagocitosi) %>%
  rename(median_value_neg = median_value)
pos_reg_fagocitosi <- calc_gene_set_median(tpm_mat, genes_positive_reg_fagocitosi) %>%
  rename(median_value_pos = median_value)

phago_df <- neg_reg_fagocitosi %>%
  left_join(pos_reg_fagocitosi, by="SampleID") %>%
  mutate(phago_ratio = median_value_neg / (median_value_pos + median_value_neg +1e-6))


genes_up_cd8_exhaustion <- c("PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","CD244","ENTPD1",
              "TOX","NR4A1","NR4A2","NR4A3","BATF","IRF4","PRDM1","EOMES","CXCL13")
genes_down_cd8_exhaustion <- c("TCF7","IL7R","CCR7","SELL","LEF1","CD28","IFNG","IL2","TNF","GZMB","GNLY")

compute_CD8_TI <- function(expr_mat){
  up_present   <- intersect(genes_up_cd8_exhaustion, rownames(expr_mat))
  down_present <- intersect(genes_down_cd8_exhaustion, rownames(expr_mat))
  up_score   <- colMeans(expr_mat[genes_up_cd8_exhaustion, , drop = FALSE], na.rm = TRUE)
  down_score <- colMeans(expr_mat[genes_down_cd8_exhaustion, , drop = FALSE], na.rm = TRUE)
  up_score / (down_score + up_score + 1e-6)
}

cd8_ti <- compute_CD8_TI(tpm_mat)
cd8_df <- tibble(SampleID = names(cd8_ti), cd8_exhaustion = as.numeric(cd8_ti))

sig_per_sample <- sample_meta %>%
  select(SampleID, Group_immune, Group_tissue) %>%
  left_join(cd8_df, by="SampleID") %>%
  left_join(phago_df %>% select(SampleID, median_value_neg, phago_ratio), by="SampleID") %>%
  left_join(df_corr %>% distinct(SampleID, NK_inactivation), by="SampleID") %>%
  filter(!is.na(Group_immune)) %>%                         ## immune plots only
  mutate(Group_immune = factor(Group_immune, levels=names(group_levels_immune)))

write.csv(sig_per_sample, file.path(cfg$out_tables, "sig_per_sample.csv"),
          row.names=FALSE, quote=FALSE)

## =========================================================
## 16) Violin plots: add BH-adjusted significance stars
## =========================================================

make_pairwise_stats <- function(df, group_col, value_col, comparisons, p_adjust="BH"){

  df2 <- df %>%
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[value_col]])) %>%
    dplyr::mutate(.grp = factor(.data[[group_col]], levels = levels(.data[[group_col]]))) %>%
    droplevels()

  lv <- levels(df2$.grp)
  comps_ok <- Filter(function(x) all(x %in% lv), comparisons)
  if (length(comps_ok) == 0) return(NULL)

  st <- ggpubr::compare_means(
    formula = as.formula(paste0(value_col, " ~ .grp")),
    data = df2,
    method = "wilcox.test",
    comparisons = comps_ok,
    p.adjust.method = p_adjust
  )

  ## ggpubr returns p.adj reliably if p.adjust.method is set
  if (!("p.adj" %in% colnames(st))) {
    ## fallback (should rarely happen)
    st$p.adj <- p.adjust(st$p, method = p_adjust)
  }

  star <- function(p){
    ifelse(is.na(p), "",
           ifelse(p <= 0.001, "***",
                  ifelse(p <= 0.01, "**",
                         ifelse(p <= 0.05, "*", ""))))
  }

  ## y positions: above max (in original scale)
  ymax <- max(df2[[value_col]], na.rm = TRUE)
  step <- 0.12 * ymax + 1e-6

  st %>%
    dplyr::mutate(
      y.position = ymax + seq_len(n()) * step,
      label = star(p.adj),
      group1 = as.character(group1),
      group2 = as.character(group2)
    )
}

comparisons_immune <- list(
  c("EC - TIME", "Healthy - EC Adjacent Healthy Tissue"),
  c("EC - TIME", "Decidua - MFI"),
  c("Decidua - MFI", "Healthy - Endometrium Healthy Tissue")
)

st_cd8 <- make_pairwise_stats(sig_per_sample, "Group_immune", "cd8_exhaustion", comparisons_immune)
p_cd8 <- ggplot(sig_per_sample, aes(x=Group_immune, y=cd8_exhaustion, fill=Group_immune)) +
  geom_violin(trim=FALSE, alpha=0.75, linewidth=0.25) +
  geom_boxplot(width=0.12, outlier.shape=NA, color = 'white',linewidth=0.25) +
  geom_jitter(width=0.15, size=0.8, alpha=0.4) +
  scale_fill_manual(values=group_colors_immune, drop=FALSE) +
  scale_y_continuous(trans="log2") +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_blank(),legend.position = 'none',
        axis.title.y.left = element_text('Score'),
        axis.ticks.y = element_blank(),
        axis.text.y =  element_blank()) +
  coord_cartesian(clip="off")


if(!is.null(st_cd8)) p_cd8 <- p_cd8 + ggpubr::stat_pvalue_manual(st_cd8, label="label",
                                                                 y.position="y.position",
                                                                 tip.length=0.01, size=4)

ggsave(file.path(cfg$out_figs,"CD8_exhaustion_violin_BH.png"),
       p_cd8, dpi=300, width=130, height=70, units="mm")

st_phago <- make_pairwise_stats(sig_per_sample, "Group_immune", "phago_ratio", comparisons_immune)
p_phago <- ggplot(sig_per_sample, aes(x=Group_immune, y=phago_ratio, fill=Group_immune)) +
  geom_violin(trim=FALSE, alpha=0.75, color="black", linewidth=0.25) +
  geom_boxplot(width=0.12, outlier.shape=NA, color = 'white',linewidth=0.25) +
  geom_jitter(width=0.15, size=0.8, alpha=0.4) +
  scale_fill_manual(values=group_colors_immune, drop=FALSE) +
  scale_y_continuous(trans="log2") +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_blank(),legend.position = 'none',
        axis.title.y.left = element_text('Score'),
        axis.ticks.y = element_blank(),
        axis.text.y =  element_blank()) +
  coord_cartesian(clip="off") 

if(!is.null(st_phago)) p_phago <- p_phago + ggpubr::stat_pvalue_manual(st_phago, label="label",
                                                                       y.position="y.position",
                                                                       tip.length=0.01, size=4)

ggsave(file.path(cfg$out_figs,"Phagocytosis_violin_BH.png"),
       p_phago, dpi=300, width=130, height=70, units="mm")

## =========================================================
## 17) Correlation panels: x=ICI (log2), y=signature
## BH adjust p-values across number of ICI used per signature
## =========================================================
ici_keep <- c(
  "CD112 - CD112R",
  "CD112/CD155 - TIGIT",
  "CD47 - SIRPG/SIRPA",
  "LGALS9.CEACAM1 - TIM3",
  "HLA-G - KIR2DL4",
  "HLA-G/HLA-F - LILRB1/LILRB2"
)

ici_for_corr <- df_all %>%
  select(SampleID, Complex, sig_value_log2, Group_immune) %>%
  distinct() %>%
  filter(Complex %in% ici_keep, !is.na(Group_immune))

corr_long <- ici_for_corr %>%
  left_join(sig_per_sample, by=c("SampleID","Group_immune")) %>%
  pivot_longer(
    cols = c(cd8_exhaustion, median_value_neg, NK_inactivation, phago_ratio),
    names_to = "Signature",
    values_to = "SignatureValue"
  ) %>%
  filter(!is.na(SignatureValue), !is.na(sig_value_log2), !is.na(Complex)) %>%
  mutate(
    Signature = recode(Signature,
      cd8_exhaustion   = "CD8 exhaustion",
      median_value_neg = "Phagocytosis NEG regulation (median TPM)",
      NK_inactivation  = "NK inactivation",
      phago_ratio      = "Phagocytosis ratio (neg/(pos+neg))"
    ),
    Group_immune = factor(Group_immune, levels = names(group_levels_immune))
  )

## compute spearman per facet + BH adjust across complexes per signature
facet_stats <- corr_long %>%
  group_by(Signature, Complex) %>%
  summarise(
    rho = suppressWarnings(cor(sig_value_log2, SignatureValue, method="spearman", use="pairwise.complete.obs")),
    p   = suppressWarnings(cor.test(sig_value_log2, SignatureValue, method="spearman")$p.value),
    .groups="drop"
  ) %>%
  group_by(Signature) %>%
  mutate(p_adj = p.adjust(p, method="BH")) %>%
  ungroup() %>%
  mutate(
    lab = paste0("rho=", sprintf("%.2f", rho), "\nBH p=", format.pval(p_adj, digits=2, eps=1e-3))
  )

## label positions per facet (top-left)
label_pos <- corr_long %>%
  group_by(Signature, Complex) %>%
  summarise(
    x = quantile(sig_value_log2, 0.05, na.rm=TRUE),
    y = quantile(SignatureValue, 0.95, na.rm=TRUE),
    .groups="drop"
  ) %>%
  left_join(facet_stats, by=c("Signature","Complex"))

p_corr_panel <- ggplot(corr_long, aes(x = sig_value_log2, y = SignatureValue, color = Group_immune)) +
  geom_point(alpha=0.5, size=1.5) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.7) +
  facet_grid(Signature ~ Complex, scales="free") +
  geom_text(data = label_pos, aes(x=x, y=y, label=lab),
            inherit.aes = FALSE, hjust=0, vjust=1, size=2.6) +
  scale_color_manual(values = group_colors_immune, drop=FALSE) +
  labs(x="ICI score (log2)", y=NULL, color="Group",
       title="Correlations: ICI scores (x) vs signatures (y) — BH adjusted across ICI per signature") +
  theme_classic() +
  theme(legend.position="bottom",
        strip.text.x = element_text(size=9, face="bold"),
        strip.text.y = element_text(size=9, face="bold"))

ggsave(file.path(cfg$out_figs, "corr_panels_ICI_vs_signatures_BH.png"),
       p_corr_panel, dpi=350, width=14, height=9, units="in")

## =========================================================
## 18) DE heatmap (alignment fix included)
## =========================================================
make_de_heatmap <- function(dds, sample_meta, group_col="Group_tissue", outfile, top_n=100){

  vsd <- DESeq2::vst(dds, blind=TRUE)
  mat <- SummarizedExperiment::assay(vsd)

  smp <- sample_meta %>%
    select(SampleID, all_of(group_col)) %>%
    filter(!is.na(.data[[group_col]])) %>%
    distinct(SampleID, .keep_all = TRUE)

  common <- intersect(colnames(mat), smp$SampleID)
  if (length(common) < 2) stop("Too few overlapping samples between dds and sample_meta")

  mat <- mat[, common, drop = FALSE]
  smp <- smp[match(common, smp$SampleID), , drop = FALSE]
  stopifnot(all(smp$SampleID == colnames(mat)))

  rv <- matrixStats::rowVars(mat)
  top <- names(sort(rv, decreasing=TRUE))[seq_len(min(top_n, length(rv)))]
  mat_top <- mat[top, , drop=FALSE]

  grp <- smp[[group_col]]
  names(grp) <- smp$SampleID

  ha <- ComplexHeatmap::HeatmapAnnotation(
    Group = grp,
    col = list(Group = group_colors_tissue),
    which = "column"
  )

  ht <- ComplexHeatmap::Heatmap(
    mat_top, name="vst",
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE
  )

  png(outfile, width=2200, height=1800, res=300)
  ComplexHeatmap::draw(ht)
  dev.off()
  invisible(ht)
}

make_de_heatmap(ddsA, sample_meta,
                group_col="Group_tissue",
                outfile=file.path(cfg$out_figs,"DE_heatmap_ddsA.png"),
                top_n=100)

message("PIPELINE COMPLETED. Outputs written to: ", cfg$out_dir)
