
############################################################
## NETWORK ANALYSIS PIPELINE (single chunk, robust)
## - CLI: --groups --tpm --counts [--out_prefix --overwrite]
## - Builds: TPM/count matrices, sample_meta, group matrices
## - Defines: make_diff_graph(), run_or_load_diff_graph(),
##            analyze_one_graph(), compare_networks_plots()
## - Fixes:
##   * cache signature uses digest (no integer overflow, no NA)
##   * LCC selection robust to igraph 2.x membership types
##   * uses delete_vertices() (igraph 2.x)
##   * module plots use layout_with_fr() (no layout_with_dh weights arg)
##   * compare_networks_plots hardened for empty graphs
## - Adds:
##   * oncology-oriented module ranking + ModuleType classification
##   * module meta-graph plot (modules as nodes, inter-module edges)
##   * GO:BP enrichment + Hallmark (MSigDB H) enrichment via msigdbr
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readxl)
  library(matrixStats)
})


## =========================================================
## Helpers
## =========================================================
clean_sample_id <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x
}

read_salmon_merged <- function(path){
  # expects: rownames = feature id; one column named gene_name; remaining = samples
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
    warning("Group '", group_name, "' has <3 samples after filtering; correlations may be unstable.")
  }

  m <- mat[, ids, drop = FALSE]
  m <- drop_na_cols(m)

  if (ncol(m) == 0) return(matrix(numeric(0), nrow = 0, ncol = 0))

  m <- m[rowMeans(m, na.rm = TRUE) > 1, , drop = FALSE]
  m
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

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
## Load groups.xlsx (audit)
## =========================================================
groups_xlsx_df <- readxl::read_excel(groups_xlsx) %>% as.data.frame()
stopifnot("ID" %in% colnames(groups_xlsx_df))

groups_mfp <- groups_xlsx_df %>%
  transmute(
    SampleID_raw = as.character(ID),
    SampleID     = clean_sample_id(as.character(ID))
  ) %>%
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
## analyze_one_graph (GO + Hallmark + ranking + module meta-graph)
## =========================================================
analyze_one_graph <- function(g,
                              label,
                              outdir = "DCN_outputs",
                              orgDb = NULL,
                              keyType = "SYMBOL",
                              min_module_size = 20,
                              top_terms_per_module = 10,
                              top_n_nodes = 20) {

  if (!requireNamespace("igraph", quietly = TRUE)) stop("Install 'igraph'")
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Install 'openxlsx'")
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Install 'clusterProfiler' (Bioconductor)")
  if (!requireNamespace("circlize", quietly = TRUE)) stop("Install 'circlize'")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Install 'org.Hs.eg.db' (Bioconductor)")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install 'ggplot2'")
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    message("[note] 'msigdbr' not installed -> Hallmark enrichment will be skipped (install.packages('msigdbr')).")
  }

  suppressPackageStartupMessages({
    library(igraph)
    library(openxlsx)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(circlize)
    library(ggplot2)
  })

  if (is.null(orgDb)) orgDb <- org.Hs.eg.db::org.Hs.eg.db

  stopifnot(inherits(g, "igraph"))

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  gdir      <- file.path(outdir, label)
  plots_dir <- file.path(gdir, "module_plots")
  dir.create(gdir,      showWarnings = FALSE, recursive = TRUE)
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  ## ---- Clean + LCC (robust to igraph 2.x membership type) ----
  g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])

  if (igraph::ecount(g) == 0 || igraph::vcount(g) < 2) {
    warning(sprintf("[%s] No edges after cleaning.", label))
    return(invisible(NULL))
  }

  cmp <- igraph::components(g)
  memb <- cmp$membership
  if (!is.atomic(memb)) memb <- unlist(memb, use.names = FALSE)
  memb <- as.integer(memb)

  csize <- cmp$csize
  if (!is.atomic(csize)) csize <- unlist(csize, use.names = FALSE)
  csize <- as.numeric(csize)

  lcc_id <- which.max(csize)
  lcc_vids <- which(memb == lcc_id)
  g <- igraph::induced_subgraph(g, vids = lcc_vids)

  ## ---- Communities (Louvain) on |weight| ----
  igraph::E(g)$w_abs <- abs(igraph::E(g)$weight %||% 0)
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$w_abs)
  mb <- igraph::membership(cl)
  if (!is.atomic(mb)) mb <- unlist(mb, use.names = TRUE)

  igraph::V(g)$module <- as.integer(mb[igraph::V(g)$name])
  mod_sizes <- igraph::sizes(cl)
  modules <- split(names(mb), mb)

  ## ---- Node summaries ----
  str_abs <- igraph::strength(g, weights = abs(igraph::E(g)$weight %||% 0))
  deg     <- igraph::degree(g)
  w_len   <- 1 / pmax(1e-9, abs(igraph::E(g)$weight %||% 0))
  btw     <- igraph::betweenness(g, directed = FALSE, weights = w_len, normalized = FALSE)

  ord_hub <- order(str_abs, decreasing = TRUE, na.last = TRUE)[seq_len(min(top_n_nodes, igraph::vcount(g)))]
  ord_btw <- order(btw,     decreasing = TRUE, na.last = TRUE)[seq_len(min(top_n_nodes, igraph::vcount(g)))]

  df_hubs <- data.frame(
    Gene        = igraph::V(g)$name[ord_hub],
    StrengthAbs = as.numeric(str_abs[ord_hub]),
    Degree      = as.integer(deg[ord_hub]),
    stringsAsFactors = FALSE
  )
  df_btw <- data.frame(
    Gene        = igraph::V(g)$name[ord_btw],
    Betweenness = as.numeric(btw[ord_btw]),
    StrengthAbs = as.numeric(str_abs[ord_btw]),
    Degree      = as.integer(deg[ord_btw]),
    stringsAsFactors = FALSE
  )

  ## ---- Graph-level topology ----
  n_nodes <- igraph::vcount(g); n_edges <- igraph::ecount(g)
  comps   <- igraph::components(g)
  m2 <- comps$membership
  if (!is.atomic(m2)) m2 <- unlist(m2, use.names = FALSE)
  m2 <- as.integer(m2)
  cs2 <- comps$csize
  if (!is.atomic(cs2)) cs2 <- unlist(cs2, use.names = FALSE)
  cs2 <- as.numeric(cs2)

  lcc_vids2 <- which(m2 == which.max(cs2))
  g_lcc <- igraph::induced_subgraph(g, lcc_vids2)

  topo_metrics <- data.frame(
    Label = label,
    Nodes = n_nodes,
    Edges = n_edges,
    Components = comps$no,
    LCC_Nodes = igraph::vcount(g_lcc),
    LCC_Edges = igraph::ecount(g_lcc),
    Density = igraph::edge_density(g, loops = FALSE),
    Avg_Degree = mean(deg),
    Transitivity = igraph::transitivity(g, type = "global", isolates = "zero"),
    Assortativity_Degree = igraph::assortativity_degree(g, directed = FALSE),
    Vertex_Connectivity = tryCatch(igraph::vertex_connectivity(g), error = function(e) NA_real_),
    Edge_Connectivity   = tryCatch(igraph::edge_connectivity(g),   error = function(e) NA_real_),
    Avg_Path_Length_LCC = tryCatch(igraph::average.path.length(g_lcc, directed = FALSE), error = function(e) NA_real_),
    Modularity_Louvain  = tryCatch(igraph::modularity(g, membership = igraph::V(g)$module, weights = igraph::E(g)$w_abs), error = function(e) NA_real_),
    stringsAsFactors = FALSE
  )

  ## =========================================================
  ## Oncology-oriented module ranking + typing + top modules
  ## =========================================================
  module_ids <- sort(unique(igraph::V(g)$module))
  mod_rank_list <- lapply(module_ids, function(m) {
    gm <- igraph::induced_subgraph(g, vids = igraph::V(g)[module == m])
    n <- igraph::vcount(gm); e <- igraph::ecount(gm)
    if (n < 2 || e == 0) {
      return(data.frame(
        Module = as.integer(m),
        N = n, E = e,
        Density = 0,
        MeanAbsStrength = 0,
        MaxBetweenness = 0,
        HubConcentration = NA_real_,
        MeanAbsEdgeWeight = NA_real_,
        EdgeWeightSD = NA_real_,
        Clustering = NA_real_,
        CrossEdgeCount = NA_real_,
        CrossEdgeAbsSum = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    dens <- igraph::edge_density(gm, loops = FALSE)
    st <- igraph::strength(gm, weights = abs(igraph::E(gm)$weight %||% 0))
    mean_st <- mean(st, na.rm = TRUE)

    w_len_m <- 1 / pmax(1e-9, abs(igraph::E(gm)$weight %||% 0))
    btw_m <- igraph::betweenness(gm, directed = FALSE, weights = w_len_m, normalized = FALSE)
    max_btw <- max(btw_m, na.rm = TRUE)

    hub_conc <- if (sum(st, na.rm = TRUE) > 0) max(st, na.rm = TRUE) / sum(st, na.rm = TRUE) else NA_real_

    wabs <- abs(igraph::E(gm)$weight %||% 0)
    mean_w <- mean(wabs, na.rm = TRUE)
    sd_w <- stats::sd(wabs, na.rm = TRUE)

    clust <- tryCatch(igraph::transitivity(gm, type = "average"), error = function(e) NA_real_)

    # cross-module edges stats computed from full graph edge list
    ed <- as.data.frame(igraph::as_edgelist(g, names = TRUE), stringsAsFactors = FALSE)
    colnames(ed) <- c("v1","v2")
    ed$w <- igraph::E(g)$weight %||% 0
    ed$m1 <- igraph::V(g)$module[match(ed$v1, igraph::V(g)$name)]
    ed$m2 <- igraph::V(g)$module[match(ed$v2, igraph::V(g)$name)]
    cross_m <- ed[(ed$m1 == m & ed$m2 != m) | (ed$m2 == m & ed$m1 != m), , drop=FALSE]
    cross_n <- nrow(cross_m)
    cross_abs_sum <- if (cross_n) sum(abs(cross_m$w), na.rm = TRUE) else 0

    data.frame(
      Module = as.integer(m),
      N = n, E = e,
      Density = dens,
      MeanAbsStrength = mean_st,
      MaxBetweenness = max_btw,
      HubConcentration = hub_conc,
      MeanAbsEdgeWeight = mean_w,
      EdgeWeightSD = sd_w,
      Clustering = clust,
      CrossEdgeCount = cross_n,
      CrossEdgeAbsSum = cross_abs_sum,
      stringsAsFactors = FALSE
    )
  })
  module_rank_df <- do.call(rbind, mod_rank_list)
  module_rank_df$PassSize <- module_rank_df$N >= min_module_size

  z <- function(x) {
    x <- as.numeric(x)
    if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  }

  module_rank_df$Score <- with(module_rank_df,
    1.2*z(Density) +
    1.2*z(MeanAbsStrength) +
    0.8*z(MaxBetweenness) +
    0.4*z(HubConcentration) +
    0.6*z(Clustering) +
    0.6*z(CrossEdgeAbsSum)
  )

  module_rank_df <- module_rank_df %>%
    arrange(desc(PassSize), desc(Score), desc(N), desc(Density))

  # thresholds for typing
  qN75 <- stats::quantile(module_rank_df$N[module_rank_df$PassSize], 0.75, na.rm = TRUE)
  qD75 <- stats::quantile(module_rank_df$Density[module_rank_df$PassSize], 0.75, na.rm = TRUE)
  qD25 <- stats::quantile(module_rank_df$Density[module_rank_df$PassSize], 0.25, na.rm = TRUE)
  qX75 <- stats::quantile(module_rank_df$CrossEdgeCount[module_rank_df$PassSize], 0.75, na.rm = TRUE)

  module_rank_df <- module_rank_df %>%
    mutate(
      ModuleType = case_when(
        !PassSize ~ "Too small / unstable",
        HubConcentration >= 0.25 & N <= 60 ~ "Small, hub-dominated (driver-centric)",
        is.finite(qD75) & is.finite(qX75) & Density >= qD75 & CrossEdgeCount >= qX75 ~ "Dense + rewired (highest priority)",
        is.finite(qN75) & is.finite(qD75) & N >= qN75 & Density >= qD75 ~ "Large, dense (pathway-level)",
        is.finite(qN75) & is.finite(qD25) & N >= qN75 & Density <= qD25 ~ "Sparse, large (usually noise)",
        TRUE ~ "Moderate program"
      )
    )

  top_modules <- head(module_rank_df$Module[module_rank_df$PassSize], 3)
  if (length(top_modules) == 0) top_modules <- head(module_rank_df$Module, 3)

  ## =========================================================
  ## Per-module plots (layout_with_fr; robust)
  ## =========================================================
  pal <- grDevices::hcl.colors(length(unique(igraph::V(g)$module)), "Dark 3")
  names(pal) <- as.character(sort(unique(igraph::V(g)$module)))

  for (m in sort(unique(igraph::V(g)$module))) {
    gm <- igraph::induced_subgraph(g, vids = igraph::V(g)[module == m])
    if (igraph::vcount(gm) < 2 || igraph::ecount(gm) == 0) next
    if (igraph::vcount(gm) < min_module_size) next

    wabs <- abs(igraph::E(gm)$weight %||% 0)
    igraph::E(gm)$width <- 0.5
    igraph::E(gm)$color <- ifelse((igraph::E(gm)$weight %||% 0) >= 0, "forestgreen", "firebrick")
    lay <- igraph::layout_with_fr(gm, weights = wabs)

    degm <- igraph::degree(gm)
    n_max <- max(degm, na.rm = TRUE)
    vsize <- if (n_max > 0) (degm * 6) / n_max + 2 else rep(3, length(degm))
    igraph::V(gm)$color <- pal[as.character(m)]

    igraph::V(gm)$label <- NA
    top5 <- order(degm, decreasing = TRUE)[seq_len(min(5, length(degm)))]
    igraph::V(gm)$label[top5] <- igraph::V(gm)$name[top5]

    fpng <- file.path(plots_dir, sprintf("module_%s_n%d.png", m, igraph::vcount(gm)))
    grDevices::png(fpng, width = 1600, height = 1300, res = 180)
    plot(gm, layout = lay,
         vertex.size = vsize,
         vertex.label = igraph::V(gm)$label,
         vertex.label.cex = 0.8,
         vertex.label.color = "black",
         main = sprintf("%s — Module %s (n=%d) — node size=degree; labels=top5 hubs", label, m, igraph::vcount(gm)))
    grDevices::dev.off()
  }

  ## =========================================================
  ## Circlize cross-module chords (optional visualization)
  ## =========================================================
  edf <- as.data.frame(igraph::as_edgelist(g, names = TRUE), stringsAsFactors = FALSE)
  colnames(edf) <- c("v1","v2")
  edf$w <- igraph::E(g)$weight %||% 0
  edf$mod1 <- igraph::V(g)$module[match(edf$v1, igraph::V(g)$name)]
  edf$mod2 <- igraph::V(g)$module[match(edf$v2, igraph::V(g)$name)]
  cross <- edf[edf$mod1 != edf$mod2, , drop = FALSE]

  if (nrow(cross)) {
    key <- apply(cross[,c("mod1","mod2")], 1, function(z) paste(sort(z), collapse = "_"))
    w_counts <- as.data.frame(table(key), stringsAsFactors = FALSE)
    parts <- strsplit(w_counts$key, "_", fixed = TRUE)
    links_df <- data.frame(
      from = vapply(parts, `[`, "", 1),
      to   = vapply(parts, `[`, "", 2),
      value = as.numeric(w_counts$Freq),
      stringsAsFactors = FALSE
    )

    mod_sizes_vec <- setNames(as.integer(mod_sizes), names(mod_sizes))
    sector_ord <- as.character(names(sort(mod_sizes_vec, decreasing = TRUE)))
    grid_col <- setNames(grDevices::hcl.colors(length(sector_ord), "Dark 3"), sector_ord)

    outfile <- file.path(gdir, sprintf("%s_modules_circlize.png", label))
    grDevices::png(outfile, width = 1600, height = 1600, res = 200)
    circlize::circos.clear()
    circlize::circos.par(start.degree = 90, gap.degree = 4, track.margin = c(0.01, 0.01))
    circlize::chordDiagram(
      x = links_df,
      order = sector_ord,
      grid.col = grid_col,
      transparency = 0.25,
      directional = 0,
      annotationTrack = c("grid"),
      preAllocateTracks = list(track.height = 0.12)
    )
    grDevices::dev.off()
  } else {
    message(sprintf("[%s] No cross-module edges; circlize plot skipped.", label))
  }

  ## =========================================================
  ## Module meta-graph plot (modules as nodes)
  ## Node size = #genes; node color = ModuleType
  ## Edge width = sum |weights| across inter-module edges
  ## Edge color = median sign (+/-) of correlations across inter-module edges
  ## =========================================================
  build_module_metagraph <- function(g, module_rank_df, out_png) {
    ed <- as.data.frame(igraph::as_edgelist(g, names = TRUE), stringsAsFactors = FALSE)
    colnames(ed) <- c("v1","v2")
    ed$w <- igraph::E(g)$weight %||% 0
    ed$m1 <- igraph::V(g)$module[match(ed$v1, igraph::V(g)$name)]
    ed$m2 <- igraph::V(g)$module[match(ed$v2, igraph::V(g)$name)]
    ed <- ed[ed$m1 != ed$m2, , drop=FALSE]
    if (!nrow(ed)) return(invisible(NULL))

    key <- apply(ed[,c("m1","m2")], 1, function(z) paste(sort(z), collapse = "__"))
    ed$key <- key

    agg <- ed %>%
      group_by(key) %>%
      summarize(
        mA = min(m1, m2),
        mB = max(m1, m2),
        EdgeAbsSum = sum(abs(w), na.rm = TRUE),
        MedianW = stats::median(w, na.rm = TRUE),
        .groups = "drop"
      )

    verts <- module_rank_df %>%
      distinct(Module, N, ModuleType) %>%
      mutate(name = as.character(Module))

    gM <- igraph::graph_from_data_frame(
      d = data.frame(from = as.character(agg$mA),
                     to   = as.character(agg$mB),
                     weight = agg$EdgeAbsSum,
                     medw   = agg$MedianW,
                     stringsAsFactors = FALSE),
      directed = FALSE,
      vertices = verts
    )

    # aesthetics
    # node sizes
    N <- igraph::V(gM)$N %||% 1
    N <- as.numeric(N); N[is.na(N)] <- 1
    vsz <- 6 + 18 * (N / max(N))
    igraph::V(gM)$size <- vsz

    # node colors by ModuleType
    type_levels <- unique(igraph::V(gM)$ModuleType)
    type_cols <- setNames(grDevices::hcl.colors(max(3, length(type_levels)), "Set 2")[seq_along(type_levels)],
                          type_levels)
    igraph::V(gM)$color <- type_cols[igraph::V(gM)$ModuleType]

    # edge width by abs sum
    w <- igraph::E(gM)$weight %||% 1
    w <- as.numeric(w); w[is.na(w)] <- 0
    ew <- 0.5 + 6 * (w / max(w))
    igraph::E(gM)$width <- ew

    # edge color by median sign
    medw <- igraph::E(gM)$medw %||% 0
    igraph::E(gM)$color <- ifelse(as.numeric(medw) >= 0, "forestgreen", "firebrick")

    lay <- igraph::layout_with_fr(gM, weights = igraph::E(gM)$weight %||% 1)

    grDevices::png(out_png, width = 1800, height = 1400, res = 200)
    plot(gM, layout = lay,
         vertex.label = igraph::V(gM)$name,
         vertex.label.cex = 0.8,
         vertex.label.color = "black",
         main = sprintf("%s — Module meta-graph (node=size genes; edge=sum|w|; color=median sign)", label))
    legend("topleft", legend = names(type_cols), col = type_cols, pch = 19, bty = "n", cex = 0.85)
    grDevices::dev.off()

    invisible(gM)
  }

  meta_png <- file.path(gdir, sprintf("%s_module_metagraph.png", label))
  g_module <- build_module_metagraph(g, module_rank_df, meta_png)

  ## =========================================================
  ## Enrichment: GO:BP per module
  ## =========================================================
  go_enrichments <- lapply(names(modules), function(m) {
    genes <- modules[[m]]
    if (length(genes) < min_module_size) return(NULL)
    suppressWarnings(tryCatch({
      clusterProfiler::enrichGO(
        gene          = genes,
        OrgDb         = orgDb,
        keyType       = keyType,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      )
    }, error = function(e) NULL))
  })
  names(go_enrichments) <- paste0("Module_", names(modules))

  ## =========================================================
  ## Enrichment: Hallmark (MSigDB H) per module (optional)
  ## =========================================================
  hallmark_enrichments <- NULL
  hallmark_summary_df <- NULL
  hallmark_term2gene <- NULL

  if (requireNamespace("msigdbr", quietly = TRUE)) {
    msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    hallmark_term2gene <- msig %>% dplyr::select(gs_name, gene_symbol) %>% distinct()
    hallmark_enrichments <- lapply(names(modules), function(m) {
      genes <- modules[[m]]
      if (length(genes) < min_module_size) return(NULL)
      suppressWarnings(tryCatch({
        clusterProfiler::enricher(
          gene = genes,
          TERM2GENE = hallmark_term2gene,
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.05,
          qvalueCutoff  = 0.2
        )
      }, error = function(e) NULL))
    })
    names(hallmark_enrichments) <- paste0("Module_", names(modules))
  }

  ## ---- Summarizer ----
  make_enrichment_summary <- function(enrichments, top_n = 10) {
    out <- lapply(names(enrichments), function(mname) {
      enr <- enrichments[[mname]]
      if (is.null(enr)) return(NULL)
      df <- as.data.frame(enr); if (!nrow(df)) return(NULL)
      mod_num <- suppressWarnings(as.integer(sub("^.*?(\\d+).*?$", "\\1", mname)))
      module  <- ifelse(is.na(mod_num), mname, mod_num)
      keep_cols <- intersect(c("ID","Description","p.adjust","pvalue","qvalue",
                               "GeneRatio","BgRatio","Count","geneID"), names(df))
      df <- df[, keep_cols, drop = FALSE]
      if ("p.adjust" %in% names(df)) df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
      df <- head(df, top_n)
      data.frame(
        Module    = module,
        Term      = df$Description %||% df$ID %||% NA,
        Adj_P     = df$p.adjust   %||% NA,
        ID        = df$ID         %||% NA,
        P_value   = df$pvalue     %||% NA,
        Q_value   = df$qvalue     %||% NA,
        GeneRatio = df$GeneRatio  %||% NA,
        BgRatio   = df$BgRatio    %||% NA,
        Count     = df$Count      %||% NA,
        Genes     = df$geneID     %||% NA,
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, out)
    if (is.null(out)) {
      out <- data.frame(Module=integer(), Term=character(), Adj_P=numeric(),
                        ID=character(), P_value=numeric(), Q_value=numeric(),
                        GeneRatio=character(), BgRatio=character(), Count=integer(),
                        Genes=character(), stringsAsFactors = FALSE)
    } else {
      out <- out[order(out$Module, out$Adj_P, na.last = TRUE), ]
    }
    rownames(out) <- NULL
    out
  }

  go_summary_df <- make_enrichment_summary(go_enrichments, top_n = top_terms_per_module)
  if (!is.null(hallmark_enrichments)) {
    hallmark_summary_df <- make_enrichment_summary(hallmark_enrichments, top_n = top_terms_per_module)
  }

  ## =========================================================
  ## Excel output
  ## =========================================================
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "Modules")
  openxlsx::writeData(wb, "Modules",
                      data.frame(Gene = names(mb), Module = as.integer(mb), stringsAsFactors = FALSE))

  openxlsx::addWorksheet(wb, "Module_Ranking")
  openxlsx::writeData(wb, "Module_Ranking", module_rank_df)
  try({
    openxlsx::addStyle(wb, "Module_Ranking", openxlsx::createStyle(textDecoration="bold"),
                       rows = 1, cols = 1:ncol(module_rank_df), gridExpand = TRUE)
    openxlsx::freezePane(wb, "Module_Ranking", firstRow = TRUE)
  }, silent = TRUE)

  openxlsx::addWorksheet(wb, "Top20_Hubs_Betweenness")
  openxlsx::writeData(wb, "Top20_Hubs_Betweenness", "Top hubs (by |strength|)", startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "Top20_Hubs_Betweenness", df_hubs, startRow = 2, startCol = 1)
  start_row_btw <- nrow(df_hubs) + 4
  openxlsx::writeData(wb, "Top20_Hubs_Betweenness", "Top betweenness", startRow = start_row_btw, startCol = 1)
  openxlsx::writeData(wb, "Top20_Hubs_Betweenness", df_btw, startRow = start_row_btw + 1, startCol = 1)

  openxlsx::addWorksheet(wb, "Graph_Topology")
  openxlsx::writeData(wb, "Graph_Topology", topo_metrics)

  openxlsx::addWorksheet(wb, "GO_BP_Summary")
  openxlsx::writeData(wb, "GO_BP_Summary", go_summary_df)
  openxlsx::addStyle(wb, "GO_BP_Summary", openxlsx::createStyle(textDecoration = "bold"),
                     rows = 1, cols = 1:ncol(go_summary_df), gridExpand = TRUE)
  openxlsx::freezePane(wb, "GO_BP_Summary", firstRow = TRUE)

  # write per-module GO tables (optional; can be large)
  for (m in names(go_enrichments)) {
    enr <- go_enrichments[[m]]
    if (!is.null(enr)) {
      sh <- paste0("GO_", m)
      sh <- substr(sh, 1, 31)
      openxlsx::addWorksheet(wb, sh)
      openxlsx::writeData(wb, sh, as.data.frame(enr))
    }
  }

  if (!is.null(hallmark_summary_df)) {
    openxlsx::addWorksheet(wb, "HALLMARK_Summary")
    openxlsx::writeData(wb, "HALLMARK_Summary", hallmark_summary_df)
    openxlsx::addStyle(wb, "HALLMARK_Summary", openxlsx::createStyle(textDecoration = "bold"),
                       rows = 1, cols = 1:ncol(hallmark_summary_df), gridExpand = TRUE)
    openxlsx::freezePane(wb, "HALLMARK_Summary", firstRow = TRUE)

    # per-module hallmark tables
    for (m in names(hallmark_enrichments)) {
      enr <- hallmark_enrichments[[m]]
      if (!is.null(enr)) {
        sh <- paste0("H_", m)
        sh <- substr(sh, 1, 31)
        openxlsx::addWorksheet(wb, sh)
        openxlsx::writeData(wb, sh, as.data.frame(enr))
      }
    }
  }

  # type legend sheet (your requested semantics)
  openxlsx::addWorksheet(wb, "ModuleType_Legend")
  legend_df <- data.frame(
    ModuleType = c("Small, hub-dominated (driver-centric)",
                   "Large, dense (pathway-level)",
                   "Sparse, large (usually noise)",
                   "Dense + rewired (highest priority)"),
    BiologicalRelevance = c("Driver-centric", "Pathway-level", "Low", "Very high"),
    OncologyInterest = c("High if tumor-specific", "High if rewired", "Usually noise", "Highest priority"),
    stringsAsFactors = FALSE
  )
  openxlsx::writeData(wb, "ModuleType_Legend", legend_df)

  xlsx_path <- file.path(gdir, sprintf("%s_results.xlsx", label))
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  ## save objects
  saveRDS(g,  file = file.path(gdir, sprintf("%s_graph_cleaned.rds", label)))
  saveRDS(cl, file = file.path(gdir, sprintf("%s_louvain_object.rds", label)))
  saveRDS(go_enrichments, file = file.path(gdir, sprintf("%s_GO_enrichments.rds", label)))
  if (!is.null(hallmark_enrichments)) saveRDS(hallmark_enrichments, file = file.path(gdir, sprintf("%s_HALLMARK_enrichments.rds", label)))
  if (!is.null(g_module)) saveRDS(g_module, file = file.path(gdir, sprintf("%s_module_metagraph.rds", label)))

  message(sprintf("[%s] Done. Excel: %s | Module plots: %s | Module meta-graph: %s",
                  label, xlsx_path, plots_dir, meta_png))

  invisible(list(
    graph=g, membership=mb, sizes=mod_sizes,
    hubs=df_hubs, betweenness=df_btw, topo=topo_metrics,
    module_ranking=module_rank_df,
    go_enrichments=go_enrichments, go_summary=go_summary_df,
    hallmark_enrichments=hallmark_enrichments, hallmark_summary=hallmark_summary_df,
    excel=xlsx_path, plots_dir=plots_dir, module_metagraph_png=meta_png
  ))
}

## =========================================================
## compare_networks_plots (robust)
## =========================================================
compare_networks_plots <- function(g1, g2, label1="Net1", label2="Net2", outdir="NET_COMPARE") {
  stopifnot(inherits(g1, "igraph"), inherits(g2, "igraph"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2")
  if (!requireNamespace("reshape2", quietly = TRUE)) stop("Install reshape2")
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap")

  suppressPackageStartupMessages({
    library(igraph)
    library(ggplot2)
    library(reshape2)
    library(pheatmap)
  })

  safe_lcc <- function(g) {
    g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
    if (igraph::ecount(g) == 0 || igraph::vcount(g) < 2) return(g)
    c <- igraph::components(g)
    memb <- c$membership
    if (!is.atomic(memb)) memb <- unlist(memb, use.names = FALSE)
    memb <- as.integer(memb)
    cs <- c$csize
    if (!is.atomic(cs)) cs <- unlist(cs, use.names = FALSE)
    cs <- as.numeric(cs)
    igraph::induced_subgraph(g, which(memb == which.max(cs)))
  }

  global_metrics <- function(g) {
    Nodes <- igraph::vcount(g); Edges <- igraph::ecount(g)
    if (Nodes == 0 || Edges == 0) {
      return(data.frame(
        Nodes=Nodes, Edges=Edges,
        Density=NA_real_, Transitivity=NA_real_, Avg_Degree=NA_real_, Avg_Path_LCC=NA_real_,
        Assortativity=NA_real_, VertexConn=NA_real_, EdgeConn=NA_real_, Modularity=NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    deg <- igraph::degree(g)
    g_lcc <- safe_lcc(g)
    apl <- tryCatch(igraph::average.path.length(g_lcc, directed = FALSE), error = function(e) NA_real_)
    v_conn <- tryCatch(igraph::vertex_connectivity(g), error = function(e) NA_real_)
    e_conn <- tryCatch(igraph::edge_connectivity(g), error = function(e) NA_real_)
    mod <- tryCatch({
      cl <- igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight %||% 1))
      igraph::modularity(g, igraph::membership(cl), weights = abs(igraph::E(g)$weight %||% 1))
    }, error = function(e) NA_real_)

    data.frame(
      Nodes=Nodes, Edges=Edges,
      Density=igraph::edge_density(g, loops = FALSE),
      Transitivity=igraph::transitivity(g, type="global", isolates="zero"),
      Avg_Degree=mean(deg),
      Avg_Path_LCC=apl,
      Assortativity=igraph::assortativity_degree(g, directed=FALSE),
      VertexConn=v_conn,
      EdgeConn=e_conn,
      Modularity=mod,
      stringsAsFactors = FALSE
    )
  }

  node_metrics <- function(g) {
    if (igraph::vcount(g) == 0) return(data.frame(gene=character(), degree=integer(), strength=numeric(), btw=numeric()))
    data.frame(
      gene = igraph::V(g)$name %||% as.character(seq_len(igraph::vcount(g))),
      degree = igraph::degree(g),
      strength = igraph::strength(g, weights = abs(igraph::E(g)$weight %||% 1)),
      btw = igraph::betweenness(g, directed=FALSE, weights = 1/pmax(1e-9, abs(igraph::E(g)$weight %||% 1)), normalized=FALSE),
      stringsAsFactors = FALSE
    )
  }

  gm1 <- global_metrics(g1); gm2 <- global_metrics(g2)
  gm1$Network <- label1; gm2$Network <- label2
  gmtab <- rbind(gm1, gm2)

  nm1 <- node_metrics(g1); nm1$Network <- label1
  nm2 <- node_metrics(g2); nm2$Network <- label2
  nodes_merged <- merge(nm1, nm2, by="gene", all=FALSE, suffixes=c(".1",".2"))

  gml <- reshape2::melt(gmtab, id.vars="Network")
  gml <- subset(gml, variable %in% c("Nodes","Edges","Density","Transitivity",
                                     "Avg_Degree","Avg_Path_LCC","Assortativity","Modularity"))
  p_global <- ggplot(gml, aes(x=variable, y=value, fill=Network)) +
    geom_col(position=position_dodge(width=0.7), width=0.6) +
    coord_flip() +
    labs(x=NULL, y=NULL, title="Global topology comparison") +
    theme_minimal(base_size=12)
  ggsave(file.path(outdir, "01_global_topology.png"), p_global, width=10, height=6, dpi=200)

  ecdf_long <- rbind(
    data.frame(Network=label1, degree=nm1$degree, strength=nm1$strength, btw=nm1$btw),
    data.frame(Network=label2, degree=nm2$degree, strength=nm2$strength, btw=nm2$btw)
  )
  plot_ecdf <- function(col, ttl) {
    ggplot(ecdf_long, aes(x=log1p(.data[[col]]), color=Network)) +
      stat_ecdf(geom="step", linewidth=1) +
      labs(x=paste0("log1p(",ttl,")"), y="ECDF", title=paste0("ECDF of ", ttl)) +
      theme_minimal(base_size=12)
  }
  ggsave(file.path(outdir,"02_ecdf_degree.png"),      plot_ecdf("degree","degree"),       width=6, height=4, dpi=200)
  ggsave(file.path(outdir,"03_ecdf_strength.png"),    plot_ecdf("strength","strength"),   width=6, height=4, dpi=200)
  ggsave(file.path(outdir,"04_ecdf_betweenness.png"), plot_ecdf("btw","betweenness"),     width=6, height=4, dpi=200)

  if (nrow(nodes_merged) > 0) {
    scat <- function(xcol, ycol, ttl) {
      ggplot(nodes_merged, aes(x=log1p(.data[[xcol]]), y=log1p(.data[[ycol]]))) +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="gray60") +
        geom_point(alpha=0.5, size=1.2) +
        labs(x=paste0(label1," log1p(",ttl,")"), y=paste0(label2," log1p(",ttl,")"),
             title=paste0(ttl, ": ", label1, " vs ", label2)) +
        theme_minimal(base_size=12)
    }
    ggsave(file.path(outdir,"05_scatter_degree.png"),   scat("degree.1","degree.2","Degree"),   width=6, height=5, dpi=200)
    ggsave(file.path(outdir,"06_scatter_strength.png"), scat("strength.1","strength.2","Strength"), width=6, height=5, dpi=200)
  }

  # Module correspondence (Jaccard) — only if both graphs non-empty
  jaccard_modules <- function(gA, gB) {
    clA <- igraph::cluster_louvain(gA, weights = abs(igraph::E(gA)$weight %||% 1))
    clB <- igraph::cluster_louvain(gB, weights = abs(igraph::E(gB)$weight %||% 1))
    setsA <- split(igraph::V(gA)$name, igraph::membership(clA))
    setsB <- split(igraph::V(gB)$name, igraph::membership(clB))
    A <- names(setsA); B <- names(setsB)
    J <- matrix(0, nrow=length(A), ncol=length(B), dimnames=list(paste0("A",A), paste0("B",B)))
    for (i in seq_along(setsA)) for (j in seq_along(setsB)) {
      inter <- length(intersect(setsA[[i]], setsB[[j]]))
      uni   <- length(union(setsA[[i]], setsB[[j]]))
      J[i,j] <- if (uni>0) inter/uni else 0
    }
    J
  }
  if (igraph::vcount(g1)>0 && igraph::vcount(g2)>0 && igraph::ecount(g1)>0 && igraph::ecount(g2)>0) {
    J <- jaccard_modules(g1, g2)
    pheatmap::pheatmap(J, cluster_rows = TRUE, cluster_cols = TRUE,
                       main = paste("Module overlap (Jaccard):", label1, "vs", label2),
                       filename = file.path(outdir, "07_module_overlap_jaccard.png"),
                       width = 8, height = 6)
  }

  edge_key <- function(g) {
    if (igraph::ecount(g)==0) return(character())
    el <- igraph::as_edgelist(g, names=TRUE)
    apply(el, 1, function(z) paste(sort(z), collapse="__"))
  }
  e1 <- edge_key(g1); e2 <- edge_key(g2)
  edge_df <- data.frame(
    Class = c(paste0("Overlap (",label1,"∩",label2,")"), paste0("Only ",label1), paste0("Only ",label2)),
    Count = c(length(intersect(e1,e2)), length(setdiff(e1,e2)), length(setdiff(e2,e1))),
    stringsAsFactors = FALSE
  )
  p_edges <- ggplot(edge_df, aes(x=Class, y=Count, fill=Class)) +
    geom_col(width=0.6) + theme_minimal(base_size=12) +
    guides(fill="none") + labs(x=NULL, y="Edges", title="Edge overlap / rewiring")
  ggsave(file.path(outdir, "08_edge_overlap.png"), p_edges, width=7, height=4, dpi=200)

  invisible(list(global=gmtab, nodes_merged=nodes_merged))
}

## =========================================================
## Run DCN (cached)
## =========================================================
dcn_params <- list(
  alpha = 0.01,
  rho = NULL,
  t_percentile = 0.90,
  min_overlap = 3,
  require_both_above_t = FALSE,
  diff_alpha = 0.01,
  diff_t = 0.3,
  fdr_method = "BH",
  diff_fdr_method = "BH"
)

res_EC_vs_Healthy <- run_or_load_diff_graph(
  label1 = "EC", label2 = "Healthy",
  M1 = gene_expression_matrix_EC,
  M2 = gene_expression_matrix_healthy,
  out_prefix = out_prefix,
  params = dcn_params,
  overwrite = overwrite
)

res_Dec_vs_Healthy <- run_or_load_diff_graph(
  label1 = "Decidua", label2 = "Healthy",
  M1 = gene_expression_matrix_decidua,
  M2 = gene_expression_matrix_healthy,
  out_prefix = out_prefix,
  params = dcn_params,
  overwrite = overwrite
)

res_EC_vs_Dec <- run_or_load_diff_graph(
  label1 = "EC", label2 = "Decidua",
  M1 = gene_expression_matrix_EC,
  M2 = gene_expression_matrix_decidua,
  out_prefix = out_prefix,
  params = dcn_params,
  overwrite = overwrite
)

## legacy filenames (optional)
saveRDS(res_EC_vs_Healthy$graph_group1,     file = "graph_groupEC_rho06.rds")
saveRDS(res_EC_vs_Healthy$graph_group2,     file = "graph_groupHealthy_rho06.rds")
saveRDS(res_EC_vs_Healthy$graph_difference, file = "graph_difference_ECvsHealthy.rds")

saveRDS(res_EC_vs_Dec$graph_group1,         file = "graph_groupEC_rho06_ECvsDecidua.rds")
saveRDS(res_EC_vs_Dec$graph_group2,         file = "graph_groupDecidua_rho06.rds")
saveRDS(res_EC_vs_Dec$graph_difference,     file = "graph_difference_ECvsDecidua.rds")

## =========================================================
## Analyze graphs + compare
## =========================================================
gEC  <- res_EC_vs_Healthy$graph_group1
gH   <- res_EC_vs_Healthy$graph_group2
gDec <- res_EC_vs_Dec$graph_group2

res_EC  <- analyze_one_graph(gEC,  label="EC",      outdir="DCN_outputs",
                            min_module_size=20, top_terms_per_module=10, top_n_nodes=20)
res_H   <- analyze_one_graph(gH,   label="Healthy", outdir="DCN_outputs",
                            min_module_size=20, top_terms_per_module=10, top_n_nodes=20)
res_Dec <- analyze_one_graph(gDec, label="Decidua", outdir="DCN_outputs",
                            min_module_size=20, top_terms_per_module=10, top_n_nodes=20)

compare_networks_plots(gEC, gDec, label1="EC", label2="Decidua", outdir="EC_vs_Decidua_compare")
compare_networks_plots(gEC, gH,   label1="EC", label2="Healthy", outdir="EC_vs_Healthy_compare")
compare_networks_plots(gDec,gH,   label1="Decidua", label2="Healthy", outdir="Decidua_vs_Healthy_compare")

message("DONE.")

## =========================================================
## Module meta-graph plot (modules as nodes)
## Node size = #genes; node color = ModuleType
## Edge width = sum |weights| across inter-module edges
## Edge color = median sign (+/-) of correlations across inter-module edges
## =========================================================
build_module_metagraph <- function(g, module_rank_df, out_png) {
    ed <- as.data.frame(igraph::as_edgelist(g, names = TRUE), stringsAsFactors = FALSE)
    colnames(ed) <- c("v1","v2")
    ed$w <- igraph::E(g)$weight %||% 0
    ed$m1 <- igraph::V(g)$module[match(ed$v1, igraph::V(g)$name)]
    ed$m2 <- igraph::V(g)$module[match(ed$v2, igraph::V(g)$name)]
    ed <- ed[ed$m1 != ed$m2, , drop=FALSE]
    if (!nrow(ed)) return(invisible(NULL))

    key <- apply(ed[,c("m1","m2")], 1, function(z) paste(sort(z), collapse = "__"))
    ed$key <- key

    agg <- ed %>%
      group_by(key) %>%
      summarize(
        mA = min(m1, m2),
        mB = max(m1, m2),
        EdgeAbsSum = sum(abs(w), na.rm = TRUE),
        MedianW = stats::median(w, na.rm = TRUE),
        .groups = "drop"
      )

    verts <- module_rank_df %>%
      distinct(Module, N, ModuleType) %>%
      mutate(name = as.character(Module))

    gM <- igraph::graph_from_data_frame(
      d = data.frame(from = as.character(agg$mA),
                     to   = as.character(agg$mB),
                     weight = agg$EdgeAbsSum,
                     medw   = agg$MedianW,
                     stringsAsFactors = FALSE),
      directed = FALSE,
      vertices = verts
    )

    # aesthetics
    # node sizes
    N <- igraph::V(gM)$N %||% 1
    N <- as.numeric(N); N[is.na(N)] <- 1
    vsz <- 6 + 18 * (N / max(N))
    igraph::V(gM)$size <- vsz

    # node colors by ModuleType
    type_levels <- unique(igraph::V(gM)$ModuleType)
    type_cols <- setNames(grDevices::hcl.colors(max(3, length(type_levels)), "Set 2")[seq_along(type_levels)],
                          type_levels)
    igraph::V(gM)$color <- type_cols[igraph::V(gM)$ModuleType]

    # edge width by abs sum
    w <- igraph::E(gM)$weight %||% 1
    w <- as.numeric(w); w[is.na(w)] <- 0
    ew <- 0.5 + 6 * (w / max(w))
    igraph::E(gM)$width <- ew

    # edge color by median sign
    medw <- igraph::E(gM)$medw %||% 0
    igraph::E(gM)$color <- ifelse(as.numeric(medw) >= 0, "forestgreen", "firebrick")

    lay <- igraph::layout_with_fr(gM, weights = igraph::E(gM)$weight %||% 1)

    grDevices::png(out_png, width = 1800, height = 1400, res = 200)
    plot(gM, layout = lay,
         vertex.label = igraph::V(gM)$name,
         vertex.label.cex = 0.8,
         vertex.label.color = "black",
         main = sprintf("%s — Module meta-graph (node=size genes; edge=sum|w|; color=median sign)", label))
    legend("topleft", legend = names(type_cols), col = type_cols, pch = 19, bty = "n", cex = 0.85)
    grDevices::dev.off()

    invisible(gM)
  }

  meta_png <- file.path(gdir, sprintf("%s_module_metagraph.png", label))
  g_module <- build_module_metagraph(g, module_rank_df, meta_png)

  ## =========================================================
  ## Enrichment: GO:BP per module
  ## =========================================================
  go_enrichments <- lapply(names(modules), function(m) {
    genes <- modules[[m]]
    if (length(genes) < min_module_size) return(NULL)
    suppressWarnings(tryCatch({
      clusterProfiler::enrichGO(
        gene          = genes,
        OrgDb         = orgDb,
        keyType       = keyType,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      )
    }, error = function(e) NULL))
  })
  names(go_enrichments) <- paste0("Module_", names(modules))

  ## =========================================================
  ## Enrichment: Hallmark (MSigDB H) per module (optional)
  ## =========================================================
  hallmark_enrichments <- NULL
  hallmark_summary_df <- NULL
  hallmark_term2gene <- NULL

  if (requireNamespace("msigdbr", quietly = TRUE)) {
    msig <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    hallmark_term2gene <- msig %>% dplyr::select(gs_name, gene_symbol) %>% distinct()
    hallmark_enrichments <- lapply(names(modules), function(m) {
      genes <- modules[[m]]
      if (length(genes) < min_module_size) return(NULL)
      suppressWarnings(tryCatch({
        clusterProfiler::enricher(
          gene = genes,
          TERM2GENE = hallmark_term2gene,
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.05,
          qvalueCutoff  = 0.2
        )
      }, error = function(e) NULL))
    })
    names(hallmark_enrichments) <- paste0("Module_", names(modules))
  }

  ## ---- Summarizer ----
  make_enrichment_summary <- function(enrichments, top_n = 10) {
    out <- lapply(names(enrichments), function(mname) {
      enr <- enrichments[[mname]]
      if (is.null(enr)) return(NULL)
      df <- as.data.frame(enr); if (!nrow(df)) return(NULL)
      mod_num <- suppressWarnings(as.integer(sub("^.*?(\\d+).*?$", "\\1", mname)))
      module  <- ifelse(is.na(mod_num), mname, mod_num)
      keep_cols <- intersect(c("ID","Description","p.adjust","pvalue","qvalue",
                               "GeneRatio","BgRatio","Count","geneID"), names(df))
      df <- df[, keep_cols, drop = FALSE]
      if ("p.adjust" %in% names(df)) df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
      df <- head(df, top_n)
      data.frame(
        Module    = module,
        Term      = df$Description %||% df$ID %||% NA,
        Adj_P     = df$p.adjust   %||% NA,
        ID        = df$ID         %||% NA,
        P_value   = df$pvalue     %||% NA,
        Q_value   = df$qvalue     %||% NA,
        GeneRatio = df$GeneRatio  %||% NA,
        BgRatio   = df$BgRatio    %||% NA,
        Count     = df$Count      %||% NA,
        Genes     = df$geneID     %||% NA,
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, out)
    if (is.null(out)) {
      out <- data.frame(Module=integer(), Term=character(), Adj_P=numeric(),
                        ID=character(), P_value=numeric(), Q_value=numeric(),
                        GeneRatio=character(), BgRatio=character(), Count=integer(),
                        Genes=character(), stringsAsFactors = FALSE)
    } else {
      out <- out[order(out$Module, out$Adj_P, na.last = TRUE), ]
    }
    rownames(out) <- NULL
    out
  }

  go_summary_df <- make_enrichment_summary(go_enrichments, top_n = top_terms_per_module)
  if (!is.null(hallmark_enrichments)) {
    hallmark_summary_df <- make_enrichment_summary(hallmark_enrichments, top_n = top_terms_per_module)
  }

  ## =========================================================
  ## Excel output
  ## =========================================================
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "Modules")
  openxlsx::writeData(wb, "Modules",
                      data.frame(Gene = names(mb), Module = as.integer(mb), stringsAsFactors = FALSE))

  openxlsx::addWorksheet(wb, "Module_Ranking")
  openxlsx::writeData(wb, "Module_Ranking", module_rank_df)
  try({
    openxlsx::addStyle(wb, "Module_Ranking", openxlsx::createStyle(textDecoration="bold"),
                       rows = 1, cols = 1:ncol(module_rank_df), gridExpand = TRUE)
    openxlsx::freezePane(wb, "Module_Ranking", firstRow = TRUE)
  }, silent = TRUE)

  openxlsx::addWorksheet(wb, "Top20_Hubs_Betweenness")
  openxlsx::writeData(wb, "Top20_Hubs_Betweenness", "Top hubs (by |strength|)", startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "Top20_Hubs_Betweenness", df_hubs, startRow = 2, startCol = 1)
  start_row_btw <- nrow(df_hubs) + 4
  openxlsx::writeData(wb, "Top20_Hubs_Betweenness", "Top betweenness", startRow = start_row_btw, startCol = 1)
  openxlsx::writeData(wb, "Top20_Hubs_Betweenness", df_btw, startRow = start_row_btw + 1, startCol = 1)

  openxlsx::addWorksheet(wb, "Graph_Topology")
  openxlsx::writeData(wb, "Graph_Topology", topo_metrics)

  openxlsx::addWorksheet(wb, "GO_BP_Summary")
  openxlsx::writeData(wb, "GO_BP_Summary", go_summary_df)
  openxlsx::addStyle(wb, "GO_BP_Summary", openxlsx::createStyle(textDecoration = "bold"),
                     rows = 1, cols = 1:ncol(go_summary_df), gridExpand = TRUE)
  openxlsx::freezePane(wb, "GO_BP_Summary", firstRow = TRUE)

  # write per-module GO tables (optional; can be large)
  for (m in names(go_enrichments)) {
    enr <- go_enrichments[[m]]
    if (!is.null(enr)) {
      sh <- paste0("GO_", m)
      sh <- substr(sh, 1, 31)
      openxlsx::addWorksheet(wb, sh)
      openxlsx::writeData(wb, sh, as.data.frame(enr))
    }
  }

  if (!is.null(hallmark_summary_df)) {
    openxlsx::addWorksheet(wb, "HALLMARK_Summary")
    openxlsx::writeData(wb, "HALLMARK_Summary", hallmark_summary_df)
    openxlsx::addStyle(wb, "HALLMARK_Summary", openxlsx::createStyle(textDecoration = "bold"),
                       rows = 1, cols = 1:ncol(hallmark_summary_df), gridExpand = TRUE)
    openxlsx::freezePane(wb, "HALLMARK_Summary", firstRow = TRUE)

    # per-module hallmark tables
    for (m in names(hallmark_enrichments)) {
      enr <- hallmark_enrichments[[m]]
      if (!is.null(enr)) {
        sh <- paste0("H_", m)
        sh <- substr(sh, 1, 31)
        openxlsx::addWorksheet(wb, sh)
        openxlsx::writeData(wb, sh, as.data.frame(enr))
      }
    }
  }

  # type legend sheet (your requested semantics)
  openxlsx::addWorksheet(wb, "ModuleType_Legend")
  legend_df <- data.frame(
    ModuleType = c("Small, hub-dominated (driver-centric)",
                   "Large, dense (pathway-level)",
                   "Sparse, large (usually noise)",
                   "Dense + rewired (highest priority)"),
    BiologicalRelevance = c("Driver-centric", "Pathway-level", "Low", "Very high"),
    OncologyInterest = c("High if tumor-specific", "High if rewired", "Usually noise", "Highest priority"),
    stringsAsFactors = FALSE
  )
  openxlsx::writeData(wb, "ModuleType_Legend", legend_df)

  xlsx_path <- file.path(gdir, sprintf("%s_results.xlsx", label))
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  ## save objects
  saveRDS(g,  file = file.path(gdir, sprintf("%s_graph_cleaned.rds", label)))
  saveRDS(cl, file = file.path(gdir, sprintf("%s_louvain_object.rds", label)))
  saveRDS(go_enrichments, file = file.path(gdir, sprintf("%s_GO_enrichments.rds", label)))
  if (!is.null(hallmark_enrichments)) saveRDS(hallmark_enrichments, file = file.path(gdir, sprintf("%s_HALLMARK_enrichments.rds", label)))
  if (!is.null(g_module)) saveRDS(g_module, file = file.path(gdir, sprintf("%s_module_metagraph.rds", label)))

  message(sprintf("[%s] Done. Excel: %s | Module plots: %s | Module meta-graph: %s",
                  label, xlsx_path, plots_dir, meta_png))

  invisible(list(
    graph=g, membership=mb, sizes=mod_sizes,
    hubs=df_hubs, betweenness=df_btw, topo=topo_metrics,
    module_ranking=module_rank_df,
    go_enrichments=go_enrichments, go_summary=go_summary_df,
    hallmark_enrichments=hallmark_enrichments, hallmark_summary=hallmark_summary_df,
    excel=xlsx_path, plots_dir=plots_dir, module_metagraph_png=meta_png
  ))
}
## =========================================================
## compare_networks_plots (robust)
## =========================================================
compare_networks_plots <- function(g1, g2, label1="Net1", label2="Net2", outdir="NET_COMPARE") {
  stopifnot(inherits(g1, "igraph"), inherits(g2, "igraph"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2")
  if (!requireNamespace("reshape2", quietly = TRUE)) stop("Install reshape2")
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap")

  suppressPackageStartupMessages({
    library(igraph)
    library(ggplot2)
    library(reshape2)
    library(pheatmap)
  })

  safe_lcc <- function(g) {
    g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
    if (igraph::ecount(g) == 0 || igraph::vcount(g) < 2) return(g)
    c <- igraph::components(g)
    memb <- c$membership
    if (!is.atomic(memb)) memb <- unlist(memb, use.names = FALSE)
    memb <- as.integer(memb)
    cs <- c$csize
    if (!is.atomic(cs)) cs <- unlist(cs, use.names = FALSE)
    cs <- as.numeric(cs)
    igraph::induced_subgraph(g, which(memb == which.max(cs)))
  }

  global_metrics <- function(g) {
    Nodes <- igraph::vcount(g); Edges <- igraph::ecount(g)
    if (Nodes == 0 || Edges == 0) {
      return(data.frame(
        Nodes=Nodes, Edges=Edges,
        Density=NA_real_, Transitivity=NA_real_, Avg_Degree=NA_real_, Avg_Path_LCC=NA_real_,
        Assortativity=NA_real_, VertexConn=NA_real_, EdgeConn=NA_real_, Modularity=NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    deg <- igraph::degree(g)
    g_lcc <- safe_lcc(g)
    apl <- tryCatch(igraph::average.path.length(g_lcc, directed = FALSE), error = function(e) NA_real_)
    v_conn <- tryCatch(igraph::vertex_connectivity(g), error = function(e) NA_real_)
    e_conn <- tryCatch(igraph::edge_connectivity(g), error = function(e) NA_real_)
    mod <- tryCatch({
      cl <- igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight %||% 1))
      igraph::modularity(g, igraph::membership(cl), weights = abs(igraph::E(g)$weight %||% 1))
    }, error = function(e) NA_real_)

    data.frame(
      Nodes=Nodes, Edges=Edges,
      Density=igraph::edge_density(g, loops = FALSE),
      Transitivity=igraph::transitivity(g, type="global", isolates="zero"),
      Avg_Degree=mean(deg),
      Avg_Path_LCC=apl,
      Assortativity=igraph::assortativity_degree(g, directed=FALSE),
      VertexConn=v_conn,
      EdgeConn=e_conn,
      Modularity=mod,
      stringsAsFactors = FALSE
    )
  }

  node_metrics <- function(g) {
    if (igraph::vcount(g) == 0) return(data.frame(gene=character(), degree=integer(), strength=numeric(), btw=numeric()))
    data.frame(
      gene = igraph::V(g)$name %||% as.character(seq_len(igraph::vcount(g))),
      degree = igraph::degree(g),
      strength = igraph::strength(g, weights = abs(igraph::E(g)$weight %||% 1)),
      btw = igraph::betweenness(g, directed=FALSE, weights = 1/pmax(1e-9, abs(igraph::E(g)$weight %||% 1)), normalized=FALSE),
      stringsAsFactors = FALSE
    )
  }

  gm1 <- global_metrics(g1); gm2 <- global_metrics(g2)
  gm1$Network <- label1; gm2$Network <- label2
  gmtab <- rbind(gm1, gm2)

  nm1 <- node_metrics(g1); nm1$Network <- label1
  nm2 <- node_metrics(g2); nm2$Network <- label2
  nodes_merged <- merge(nm1, nm2, by="gene", all=FALSE, suffixes=c(".1",".2"))

  gml <- reshape2::melt(gmtab, id.vars="Network")
  gml <- subset(gml, variable %in% c("Nodes","Edges","Density","Transitivity",
                                     "Avg_Degree","Avg_Path_LCC","Assortativity","Modularity"))
  p_global <- ggplot(gml, aes(x=variable, y=value, fill=Network)) +
    geom_col(position=position_dodge(width=0.7), width=0.6) +
    coord_flip() +
    labs(x=NULL, y=NULL, title="Global topology comparison") +
    theme_minimal(base_size=12)
  ggsave(file.path(outdir, "01_global_topology.png"), p_global, width=10, height=6, dpi=200)

  ecdf_long <- rbind(
    data.frame(Network=label1, degree=nm1$degree, strength=nm1$strength, btw=nm1$btw),
    data.frame(Network=label2, degree=nm2$degree, strength=nm2$strength, btw=nm2$btw)
  )
  plot_ecdf <- function(col, ttl) {
    ggplot(ecdf_long, aes(x=log1p(.data[[col]]), color=Network)) +
      stat_ecdf(geom="step", linewidth=1) +
      labs(x=paste0("log1p(",ttl,")"), y="ECDF", title=paste0("ECDF of ", ttl)) +
      theme_minimal(base_size=12)
  }
  ggsave(file.path(outdir,"02_ecdf_degree.png"),      plot_ecdf("degree","degree"),       width=6, height=4, dpi=200)
  ggsave(file.path(outdir,"03_ecdf_strength.png"),    plot_ecdf("strength","strength"),   width=6, height=4, dpi=200)
  ggsave(file.path(outdir,"04_ecdf_betweenness.png"), plot_ecdf("btw","betweenness"),     width=6, height=4, dpi=200)

  if (nrow(nodes_merged) > 0) {
    scat <- function(xcol, ycol, ttl) {
      ggplot(nodes_merged, aes(x=log1p(.data[[xcol]]), y=log1p(.data[[ycol]]))) +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="gray60") +
        geom_point(alpha=0.5, size=1.2) +
        labs(x=paste0(label1," log1p(",ttl,")"), y=paste0(label2," log1p(",ttl,")"),
             title=paste0(ttl, ": ", label1, " vs ", label2)) +
        theme_minimal(base_size=12)
    }
    ggsave(file.path(outdir,"05_scatter_degree.png"),   scat("degree.1","degree.2","Degree"),   width=6, height=5, dpi=200)
    ggsave(file.path(outdir,"06_scatter_strength.png"), scat("strength.1","strength.2","Strength"), width=6, height=5, dpi=200)
  }

  # Module correspondence (Jaccard) — only if both graphs non-empty
  jaccard_modules <- function(gA, gB) {
    clA <- igraph::cluster_louvain(gA, weights = abs(igraph::E(gA)$weight %||% 1))
    clB <- igraph::cluster_louvain(gB, weights = abs(igraph::E(gB)$weight %||% 1))
    setsA <- split(igraph::V(gA)$name, igraph::membership(clA))
    setsB <- split(igraph::V(gB)$name, igraph::membership(clB))
    A <- names(setsA); B <- names(setsB)
    J <- matrix(0, nrow=length(A), ncol=length(B), dimnames=list(paste0("A",A), paste0("B",B)))
    for (i in seq_along(setsA)) for (j in seq_along(setsB)) {
      inter <- length(intersect(setsA[[i]], setsB[[j]]))
      uni   <- length(union(setsA[[i]], setsB[[j]]))
      J[i,j] <- if (uni>0) inter/uni else 0
    }
    J
  }
  if (igraph::vcount(g1)>0 && igraph::vcount(g2)>0 && igraph::ecount(g1)>0 && igraph::ecount(g2)>0) {
    J <- jaccard_modules(g1, g2)
    pheatmap::pheatmap(J, cluster_rows = TRUE, cluster_cols = TRUE,
                       main = paste("Module overlap (Jaccard):", label1, "vs", label2),
                       filename = file.path(outdir, "07_module_overlap_jaccard.png"),
                       width = 8, height = 6)
  }

  edge_key <- function(g) {
    if (igraph::ecount(g)==0) return(character())
    el <- igraph::as_edgelist(g, names=TRUE)
    apply(el, 1, function(z) paste(sort(z), collapse="__"))
  }
  e1 <- edge_key(g1); e2 <- edge_key(g2)
  edge_df <- data.frame(
    Class = c(paste0("Overlap (",label1,"∩",label2,")"), paste0("Only ",label1), paste0("Only ",label2)),
    Count = c(length(intersect(e1,e2)), length(setdiff(e1,e2)), length(setdiff(e2,e1))),
    stringsAsFactors = FALSE
  )
  p_edges <- ggplot(edge_df, aes(x=Class, y=Count, fill=Class)) +
    geom_col(width=0.6) + theme_minimal(base_size=12) +
    guides(fill="none") + labs(x=NULL, y="Edges", title="Edge overlap / rewiring")
  ggsave(file.path(outdir, "08_edge_overlap.png"), p_edges, width=7, height=4, dpi=200)

  invisible(list(global=gmtab, nodes_merged=nodes_merged))
}

## =========================================================
## Run DCN (cached)
## =========================================================
dcn_params <- list(
  alpha = 0.01,
  rho = NULL,
  t_percentile = 0.90,
  min_overlap = 3,
  require_both_above_t = FALSE,
  diff_alpha = 0.01,
  diff_t = 0.3,
  fdr_method = "BH",
  diff_fdr_method = "BH"
)

res_EC_vs_Healthy <- gd2
res_EC_vs_Dec <- gd1

## legacy filenames (optional)
saveRDS(res_EC_vs_Healthy$graph_group1,     file = "graph_groupEC_rho06.rds")
saveRDS(res_EC_vs_Healthy$graph_group2,     file = "graph_groupHealthy_rho06.rds")
saveRDS(res_EC_vs_Healthy$graph_difference, file = "graph_difference_ECvsHealthy.rds")

saveRDS(res_EC_vs_Dec$graph_group1,         file = "graph_groupEC_rho06_ECvsDecidua.rds")
saveRDS(res_EC_vs_Dec$graph_group2,         file = "graph_groupDecidua_rho06.rds")
saveRDS(res_EC_vs_Dec$graph_difference,     file = "graph_difference_ECvsDecidua.rds")

## =========================================================
## Analyze graphs + compare
## =========================================================
gEC  <- res_EC_vs_Healthy$graph_group1
gH   <- res_EC_vs_Healthy$graph_group2
gDec <- res_EC_vs_Dec$graph_group2

res_EC  <- analyze_one_graph(g1,  label="EC",      outdir="DCN_outputs",
                            min_module_size=20, top_terms_per_module=10, top_n_nodes=20)
res_H   <- analyze_one_graph(g3,   label="Healthy", outdir="DCN_outputs",
                            min_module_size=20, top_terms_per_module=10, top_n_nodes=20)
res_Dec <- analyze_one_graph(g2, label="Decidua", outdir="DCN_outputs",
                            min_module_size=20, top_terms_per_module=10, top_n_nodes=20)

compare_networks_plots(gEC, gDec, label1="EC", label2="Decidua", outdir="EC_vs_Decidua_compare")
compare_networks_plots(gEC, gH,   label1="EC", label2="Healthy", outdir="EC_vs_Healthy_compare")
compare_networks_plots(gDec,gH,   label1="Decidua", label2="Healthy", outdir="Decidua_vs_Healthy_compare")

message("DONE.")
