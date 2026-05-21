###############################################################################
## 0) Global configuration (single source of truth for paths)
###############################################################################

FILEDIR <- file.path("network", "network_06", "corr06")  # <<<<<< requested
dir.create(FILEDIR, showWarnings = FALSE, recursive = TRUE)

# Output roots (all under FILEDIR)
OUT_ANALYSIS   <- file.path(FILEDIR, "DCN_outputs")
OUT_COMPARE    <- file.path(FILEDIR, "NET_compare")
OUT_DIFF       <- file.path(FILEDIR, "DCN_diff_outputs")

dir.create(OUT_ANALYSIS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_COMPARE,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIFF,     showWarnings = FALSE, recursive = TRUE)

###############################################################################
## 1) Dependency checks (do this once; avoid repeated library() calls)
###############################################################################

need_pkg <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- if (bioc) {
      paste0("Missing package '", pkg, "'. Install via BiocManager::install('", pkg, "').")
    } else {
      paste0("Missing package '", pkg, "'. Install via install.packages('", pkg, "').")
    }
    stop(msg, call. = FALSE)
  }
}

# Core
need_pkg("igraph")
need_pkg("ggplot2")
need_pkg("reshape2")
need_pkg("pheatmap")

# Excel output
need_pkg("openxlsx")

# Enrichment
need_pkg("clusterProfiler", bioc = TRUE)
need_pkg("org.Hs.eg.db", bioc = TRUE)

# Optional visualization
HAS_CIRCLIZE <- requireNamespace("circlize", quietly = TRUE)

###############################################################################
## 2) Utilities (centralized helpers; avoid duplicates everywhere)
###############################################################################

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- robust numeric vector conversion ----
as_numeric_safe <- function(x, default = 0) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  x <- suppressWarnings(as.numeric(x))
  x[!is.finite(x)] <- default
  x
}

# ---- ensure igraph edge weights exist and are numeric ----
sanitize_graph_weights <- function(g, default = 0) {
  stopifnot(inherits(g, "igraph"))
  if (igraph::ecount(g) == 0) return(g)

  w <- igraph::E(g)$weight
  if (is.null(w)) {
    igraph::E(g)$weight <- rep(default, igraph::ecount(g))
    return(g)
  }

  w <- as_numeric_safe(w, default = default)
  if (length(w) != igraph::ecount(g)) w <- rep(default, igraph::ecount(g))

  igraph::E(g)$weight <- w
  g
}

# ---- largest connected component after removing isolates ----
safe_lcc <- function(g) {
  stopifnot(inherits(g, "igraph"))
  if (igraph::vcount(g) == 0) return(g)

  g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
  if (igraph::ecount(g) == 0 || igraph::vcount(g) < 2) return(g)

  cmp <- igraph::components(g)

  memb <- cmp$membership
  if (!is.atomic(memb)) memb <- unlist(memb, use.names = FALSE)
  memb <- as.integer(memb)

  cs <- cmp$csize
  if (!is.atomic(cs)) cs <- unlist(cs, use.names = FALSE)
  cs <- as.numeric(cs)

  lcc_id <- which.max(cs)
  igraph::induced_subgraph(g, vids = which(memb == lcc_id))
}

# ---- edge key for undirected overlap measures ----
edge_key <- function(g) {
  if (igraph::ecount(g) == 0) return(character())
  el <- igraph::as_edgelist(g, names = TRUE)
  apply(el, 1, function(z) paste(sort(z), collapse = "__"))
}

# ---- consistent plotting save ----
save_png <- function(path, plot_obj, width, height, dpi = 200) {
  ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, dpi = dpi)
}


###############################################################################
## 3) Analyze a single network (single canonical function)
###############################################################################

analyze_one_graph <- function(g,
                              label,
                              outdir = OUT_ANALYSIS,
                              orgDb = org.Hs.eg.db::org.Hs.eg.db,
                              keyType = "SYMBOL",
                              min_module_size = 20,
                              top_terms_per_module = 10,
                              top_n_nodes = 20,
                              make_circlize = TRUE,
                              make_module_metagraph = TRUE) {
  stopifnot(inherits(g, "igraph"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Per-graph output dirs
  gdir      <- file.path(outdir, label)
  plots_dir <- file.path(gdir, "module_plots")
  dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  # Sanitize weights early
  g <- sanitize_graph_weights(g, default = 0)

  # ---- Clean: drop isolates; keep LCC ----
  g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
  if (igraph::ecount(g) == 0 || igraph::vcount(g) < 2) {
    warning(sprintf("[%s] No edges after cleaning.", label))
    return(invisible(NULL))
  }
  g <- safe_lcc(g)

  # ---- Communities: Louvain on |weight| ----
  igraph::E(g)$w_abs <- abs(igraph::E(g)$weight %||% 0)
  cl <- igraph::cluster_louvain(g, weights = igraph::E(g)$w_abs)

  mb <- igraph::membership(cl)
  if (!is.atomic(mb)) mb <- unlist(mb, use.names = TRUE)

  # Map membership to vertices by name
  igraph::V(g)$module <- as.integer(mb[igraph::V(g)$name])
  mod_sizes <- igraph::sizes(cl)
  modules   <- split(names(mb), mb)

  # ---- Node-level metrics (for hubs & betweenness) ----
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

  # ---- Graph-level topology metrics ----
  comps <- igraph::components(g)
  g_lcc <- safe_lcc(g)

  topo_metrics <- data.frame(
    Label = label,
    Nodes = igraph::vcount(g),
    Edges = igraph::ecount(g),
    Components = comps$no,
    LCC_Nodes = igraph::vcount(g_lcc),
    LCC_Edges = igraph::ecount(g_lcc),
    Density = igraph::edge_density(g, loops = FALSE),
    Avg_Degree = mean(as.numeric(deg)),
    Transitivity = igraph::transitivity(g, type = "global", isolates = "zero"),
    Assortativity_Degree = tryCatch(igraph::assortativity_degree(g, directed = FALSE),
                                    error = function(e) NA_real_),
    Vertex_Connectivity = tryCatch(igraph::vertex_connectivity(g), error = function(e) NA_real_),
    Edge_Connectivity   = tryCatch(igraph::edge_connectivity(g),   error = function(e) NA_real_),
    Avg_Path_Length_LCC = tryCatch(igraph::average.path.length(g_lcc, directed = FALSE),
                                   error = function(e) NA_real_),
    Modularity_Louvain  = tryCatch(igraph::modularity(g, membership = igraph::V(g)$module,
                                                     weights = igraph::E(g)$w_abs),
                                   error = function(e) NA_real_),
    stringsAsFactors = FALSE
  )

  # ---- Module ranking (topology-based composite score) ----
  z <- function(x) {
    x <- as.numeric(x)
    if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  }

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
      stringsAsFactors = FALSE
    )
  })
  module_rank_df <- do.call(rbind, mod_rank_list)
  module_rank_df$PassSize <- module_rank_df$N >= min_module_size

  module_rank_df$Score <- with(module_rank_df,
    1.2 * z(Density) +
    1.2 * z(MeanAbsStrength) +
    0.8 * z(MaxBetweenness) +
    0.4 * z(HubConcentration) +
    0.6 * z(Clustering)
  )

  # Order: size-pass first, then score
  module_rank_df <- module_rank_df[order(-as.integer(module_rank_df$PassSize),
                                         -module_rank_df$Score,
                                         -module_rank_df$N,
                                         -module_rank_df$Density), , drop = FALSE]

  # Lightweight module typing for interpretation/figures
  dens_q75 <- stats::quantile(module_rank_df$Density[module_rank_df$PassSize], 0.75, na.rm = TRUE)
  module_rank_df$ModuleType <- with(module_rank_df, ifelse(
    !PassSize, "Too small / unstable",
    ifelse(HubConcentration >= 0.25 & N <= 60, "Hub-dominated driver-like core",
           ifelse(is.finite(dens_q75) & Density >= dens_q75, "Dense coherent program",
                  "Moderate program"))
  ))

  # ---- Per-module graph plots (PNG): size by degree; label top-5 hubs ----
  pal <- grDevices::hcl.colors(length(module_ids), "Dark 3")
  names(pal) <- as.character(module_ids)

  for (m in module_ids) {
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
         main = sprintf("%s — Module %s (n=%d) — node size=degree; labels=top5 hubs",
                        label, m, igraph::vcount(gm)))
    grDevices::dev.off()
  }

  # ---- Optional: circlize cross-module chord plot ----
  circlize_png <- NA_character_
  if (make_circlize && HAS_CIRCLIZE) {
    # Build cross-module edge counts (undirected)
    edf <- as.data.frame(igraph::as_edgelist(g, names = TRUE), stringsAsFactors = FALSE)
    colnames(edf) <- c("v1","v2")
    edf$m1 <- igraph::V(g)$module[match(edf$v1, igraph::V(g)$name)]
    edf$m2 <- igraph::V(g)$module[match(edf$v2, igraph::V(g)$name)]
    cross <- edf[edf$m1 != edf$m2, , drop = FALSE]

    if (nrow(cross) > 0) {
      key <- apply(cross[, c("m1","m2")], 1, function(z) paste(sort(z), collapse = "_"))
      w_counts <- as.data.frame(table(key), stringsAsFactors = FALSE)
      parts <- strsplit(w_counts$key, "_", fixed = TRUE)

      links_df <- data.frame(
        from  = vapply(parts, `[`, "", 1),
        to    = vapply(parts, `[`, "", 2),
        value = as.numeric(w_counts$Freq),
        stringsAsFactors = FALSE
      )

      mod_sizes_vec <- setNames(as.integer(mod_sizes), names(mod_sizes))
      sector_ord <- as.character(names(sort(mod_sizes_vec, decreasing = TRUE)))
      grid_col <- setNames(grDevices::hcl.colors(length(sector_ord), "Dark 3"), sector_ord)

      circlize_png <- file.path(gdir, sprintf("%s_modules_circlize.png", label))
      grDevices::png(circlize_png, width = 1600, height = 1600, res = 200)
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
    }
  }

  # ---- Optional: module meta-graph (modules as nodes) ----
  module_metagraph_png <- NA_character_
  g_module <- NULL

  if (make_module_metagraph) {
    build_module_metagraph <- function(g, module_rank_df, out_png) {
      ed <- as.data.frame(igraph::as_edgelist(g, names = TRUE), stringsAsFactors = FALSE)
      colnames(ed) <- c("v1","v2")
      ed$w <- igraph::E(g)$weight %||% 0
      ed$m1 <- igraph::V(g)$module[match(ed$v1, igraph::V(g)$name)]
      ed$m2 <- igraph::V(g)$module[match(ed$v2, igraph::V(g)$name)]
      ed <- ed[ed$m1 != ed$m2, , drop = FALSE]
      if (!nrow(ed)) return(NULL)

      key <- apply(ed[, c("m1","m2")], 1, function(z) paste(sort(z), collapse = "__"))
      ed$key <- key

      # Prefer dplyr if available; otherwise do base aggregation
      if (requireNamespace("dplyr", quietly = TRUE)) {
        agg <- dplyr::as_tibble(ed) |>
          dplyr::group_by(.data$key) |>
          dplyr::summarise(
            mA = min(.data$m1, .data$m2),
            mB = max(.data$m1, .data$m2),
            EdgeAbsSum = sum(abs(.data$w), na.rm = TRUE),
            MedianW = stats::median(.data$w, na.rm = TRUE),
            .groups = "drop"
          )
        agg <- as.data.frame(agg)
      } else {
        # Base fallback aggregation
        split_ed <- split(ed, ed$key)
        agg <- do.call(rbind, lapply(split_ed, function(df) {
          data.frame(
            mA = min(df$m1, df$m2),
            mB = max(df$m1, df$m2),
            EdgeAbsSum = sum(abs(df$w), na.rm = TRUE),
            MedianW = stats::median(df$w, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }))
        agg$key <- rownames(agg); rownames(agg) <- NULL
      }

      verts <- module_rank_df[, c("Module","N","ModuleType"), drop = FALSE]
      verts$name <- as.character(verts$Module)

      gM <- igraph::graph_from_data_frame(
        d = data.frame(from = as.character(agg$mA),
                       to   = as.character(agg$mB),
                       weight = agg$EdgeAbsSum,
                       medw   = agg$MedianW,
                       stringsAsFactors = FALSE),
        directed = FALSE,
        vertices = verts
      )

      # Node sizes ~ module size
      N <- igraph::V(gM)$N %||% 1
      N <- as.numeric(N); N[is.na(N)] <- 1
      igraph::V(gM)$size <- 6 + 18 * (N / max(N))

      # Node colors by module type
      type_levels <- unique(igraph::V(gM)$ModuleType)
      type_cols <- setNames(grDevices::hcl.colors(max(3, length(type_levels)), "Set 2")[seq_along(type_levels)],
                            type_levels)
      igraph::V(gM)$color <- type_cols[igraph::V(gM)$ModuleType]

      # Edge widths by sum(|w|)
      w <- igraph::E(gM)$weight %||% 1
      w <- as.numeric(w); w[is.na(w)] <- 0
      igraph::E(gM)$width <- 0.5 + 6 * (w / max(w))

      # Edge colors by median sign
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

      gM
    }

    module_metagraph_png <- file.path(gdir, sprintf("%s_module_metagraph.png", label))
    g_module <- build_module_metagraph(g, module_rank_df, module_metagraph_png)
  }

  # ---- Enrichment: GO:BP per module (clusterProfiler) ----
  go_enrichments <- lapply(names(modules), function(m) {
    genes <- modules[[m]]
    if (length(genes) < min_module_size) return(NULL)

    suppressWarnings(tryCatch(
      clusterProfiler::enrichGO(
        gene          = genes,
        OrgDb         = orgDb,
        keyType       = keyType,
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) NULL
    ))
  })
  names(go_enrichments) <- paste0("Module_", names(modules))

  make_enrichment_summary <- function(enrichments, top_n = 10) {
    out <- lapply(names(enrichments), function(mname) {
      enr <- enrichments[[mname]]
      if (is.null(enr)) return(NULL)
      df <- as.data.frame(enr)
      if (!nrow(df)) return(NULL)

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
      out <- data.frame(Module = integer(), Term = character(), Adj_P = numeric(),
                        ID = character(), P_value = numeric(), Q_value = numeric(),
                        GeneRatio = character(), BgRatio = character(), Count = integer(),
                        Genes = character(), stringsAsFactors = FALSE)
    } else {
      out <- out[order(out$Module, out$Adj_P, na.last = TRUE), , drop = FALSE]
      rownames(out) <- NULL
    }
    out
  }

  go_summary_df <- make_enrichment_summary(go_enrichments, top_n = top_terms_per_module)

  # ---- Excel output (single workbook; correct ordering) ----
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "Modules")
  openxlsx::writeData(wb, "Modules",
                      data.frame(Gene = names(mb), Module = as.integer(mb), stringsAsFactors = FALSE))

  openxlsx::addWorksheet(wb, "Module_Ranking")
  openxlsx::writeData(wb, "Module_Ranking", module_rank_df)
  try({
    openxlsx::addStyle(wb, "Module_Ranking", openxlsx::createStyle(textDecoration = "bold"),
                       rows = 1, cols = 1:ncol(module_rank_df), gridExpand = TRUE)
    openxlsx::freezePane(wb, "Module_Ranking", firstRow = TRUE)
  }, silent = TRUE)

  openxlsx::addWorksheet(wb, "TopHubs_Betweenness")
  openxlsx::writeData(wb, "TopHubs_Betweenness", "Top hubs (by |strength|)", startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "TopHubs_Betweenness", df_hubs, startRow = 2, startCol = 1)

  start_row_btw <- nrow(df_hubs) + 4
  openxlsx::writeData(wb, "TopHubs_Betweenness", "Top betweenness", startRow = start_row_btw, startCol = 1)
  openxlsx::writeData(wb, "TopHubs_Betweenness", df_btw, startRow = start_row_btw + 1, startCol = 1)

  openxlsx::addWorksheet(wb, "Graph_Topology")
  openxlsx::writeData(wb, "Graph_Topology", topo_metrics)

  openxlsx::addWorksheet(wb, "GO_BP_Summary")
  openxlsx::writeData(wb, "GO_BP_Summary", go_summary_df)
  try({
    openxlsx::addStyle(wb, "GO_BP_Summary", openxlsx::createStyle(textDecoration = "bold"),
                       rows = 1, cols = 1:ncol(go_summary_df), gridExpand = TRUE)
    openxlsx::freezePane(wb, "GO_BP_Summary", firstRow = TRUE)
  }, silent = TRUE)

  # Optional: per-module GO result sheets (can be large)
  for (m in names(go_enrichments)) {
    enr <- go_enrichments[[m]]
    if (!is.null(enr)) {
      sh <- substr(paste0("GO_", m), 1, 31)  # Excel sheet name limit
      openxlsx::addWorksheet(wb, sh)
      openxlsx::writeData(wb, sh, as.data.frame(enr))
    }
  }

  xlsx_path <- file.path(gdir, sprintf("%s_results.xlsx", label))
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  # ---- Save objects for reproducibility ----
  saveRDS(g,             file = file.path(gdir, sprintf("%s_graph_cleaned.rds", label)))
  saveRDS(cl,            file = file.path(gdir, sprintf("%s_louvain_object.rds", label)))
  saveRDS(go_enrichments,file = file.path(gdir, sprintf("%s_GO_enrichments.rds", label)))
  if (!is.null(g_module)) saveRDS(g_module, file = file.path(gdir, sprintf("%s_module_metagraph.rds", label)))

  message(sprintf("[%s] Done. Excel: %s | Module plots: %s",
                  label, xlsx_path, plots_dir))

  invisible(list(
    graph = g,
    membership = mb,
    sizes = mod_sizes,
    hubs = df_hubs,
    betweenness = df_btw,
    topo = topo_metrics,
    module_ranking = module_rank_df,
    go_enrichments = go_enrichments,
    go_summary = go_summary_df,
    excel = xlsx_path,
    plots_dir = plots_dir,
    circlize_png = circlize_png,
    module_metagraph_png = module_metagraph_png
  ))
}


###############################################################################
## 4) Differential DCN analysis — READS EXISTING GRAPHS (RDS)
## Outputs:
##  1) Topology barplots for ORIGINAL graphs (EC/Healthy/Decidua)
##  2) Rewiring burden barplots for DIFFERENTIAL graphs (3 contrasts)
##  3) Hallmark (MSigDB H) enrichment heatmap across differential modules
##  4) Top rewired genes per contrast (tables + optional plots)
###############################################################################

# ---- Optional packages for this section ----
HAS_TIDYR    <- requireNamespace("tidyr", quietly = TRUE)
HAS_MSIGDBR  <- requireNamespace("msigdbr", quietly = TRUE)
HAS_CP       <- requireNamespace("clusterProfiler", quietly = TRUE)
HAS_PHEATMAP <- requireNamespace("pheatmap", quietly = TRUE)

if (!HAS_TIDYR) message("[note] install.packages('tidyr') for pivot_longer/wider convenience.")
if (!HAS_MSIGDBR) message("[note] install.packages('msigdbr') to enable Hallmark enrichment.")
if (!HAS_CP) message("[note] BiocManager::install('clusterProfiler') to enable Hallmark enrichment.")
if (!HAS_PHEATMAP) message("[note] install.packages('pheatmap') to enable heatmap output.")

dir.create(OUT_DIFF, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 6.1) Helper: safe graph reader (validates igraph + sanitizes weights)
# ---------------------------------------------------------------------------
read_graph_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  g <- readRDS(path)
  if (!inherits(g, "igraph")) stop("Not an igraph object: ", path, call. = FALSE)
  sanitize_graph_weights(g, default = 0)
}

# ---------------------------------------------------------------------------
# 6.2) Configure RDS inputs (EDIT FILENAMES HERE ONCE)
#      These are *examples*. Keep them all under FILEDIR.
# ---------------------------------------------------------------------------

# Original (condition) graphs: weight = rho/correlation
f_orig <- list(
  EC      = file.path(FILEDIR, "DCN_graphs_08c4a5f8534f4497_EC_group1.rds"),
  Healthy = file.path(FILEDIR, "DCN_graphs_08c4a5f8534f4497_Healthy_group2.rds"),
  Decidua =  file.path(FILEDIR, "DCN_graphs_98b61808d296fb76_Decidua_group1.rds")
)


# Differential (contrast) graphs: weight = delta correlation (case - ctrl)
f_diff <- list(
  EC_vs_Healthy      = file.path(FILEDIR, "DCN_graphs_08c4a5f8534f4497_EC_vs_Healthy_diff.rds"),
  EC_vs_Decidua      = file.path(FILEDIR, "DCN_graphs_a8e99029cf08c086_EC_vs_Decidua_diff.rds"),
  Decidua_vs_Healthy = file.path(FILEDIR, "DCN_graphs_98b61808d296fb76_Decidua_vs_Healthy_diff.rds")  # rename if needed
)

g_orig <- lapply(f_orig, read_graph_safe)
g_diff <- lapply(f_diff, read_graph_safe)

# Fail fast if required graphs missing
if (any(vapply(g_orig, is.null, logical(1)))) {
  stop("Missing original graphs:\n",
       paste(names(g_orig)[vapply(g_orig, is.null, logical(1))], collapse = ", "),
       "\nCheck f_orig paths.", call. = FALSE)
}
if (any(vapply(g_diff, is.null, logical(1)))) {
  stop("Missing differential graphs:\n",
       paste(names(g_diff)[vapply(g_diff, is.null, logical(1))], collapse = ", "),
       "\nCheck f_diff paths.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# 6.3) Global topology metrics for ORIGINAL graphs
#      (Used for barplots + CSV)
# ---------------------------------------------------------------------------
global_topology <- function(g) {
  stopifnot(inherits(g, "igraph"))
  Nodes <- igraph::vcount(g)
  Edges <- igraph::ecount(g)

  if (Nodes < 2 || Edges < 1) {
    return(data.frame(
      Nodes = Nodes, Edges = Edges,
      Density = NA_real_,
      Transitivity = NA_real_,
      Avg_Degree = NA_real_,
      Mean_Distance_LCC = NA_real_,
      Assortativity = NA_real_,
      Modularity = NA_real_,
      LCC_Nodes = NA_real_,
      LCC_Edges = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  g_lcc <- safe_lcc(g)
  deg <- igraph::degree(g)

  md <- tryCatch(igraph::mean_distance(g_lcc, directed = FALSE),
                 error = function(e) NA_real_)

  wabs <- abs(igraph::E(g)$weight %||% 0)
  mod <- tryCatch({
    cl <- igraph::cluster_louvain(g, weights = wabs)
    igraph::modularity(g, igraph::membership(cl), weights = wabs)
  }, error = function(e) NA_real_)

  data.frame(
    Nodes = Nodes,
    Edges = Edges,
    Density = igraph::edge_density(g, loops = FALSE),
    Transitivity = igraph::transitivity(g, type = "global", isolates = "zero"),
    Avg_Degree = mean(as.numeric(deg)),
    Mean_Distance_LCC = md,
    Assortativity = tryCatch(igraph::assortativity_degree(g, directed = FALSE),
                             error = function(e) NA_real_),
    Modularity = mod,
    LCC_Nodes = igraph::vcount(g_lcc),
    LCC_Edges = igraph::ecount(g_lcc),
    stringsAsFactors = FALSE
  )
}

topo_df <- do.call(rbind, lapply(names(g_orig), function(nm) {
  df <- global_topology(g_orig[[nm]])
  df$Graph <- nm
  df
}))

# Save CSV
write.csv(topo_df, file.path(OUT_DIFF, "original_topology_metrics.csv"),
          row.names = FALSE, quote = FALSE)

# Barplots (facet per metric)
if (HAS_TIDYR) {
  topo_long <- tidyr::pivot_longer(
    topo_df,
    cols = c("Nodes","Edges","Density","Transitivity","Avg_Degree",
             "Mean_Distance_LCC","Assortativity","Modularity"),
    names_to = "Metric",
    values_to = "Value"
  )

  p_topo <- ggplot2::ggplot(topo_long, ggplot2::aes(x = Graph, y = Value, fill = Graph)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 3) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   legend.position = "none") +
    ggplot2::labs(title = "Original graphs — topology metrics", x = NULL, y = NULL)

  save_png(file.path(OUT_DIFF, "01_original_topology_metrics_barplots.png"),
           p_topo, width = 12, height = 7, dpi = 200)
}

# ---------------------------------------------------------------------------
# 6.4) Rewiring burden metrics for DIFFERENTIAL graphs
#      Assumes E(gd)$weight = delta (case - ctrl)
# ---------------------------------------------------------------------------
rewiring_burden <- function(gd) {
  stopifnot(inherits(gd, "igraph"))

  if (igraph::vcount(gd) < 2 || igraph::ecount(gd) < 1) {
    return(data.frame(
      DiffEdges = igraph::ecount(gd),
      SumAbsDelta = 0,
      MeanAbsDelta = NA_real_,
      FracPos = NA_real_,
      Balance = NA_real_,   # (sum_pos - sum_neg) / (sum_pos + sum_neg)
      stringsAsFactors = FALSE
    ))
  }

  d <- as.numeric(igraph::E(gd)$weight %||% numeric(0))
  d <- d[is.finite(d)]
  if (!length(d)) {
    return(data.frame(DiffEdges = igraph::ecount(gd), SumAbsDelta = 0,
                      MeanAbsDelta = NA_real_, FracPos = NA_real_, Balance = NA_real_,
                      stringsAsFactors = FALSE))
  }

  sum_pos <- sum(d[d > 0], na.rm = TRUE)
  sum_neg <- sum(abs(d[d < 0]), na.rm = TRUE)
  denom <- sum_pos + sum_neg
  bal <- if (denom > 0) (sum_pos - sum_neg) / denom else NA_real_

  data.frame(
    DiffEdges = length(d),
    SumAbsDelta = sum(abs(d), na.rm = TRUE),
    MeanAbsDelta = mean(abs(d), na.rm = TRUE),
    FracPos = mean(d > 0, na.rm = TRUE),
    Balance = bal,
    stringsAsFactors = FALSE
  )
}

rew_df <- do.call(rbind, lapply(names(g_diff), function(nm) {
  df <- rewiring_burden(g_diff[[nm]])
  df$Contrast <- nm
  df
}))

write.csv(rew_df, file.path(OUT_DIFF, "differential_rewiring_burden.csv"),
          row.names = FALSE, quote = FALSE)

if (HAS_TIDYR) {
  rew_long <- tidyr::pivot_longer(
    rew_df,
    cols = c("DiffEdges","SumAbsDelta","MeanAbsDelta","FracPos","Balance"),
    names_to = "Metric",
    values_to = "Value"
  )

  p_rew <- ggplot2::ggplot(rew_long, ggplot2::aes(x = Contrast, y = Value, fill = Contrast)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~Metric, scales = "free_y", ncol = 3) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   legend.position = "none") +
    ggplot2::labs(title = "Differential graphs — rewiring burden metrics", x = NULL, y = NULL)

  save_png(file.path(OUT_DIFF, "02_differential_rewiring_burden_barplots.png"),
           p_rew, width = 12, height = 7, dpi = 200)
}

# ---------------------------------------------------------------------------
# 6.5) Top rewired genes per contrast
#      Node rewiring score = sum(|delta|) over incident edges
# ---------------------------------------------------------------------------
node_rewiring_scores <- function(gd) {
  stopifnot(inherits(gd, "igraph"))

  vnames <- igraph::V(gd)$name
  if (is.null(vnames)) vnames <- as.character(seq_len(igraph::vcount(gd)))

  if (igraph::vcount(gd) == 0) {
    return(data.frame(Gene = character(), R_strength = numeric(), R_degree = integer(),
                      stringsAsFactors = FALSE))
  }
  if (igraph::ecount(gd) == 0) {
    return(data.frame(Gene = vnames, R_strength = 0, R_degree = 0L,
                      stringsAsFactors = FALSE))
  }

  w <- as.numeric(igraph::E(gd)$weight %||% rep(0, igraph::ecount(gd)))
  w[!is.finite(w)] <- 0
  wabs <- abs(w)

  el <- igraph::as_edgelist(gd, names = TRUE)
  df <- data.frame(v1 = el[,1], v2 = el[,2], wabs = wabs, stringsAsFactors = FALSE)

  # Sum abs weights incident to each endpoint
  rs <- rbind(
    data.frame(Gene = df$v1, wabs = df$wabs),
    data.frame(Gene = df$v2, wabs = df$wabs)
  )
  agg <- aggregate(rs$wabs, by = list(Gene = rs$Gene), FUN = sum)
  names(agg)[2] <- "R_strength"

  deg_df <- data.frame(Gene = vnames,
                       R_degree = as.integer(igraph::degree(gd)),
                       stringsAsFactors = FALSE)

  out <- merge(deg_df, agg, by = "Gene", all.x = TRUE)
  out$R_strength[is.na(out$R_strength)] <- 0
  out <- out[order(-out$R_strength, -out$R_degree), , drop = FALSE]
  rownames(out) <- NULL
  out
}

TOPN <- 50
top_rewired_list <- lapply(names(g_diff), function(nm) {
  df <- node_rewiring_scores(g_diff[[nm]])
  df$Contrast <- nm
  head(df, TOPN)
})
top_rewired_df <- do.call(rbind, top_rewired_list)

write.csv(top_rewired_df,
          file.path(OUT_DIFF, sprintf("top_rewired_genes_top%d_per_contrast.csv", TOPN)),
          row.names = FALSE, quote = FALSE)

# ---------------------------------------------------------------------------
# 6.6) Hallmark enrichment heatmap across differential modules
#      - Detect modules in each differential graph (Louvain on |delta|)
#      - Enrichment with msigdbr + clusterProfiler::enricher
#      - Heatmap score = -log10(FDR), capped for display
# ---------------------------------------------------------------------------
hallmark_enrich_by_module <- function(gd, min_module_size = 20) {
  if (!(HAS_MSIGDBR && HAS_CP)) return(NULL)

  gd <- sanitize_graph_weights(gd, default = 0)
  if (igraph::vcount(gd) < 2 || igraph::ecount(gd) < 1) return(NULL)

  # Module detection on |delta|
  wabs <- abs(igraph::E(gd)$weight %||% 0)
  cl <- igraph::cluster_louvain(gd, weights = wabs)
  mb <- igraph::membership(cl)
  if (!is.atomic(mb)) mb <- unlist(mb, use.names = TRUE)

  modules <- split(names(mb), mb)

  # MSigDB Hallmark TERM2GENE
  hall_tbl <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  term2gene <- hall_tbl[, c("gs_name","gene_symbol")]
  colnames(term2gene) <- c("term","gene")

  out <- lapply(names(modules), function(m) {
    genes <- modules[[m]]
    if (length(genes) < min_module_size) return(NULL)

    enr <- tryCatch(
      clusterProfiler::enricher(
        gene = genes,
        TERM2GENE = term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      ),
      error = function(e) NULL
    )
    if (is.null(enr)) return(NULL)

    df <- as.data.frame(enr)
    if (!nrow(df)) return(NULL)

    df$Module <- as.integer(m)
    df
  })

  out <- do.call(rbind, out)
  if (is.null(out) || !nrow(out)) return(NULL)
  out
}

if (HAS_MSIGDBR && HAS_CP && HAS_PHEATMAP && HAS_TIDYR) {
  hall_list <- lapply(names(g_diff), function(nm) {
    df <- hallmark_enrich_by_module(g_diff[[nm]], min_module_size = 20)
    if (is.null(df) || !nrow(df)) return(NULL)

    df$Contrast <- nm
    df$ColID <- paste0(nm, "::M", df$Module)
    df$NegLogFDR <- -log10(df$p.adjust %||% df$pvalue)

    df[, c("ColID","Contrast","Module","ID","Description","p.adjust","NegLogFDR","Count","geneID")]
  })
  hall_df <- do.call(rbind, hall_list)

  if (!is.null(hall_df) && nrow(hall_df) > 0) {
    # Keep top terms per column to control heatmap size
    TOP_TERMS <- 10
    hall_df2 <- hall_df[order(hall_df$ColID, hall_df$p.adjust, hall_df$NegLogFDR, decreasing = FALSE), ]
    hall_df2 <- do.call(rbind, lapply(split(hall_df2, hall_df2$ColID), function(x) head(x, TOP_TERMS)))

    # Term-by-(contrast::module) matrix
    mat_df <- aggregate(hall_df2$NegLogFDR,
                        by = list(Term = hall_df2$Description, ColID = hall_df2$ColID),
                        FUN = max, na.rm = TRUE)
    names(mat_df)[3] <- "Score"

    mat_wide <- tidyr::pivot_wider(mat_df, names_from = "ColID", values_from = "Score", values_fill = 0)
    mat <- as.data.frame(mat_wide)
    rownames(mat) <- mat$Term
    mat$Term <- NULL
    mat <- as.matrix(mat)

    # cap for display
    mat_cap <- pmin(mat, 10)

    ann_col <- data.frame(
      Contrast = sub("::M\\d+$", "", colnames(mat_cap))
    )
    rownames(ann_col) <- colnames(mat_cap)

    pheatmap::pheatmap(
      mat_cap,
      annotation_col = ann_col,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      main = "Hallmark enrichment across differential modules (score = -log10(FDR), capped at 10)",
      filename = file.path(OUT_DIFF, "03_hallmark_enrichment_heatmap_diff_modules.png"),
      width = 12,
      height = 8
    )

    write.csv(hall_df, file.path(OUT_DIFF, "hallmark_enrichment_by_diff_module_full.csv"),
              row.names = FALSE, quote = FALSE)
  } else {
    message("[note] No Hallmark enrichment results (modules too small or gene symbols mismatch).")
  }
} else {
  message("[skip] Hallmark heatmap skipped (need: tidyr + msigdbr + clusterProfiler + pheatmap).")
}

message("[done] Differential outputs written to: ", OUT_DIFF)
