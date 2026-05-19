compare_networks_plots <- function(g1, g2, label1="Net1", label2="Net2", outdir="NET_COMPARE") {
  stopifnot(inherits(g1, "igraph"), inherits(g2, "igraph"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  suppressPackageStartupMessages({
    library(igraph); library(ggplot2); library(reshape2); library(pheatmap)
  })

  # ---------- helpers ----------
  safe_lcc <- function(g) {
    g <- delete_vertices(g, degree(g)==0)
    if (ecount(g)==0 || vcount(g)<2) return(g)
    c <- components(g); induced_subgraph(g, which(c$membership==which.max(c$csize)))
  }
  global_metrics <- function(g) {
    if (vcount(g)==0 || ecount(g)==0) {
      return(data.frame(Nodes=vcount(g), Edges=ecount(g), Density=NA, Transitivity=NA,
                        Avg_Degree=NA, Avg_Path_LCC=NA, Assortativity=NA,
                        VertexConn=NA, EdgeConn=NA, Modularity=NA))
    }
    mb_dummy <- rep(1, vcount(g))  # placeholder for modularity if no communities run
    data.frame(
      Density = edge_density(g),
      Transitivity = transitivity(g, type="global", isolates="zero"),
      Assortativity = assortativity_degree(g, directed=FALSE),
      Modularity = tryCatch({
        cl <- cluster_louvain(g, weights = abs(E(g)$weight))
        modularity(g, membership(cl), weights = abs(E(g)$weight))
      }, error=function(e) NA)
    )
  }
  node_metrics <- function(g) {
    if (vcount(g)==0) return(data.frame(gene=character(), degree=integer(), strength=numeric(), btw=numeric()))
    data.frame(
      gene = V(g)$name %||% as.character(seq_len(vcount(g))),
      degree = degree(g),
      strength = strength(g, weights = abs(E(g)$weight %||% 1)),
      btw = betweenness(g, directed=FALSE,
                        weights = 1/pmax(1e-9, abs(E(g)$weight %||% 1)), normalized=FALSE)
    )
  }
  `%||%` <- function(a,b) if (!is.null(a)) a else b

  # ---------- compute ----------
  gm1 <- global_metrics(g1); gm2 <- global_metrics(g2)
  gm1$Network <- label1; gm2$Network <- label2
  gmtab <- rbind(gm1, gm2)

  nm1 <- node_metrics(g1); nm1$Network <- label1
  nm2 <- node_metrics(g2); nm2$Network <- label2
  nodes_merged <- merge(nm1, nm2, by="gene", all=FALSE, suffixes=c(".1",".2"))

  # ---------- 1) Global dashboard ----------
  # ---------- 1) Global dashboard ----------
  gml <- reshape2::melt(gmtab, id.vars = "Network")
  
  keep_vars <- c("Nodes","Edges","Density","Transitivity","Avg_Degree","Avg_Path_LCC","Assortativity","Modularity")
  
  # robust filter (works even if gml is empty)
  if (!("variable" %in% colnames(gml)) || nrow(gml) == 0) {
    warning("Global metrics melt produced 0 rows; skipping global topology plot.")
  } else {
    gml <- gml[gml$variable %in% keep_vars, , drop = FALSE]
    if (nrow(gml) == 0) {
      warning("No global metrics matched keep_vars; skipping global topology plot.")
    } else {
      p_global <- ggplot(gml, aes(x = variable, y = value, fill = Network)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        coord_flip() +
        labs(x = NULL, y = NULL, title = "Global topology comparison") +
        theme_minimal(base_size = 12)
      ggsave(file.path(outdir, "01_global_topology.png"), p_global, width = 10, height = 6, dpi = 200)
    }
  }


  # ---------- 2) ECDFs of node metrics ----------
  ecdf_long <- rbind(
    data.frame(Network=label1, degree=nm1$degree, strength=nm1$strength, btw=nm1$btw),
    data.frame(Network=label2, degree=nm2$degree, strength=nm2$strength, btw=nm2$btw)
  )
  plot_ecdf <- function(col, ttl) {
    ggplot(ecdf_long, aes_string(x=paste0("log1p(",col,")"), color="Network")) +
      stat_ecdf(geom="step", size=1) +
      labs(x=paste0("log1p(",ttl,")"), y="ECDF", title=paste0("ECDF of ", ttl)) +
      theme_minimal(base_size=12)
  }
  ggsave(file.path(outdir,"02_ecdf_degree.png"),   plot_ecdf("degree","degree"),   width=6, height=4, dpi=200)
  ggsave(file.path(outdir,"03_ecdf_strength.png"), plot_ecdf("strength","strength"), width=6, height=4, dpi=200)
  ggsave(file.path(outdir,"04_ecdf_betweenness.png"), plot_ecdf("btw","betweenness"), width=6, height=4, dpi=200)

  # ---------- 3) Node-to-node scatter (degree & strength) ----------
  scat <- function(xcol, ycol, ttl, diag=TRUE) {
    p <- ggplot(nodes_merged, aes_string(x=paste0("log1p(",xcol,")"),
                                         y=paste0("log1p(",ycol,")"))) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="gray60") +
      geom_point(alpha=0.5, size=1.2) +
      labs(x=paste0(label1," log1p(",ttl,")"), y=paste0(label2," log1p(",ttl,")"),
           title=paste0(ttl, ": ", label1, " vs ", label2)) +
      theme_minimal(base_size=12)
    p
  }
  p_deg <- scat("degree.1","degree.2","Degree")
  p_str <- scat("strength.1","strength.2","Strength")
  ggsave(file.path(outdir,"05_scatter_degree.png"), p_deg, width=6, height=5, dpi=200)
  ggsave(file.path(outdir,"06_scatter_strength.png"), p_str, width=6, height=5, dpi=200)

  # ---------- 4) Module correspondence heatmap (Jaccard) ----------
  jaccard_modules <- function(gA, gB) {
    clA <- cluster_louvain(gA, weights = abs(E(gA)$weight %||% 1))
    clB <- cluster_louvain(gB, weights = abs(E(gB)$weight %||% 1))
    setsA <- split(V(gA)$name, membership(clA))
    setsB <- split(V(gB)$name, membership(clB))
    A <- names(setsA); B <- names(setsB)
    J <- matrix(0, nrow=length(A), ncol=length(B), dimnames=list(paste0("A",A), paste0("B",B)))
    for (i in seq_along(setsA)) for (j in seq_along(setsB)) {
      inter <- length(intersect(setsA[[i]], setsB[[j]]))
      uni   <- length(union(setsA[[i]], setsB[[j]]))
      J[i,j] <- if (uni>0) inter/uni else 0
    }
    J
  }
  if (vcount(g1)>0 && vcount(g2)>0 && ecount(g1)>0 && ecount(g2)>0) {
    J <- jaccard_modules(g1, g2)
    pheatmap::pheatmap(J, cluster_rows = TRUE, cluster_cols = TRUE,
                       main = paste("Module overlap (Jaccard):", label1, "vs", label2),
                       filename = file.path(outdir, "07_module_overlap_jaccard.png"),
                       width = 8, height = 6)
  }

  # ---------- 5) Edge overlap / rewiring ----------
  edge_key <- function(g) {
    if (ecount(g)==0) return(character())
    el <- as_edgelist(g, names=TRUE)
    apply(el, 1, function(z) paste(sort(z), collapse="__"))
  }
  e1 <- edge_key(g1); e2 <- edge_key(g2)
  overlap <- length(intersect(e1,e2))
  only1   <- length(setdiff(e1,e2))
  only2   <- length(setdiff(e2,e1))
  edge_df <- data.frame(
    Class = c(paste0("Overlap (",label1,"∩",label2,")"),
              paste0("Only ",label1),
              paste0("Only ",label2)),
    Count = c(overlap, only1, only2)
  )
  p_edges <- ggplot(edge_df, aes(x=Class, y=Count, fill=Class)) +
    geom_col(width=0.6) + theme_minimal(base_size=12) +
    guides(fill="none") + labs(x=NULL, y="Edges", title="Edge overlap / rewiring")
  ggsave(file.path(outdir, "08_edge_overlap.png"), p_edges, width=7, height=4, dpi=200)

  invisible(list(global=gmtab, nodes_merged=nodes_merged))
}

analyze_one_graph <- function(g,
                              label,
                              outdir = "DCN_outputs",
                              orgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",
                              min_module_size = 20,
                              top_terms_per_module = 10,
                              top_n_nodes = 20) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  suppressPackageStartupMessages({
    library(igraph)
    library(openxlsx)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(circlize)  # <-- added
  })

  stopifnot(inherits(g, "igraph"))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  gdir      <- file.path(outdir, label)
  plots_dir <- file.path(gdir, "module_plots")
  dir.create(gdir,      showWarnings = FALSE, recursive = TRUE)
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  # 0) Clean → LCC
  g <- delete_vertices(g, degree(g) == 0)
  if (ecount(g) == 0 || vcount(g) < 2) {
    warning(sprintf("[%s] No edges after cleaning.", label))
    return(invisible(NULL))
  }
  cmp <- components(g)
  g   <- induced_subgraph(g, vids = which(cmp$membership == which.max(cmp$csize)))

  # 1) Communities (Louvain) on |weight|
  E(g)$w_abs <- abs(E(g)$weight %||% 0)
  cl <- cluster_louvain(g, weights = E(g)$w_abs)
  mb <- membership(cl)
  V(g)$module <- mb[V(g)$name]
  mod_sizes   <- sizes(cl)
  modules     <- split(names(mb), mb)
  ## =========================
  ## Module ranking (oncology-oriented)
  ## =========================
  # Per-module induced subgraphs + topological metrics
  module_ids <- sort(unique(mb))
  mod_rank_list <- lapply(module_ids, function(m) {
    vids <- V(g)[module == m]
    gm   <- induced_subgraph(g, vids = vids)

    n <- vcount(gm); e <- ecount(gm)
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

    # density
    dens <- edge_density(gm, loops = FALSE)

    # node strength (abs)
    st <- strength(gm, weights = abs(E(gm)$weight %||% 0))
    mean_st <- mean(st, na.rm = TRUE)

    # betweenness (weights as distance = 1/abs(w))
    w_len_m <- 1 / pmax(1e-9, abs(E(gm)$weight %||% 0))
    btw_m   <- betweenness(gm, directed = FALSE, weights = w_len_m, normalized = FALSE)
    max_btw <- max(btw_m, na.rm = TRUE)

    # hub concentration = top1 strength / sum strength
    hub_conc <- if (sum(st, na.rm = TRUE) > 0) max(st, na.rm = TRUE) / sum(st, na.rm = TRUE) else NA_real_

    # edge weight summaries
    wabs <- abs(E(gm)$weight %||% 0)
    mean_w <- mean(wabs, na.rm = TRUE)
    sd_w   <- stats::sd(wabs, na.rm = TRUE)

    # clustering
    clust <- tryCatch(transitivity(gm, type = "average"), error = function(e) NA_real_)

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

  # Filter: focus on interpretable modules (min size)
  module_rank_df$PassSize <- module_rank_df$N >= min_module_size

  # Oncology-oriented composite score:
  # - prefer dense + strong modules (coherent programs)
  # - also reward control points (betweenness)
  # - mild reward for hub-dominance (driver-like), but avoid extreme singletons via PassSize
  z <- function(x) { x <- as.numeric(x); if (all(is.na(x)) || stats::sd(x, na.rm=TRUE)==0) return(rep(0, length(x))); (x - mean(x, na.rm=TRUE))/stats::sd(x, na.rm=TRUE) }

  module_rank_df$Score <- with(module_rank_df,
    1.2*z(Density) +
    1.2*z(MeanAbsStrength) +
    0.8*z(MaxBetweenness) +
    0.4*z(HubConcentration) +
    0.6*z(Clustering)
  )

  module_rank_df <- module_rank_df %>%
    arrange(desc(PassSize), desc(Score), desc(N), desc(Density))

  # Annotate module "type" for figure logic
  module_rank_df$ModuleType <- with(module_rank_df, ifelse(
    !PassSize, "Too small / unstable",
    ifelse(HubConcentration >= 0.25 & N <= 60, "Hub-dominated driver-like core",
           ifelse(Density >= stats::quantile(Density[PassSize], 0.75, na.rm=TRUE), "Dense coherent program",
                  "Moderate program"))
  ))

  # pick top modules to highlight in figure
  top_modules <- head(module_rank_df$Module[module_rank_df$PassSize], 3)
  if (length(top_modules) == 0) top_modules <- head(module_rank_df$Module, 3)


  # 2) Topology summaries (for Excel)
  str_abs <- strength(g, weights = abs(E(g)$weight))
  deg     <- degree(g)
  w_len   <- 1 / pmax(1e-9, abs(E(g)$weight))
  btw     <- betweenness(g, directed = FALSE, weights = w_len, normalized = FALSE)

  ord_hub <- order(str_abs, decreasing = TRUE, na.last = TRUE)[seq_len(min(top_n_nodes, vcount(g)))]
  ord_btw <- order(btw,     decreasing = TRUE, na.last = TRUE)[seq_len(min(top_n_nodes, vcount(g)))]

  df_hubs <- data.frame(
    Gene        = V(g)$name[ord_hub],
    StrengthAbs = as.numeric(str_abs[ord_hub]),
    Degree      = as.integer(deg[ord_hub]),
    stringsAsFactors = FALSE
  )
  df_btw <- data.frame(
    Gene        = V(g)$name[ord_btw],
    Betweenness = as.numeric(btw[ord_btw]),
    StrengthAbs = as.numeric(str_abs[ord_btw]),
    Degree      = as.integer(deg[ord_btw]),
    stringsAsFactors = FALSE
  )

  # --- graph-level topology metrics for Excel ---
  n_nodes <- vcount(g); n_edges <- ecount(g)
  comps   <- components(g)
  n_comp  <- comps$no
  lcc_vids <- which(comps$membership == which.max(comps$csize))
  g_lcc    <- induced_subgraph(g, lcc_vids)
  lcc_nodes <- vcount(g_lcc); lcc_edges <- ecount(g_lcc)
  density_g      <- edge_density(g, loops = FALSE)
  avg_degree     <- mean(deg)
  transitivity_g <- transitivity(g, type = "global", isolates = "zero")
  assort_deg     <- assortativity_degree(g, directed = FALSE)
  v_conn <- tryCatch(vertex_connectivity(g), error = function(e) NA_real_)
  e_conn <- tryCatch(edge_connectivity(g),   error = function(e) NA_real_)
  apl_lcc <- tryCatch(average.path.length(g_lcc), error = function(e) NA_real_)
  modularity_g <- modularity(g, membership = mb, weights = E(g)$w_abs)

  topo_metrics <- data.frame(
    Label = label,
    Nodes = n_nodes,
    Edges = n_edges,
    Components = n_comp,
    LCC_Nodes = lcc_nodes,
    LCC_Edges = lcc_edges,
    Density = density_g,
    Avg_Degree = avg_degree,
    Transitivity = transitivity_g,
    Assortativity_Degree = assort_deg,
    Vertex_Connectivity = v_conn,
    Edge_Connectivity = e_conn,
    Avg_Path_Length_LCC = apl_lcc,
    Modularity_Louvain = modularity_g,
    stringsAsFactors = FALSE
  )

  # 3) Per-module plots (PNG) — thin edges, label top-5 hubs, size = (deg * 5)/n_max
  pal <- rainbow(length(unique(mb)))
  for (m in sort(unique(mb))) {
    vids <- V(g)[module == m]
    gm   <- induced_subgraph(g, vids = vids)
    if (vcount(gm) < 2 || ecount(gm) == 0) next
    if (vcount(gm) < min_module_size)       next

    wabs <- abs(E(gm)$weight)
    E(gm)$width <- 0.5
    E(gm)$color <- ifelse(E(gm)$weight >= 0, "forestgreen", "firebrick")
    lay         <-  layout_with_fr(gm, weights = wabs)


    degm  <- degree(gm)
    n_max <- max(degm, na.rm = TRUE)
    vsize <- if (n_max > 0) (degm * 5) / n_max else rep(2, length(degm))
    V(gm)$color <- pal[m]

    V(gm)$label <- NA
    top5 <- order(degm, decreasing = TRUE)[seq_len(min(5, length(degm)))]
    V(gm)$label[top5] <- V(gm)$name[top5]

    fpng <- file.path(plots_dir, sprintf("module_%s_n%d.png", m, vcount(gm)))
    png(fpng, width = 1600, height = 1300, res = 180)
    plot(gm, layout = lay,
         vertex.size = vsize,
         vertex.label = NULL,
         vertex.label.cex = 0.8,
         vertex.label.color = "black",
         main = sprintf("%s — Module %s (size: %d) — node size = degree; labels = top-5 hubs", label, m, vcount(gm)))
    dev.off()
  }

  # 4) Module-level circle plot with circlize (chords = cross edges; outer bars = module size)
  edf <- as.data.frame(as_edgelist(g), stringsAsFactors = FALSE)
  names(edf) <- c("v1","v2")
  edf$mod1 <- V(g)$module[ match(edf$v1, V(g)$name) ]
  edf$mod2 <- V(g)$module[ match(edf$v2, V(g)$name) ]
  cross <- subset(edf, mod1 != mod2)

  if (nrow(cross)) {
    # count cross-module edges (undirected: collapse A-B and B-A)
    key <- apply(cross[,c("mod1","mod2")], 1, function(z) paste(sort(z), collapse = "_"))
    w_counts <- as.data.frame(table(key), stringsAsFactors = FALSE)
    mods_split <- strsplit(w_counts$key, "_")
    links_df <- data.frame(
      from   = as.character(sapply(mods_split, `[`, 1)),
      to     = as.character(sapply(mods_split, `[`, 2)),
      weight = as.numeric(w_counts$Freq),
      stringsAsFactors = FALSE
    )

    # module sizes as named character vector (names must match sector labels)
    mod_sizes_vec <- setNames(as.integer(mod_sizes), names(mod_sizes))

    # sector order by size (largest first)
    sector_ord <- as.character(names(sort(mod_sizes_vec, decreasing = TRUE)))
    grid_col   <- setNames(grDevices::hcl.colors(length(sector_ord), "Dark 3"), sector_ord)

    # normalize bar heights to [0,1]
    h_max <- max(mod_sizes_vec, na.rm = TRUE)
    bar_h <- mod_sizes_vec / ifelse(h_max > 0, h_max, 1)

    # prepare links for chordDiagram
    links <- data.frame(
      from  = factor(links_df$from, levels = sector_ord),
      to    = factor(links_df$to,   levels = sector_ord),
      value = links_df$weight
    )

    outfile <- file.path(gdir, sprintf("%s_modules_circlize.png", label))
    grDevices::png(outfile, width = 1600, height = 1600, res = 200)
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 4, track.margin = c(0.01, 0.01))

    chordDiagram(
      x = links,
      order = sector_ord,
      grid.col = grid_col,
      transparency = 0.25,
      directional = 0,
      annotationTrack = c("grid"),
      preAllocateTracks = list(track.height = 0.12)
    )

    # track 1: module labels just outside the grid
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      si <- get.cell.meta.data("sector.index")
      circos.text(get.cell.meta.data("xcenter"),
                  get.cell.meta.data("ylim")[2] + 0.02,
                  labels = si, cex = 0.7, facing = "clockwise",
                  niceFacing = TRUE, adj = c(0, 0.5))
    }, bg.border = NA)

    # track 2: barplot of module sizes (height ∝ size)
    circos.trackPlotRegion(track.height = 0.12, panel.fun = function(x, y) {
      si   <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      xc   <- mean(xlim)
      w    <- diff(xlim) * 0.6
      h    <- bar_h[si]
      circos.rect(xleft = xc - w/2, ybottom = 0,
                  xright = xc + w/2, ytop = h,
                  col = 'white', border = 'black')
      circos.text(xc, h + 0.07, labels = mod_sizes_vec[si],
                  cex = 0.6, facing = "outside", adj = c(0.5, 0))
    }, bg.border = NA, ylim = c(0, 1))

    grid::grid.text(sprintf("%s — Module circle (bars = #genes)", label),
                    x = 0.5, y = 0.97)
    grDevices::dev.off()
  } else {
    message(sprintf("[%s] No cross-module edges; circlize plot skipped.", label))
  }

  # 5) Enrichment (GO:BP)
  enrichments <- lapply(names(modules), function(m) {
    genes <- modules[[m]]
    if (length(genes) < min_module_size) return(NULL)
    suppressWarnings(
      tryCatch({
        enrichGO(
          gene          = genes,
          OrgDb         = orgDb,
          keyType       = keyType,
          ont           = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff  = 0.05,
          qvalueCutoff  = 0.2,
          readable      = TRUE
        )
      }, error = function(e) NULL)
    )
  })
  names(enrichments) <- paste0("Module_", names(modules))

  make_enrichment_summary <- function(enrichments, top_n = NA) {
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
      if (!is.na(top_n)) df <- head(df, top_n)
      data.frame(
        Module    = module,
        Term      = df$Description %||% NA,
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
    rownames(out) <- NULL; out
  }
  summary_df <- make_enrichment_summary(enrichments, top_n = top_terms_per_module)

  # 6) Excel (adds topology metrics tab)
  openxlsx::addWorksheet(wb, "Module_Ranking")
  openxlsx::writeData(wb, "Module_Ranking", module_rank_df)

  # optional: conditional formatting for quick scan
  try({
    openxlsx::addStyle(wb, "Module_Ranking",
                       openxlsx::createStyle(textDecoration="bold"),
                       rows = 1, cols = 1:ncol(module_rank_df), gridExpand = TRUE)
    openxlsx::freezePane(wb, "Module_Ranking", firstRow = TRUE)
  }, silent = TRUE)

  wb <- createWorkbook()

  addWorksheet(wb, "Modules")
  writeData(wb, "Modules", data.frame(Gene = names(mb), Module = as.integer(mb)))

  addWorksheet(wb, "Top20_Hubs_Betweenness")
  writeData(wb, "Top20_Hubs_Betweenness", "Top 20 Hubs (by |strength|)", startRow = 1, startCol = 1)
  writeData(wb, "Top20_Hubs_Betweenness", df_hubs, startRow = 2, startCol = 1)
  start_row_btw <- nrow(df_hubs) + 4
  writeData(wb, "Top20_Hubs_Betweenness", "Top 20 Betweenness", startRow = start_row_btw, startCol = 1)
  writeData(wb, "Top20_Hubs_Betweenness", df_btw, startRow = start_row_btw + 1, startCol = 1)

  addWorksheet(wb, "Graph_Topology")
  writeData(wb, "Graph_Topology", topo_metrics)

  for (m in names(enrichments)) {
    enr <- enrichments[[m]]
    if (!is.null(enr)) {
      addWorksheet(wb, paste0("GO_", m))
      writeData(wb, paste0("GO_", m), as.data.frame(enr))
    }
  }

  addWorksheet(wb, "GO_BP_Summary")
  writeData(wb, "GO_BP_Summary", summary_df)
  addStyle(wb, "GO_BP_Summary", createStyle(textDecoration = "bold"),
           rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE)
  freezePane(wb, "GO_BP_Summary", firstRow = TRUE)

  xlsx_path <- file.path(gdir, sprintf("%s_results.xlsx", label))
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  saveRDS(g,      file = file.path(gdir, sprintf("%s_graph_cleaned.rds", label)))
  saveRDS(cl,     file = file.path(gdir, sprintf("%s_louvain_object.rds", label)))
  saveRDS(enrichments, file = file.path(gdir, sprintf("%s_enrichments.rds", label)))

  message(sprintf("[%s] Done. Module plots: %s | Excel: %s | Circlize: %s_modules_circlize.png",
                  label, plots_dir, xlsx_path, label))
  invisible(list(graph=g, membership=mb, sizes=mod_sizes,
                 hubs=df_hubs, betweenness=df_btw,
                 topo=topo_metrics,
                 enrichments=enrichments, enrichment_summary=summary_df,
                 excel=xlsx_path, plots_dir=plots_dir))
}


```
