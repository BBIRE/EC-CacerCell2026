#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(scales)
  library(tibble)
})

# ============================================================
# INPUTS
# ============================================================
group = 'EC'
if (group == 'EC') {
excel_file <- "network/network_06/corr06/DCN_outputs/EC/EC_results.xlsx"
rds_file   <- "network/network_06/corr06/DCN_graphs_6848ba9360c61375_EC_group1.rds"
out_pdf    <- "network/network_06/corr06/DCN_outputs/EC/module_metagraph.pdf"
out_png    <- "network/network_06/corr06/DCN_outputs/EC/module_metagraph.png"
}

if (group == 'Decidua') {
excel_file <- "network/network_06/corr06/DCN_outputs/Decidua/Decidua_results.xlsx"
rds_file   <- "network/network_06/corr06/DCN_graphs_937783c40de2b360_Decidua_group1.rds"
out_pdf    <- "network/network_06/corr06/DCN_outputs/Decidua/module_metagraph.pdf"
out_png    <- "network/network_06/corr06/DCN_outputs/Decidua/module_metagraph.png"
}

if (group == 'Healthy') {
excel_file <- "network/network_06/corr06/DCN_outputs/Healthy/Healthy_results.xlsx"
rds_file   <- "network/network_06/corr06/DCN_graphs_6848ba9360c61375_Healthy_group2.rds"
out_pdf    <- "network/network_06/corr06/DCN_outputs/Healthy/module_metagraph.pdf"
out_png    <- "network/network_06/corr06/DCN_outputs/Healthy/module_metagraph.png"
}


#!/usr/bin/env Rscript


sheet_gene_modules <- 1   # Sheet 1: Gene, Module
sheet_module_stats <- 2   # Sheet 2: Module summary with ModuleType

set.seed(123)

# ============================================================
# HELPERS
# ============================================================
clean_names_basic <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

normalize_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  x <- trimws(x)
  x
}

normalize_module_id <- function(x) {
  normalize_text(x)
}

normalize_moduletype <- function(x) {
  normalize_text(x)
}

flatten_chr <- function(x) {
  if (is.list(x)) {
    vapply(x, function(z) {
      if (is.null(z) || length(z) == 0 || all(is.na(z))) return(NA_character_)
      as.character(z[[1]])
    }, character(1))
  } else {
    as.character(x)
  }
}

flatten_num <- function(x) {
  if (is.list(x)) {
    suppressWarnings(
      vapply(x, function(z) {
        if (is.null(z) || length(z) == 0 || all(is.na(z))) return(NA_real_)
        as.numeric(z[[1]])
      }, numeric(1))
    )
  } else {
    suppressWarnings(as.numeric(x))
  }
}

rescale01_safe <- function(x, to = c(0, 1)) {
  if (length(x) == 0) return(numeric(0))
  if (all(is.na(x))) return(rep(mean(to), length(x)))
  xr <- range(x, na.rm = TRUE)
  if (xr[1] == xr[2]) return(rep(mean(to), length(x)))
  scales::rescale(x, to = to, from = xr)
}

find_first_col <- function(df, candidates, required = TRUE) {
  nm_clean <- clean_names_basic(names(df))
  cand_clean <- clean_names_basic(candidates)
  hit <- which(nm_clean %in% cand_clean)
  if (length(hit) == 0) {
    if (required) stop("Missing required column. Tried: ", paste(candidates, collapse = ", "))
    return(NULL)
  }
  names(df)[hit[1]]
}

# ============================================================
# GRAPH EXTRACTION
# ============================================================
extract_igraph_recursive <- function(x, max_depth = 6, depth = 0) {
  if (inherits(x, "igraph")) return(x)
  if (depth >= max_depth) return(NULL)

  if (is.list(x)) {
    for (nm in c("graph", "igraph", "g", "net", "network")) {
      if (!is.null(x[[nm]])) {
        out <- extract_igraph_recursive(x[[nm]], max_depth = max_depth, depth = depth + 1)
        if (!is.null(out)) return(out)
      }
    }
    for (i in seq_along(x)) {
      out <- extract_igraph_recursive(x[[i]], max_depth = max_depth, depth = depth + 1)
      if (!is.null(out)) return(out)
    }
  }
  NULL
}

extract_graph_from_rds <- function(obj) {
  g <- extract_igraph_recursive(obj)
  if (is.null(g)) stop("No igraph object found in the RDS file.")
  g
}

get_graph_gene_names <- function(g) {
  vn <- igraph::vertex_attr(g, "name")
  if (is.null(vn)) stop("Graph vertices do not have a 'name' attribute.")
  as.character(vn)
}

graph_to_edge_df <- function(g) {
  if (!inherits(g, "igraph")) stop("Input to graph_to_edge_df() is not an igraph object.")

  edf <- igraph::as_data_frame(g, what = "edges")
  if (nrow(edf) == 0) stop("Graph has no edges.")
  if (!all(c("from", "to") %in% names(edf))) {
    stop("Edge table does not contain required columns 'from' and 'to'.")
  }

  weight_col <- NULL
  for (cand in c("weight", "cor", "corr", "correlation", "value")) {
    if (cand %in% names(edf)) {
      weight_col <- cand
      break
    }
  }

  if (is.null(weight_col)) {
    edf$weight <- 1
  } else if (weight_col != "weight") {
    edf <- dplyr::rename(edf, weight = !!rlang::sym(weight_col))
  }

  edf %>%
    dplyr::mutate(
      from = as.character(from),
      to = as.character(to),
      weight = suppressWarnings(as.numeric(weight))
    )
}

# ============================================================
# NETWORK SUMMARIES
# ============================================================
compute_module_density_from_graph <- function(edge_df, gene_module_df) {
  gene_module_df <- gene_module_df %>%
    dplyr::mutate(
      gene = flatten_chr(gene),
      module = flatten_chr(module)
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(module), module != "") %>%
    dplyr::distinct(gene, module)

  edge_df <- edge_df %>%
    dplyr::mutate(
      from = as.character(from),
      to   = as.character(to)
    )

  mod_sizes <- gene_module_df %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(
      N_from_sheet1 = dplyr::n(),
      .groups = "drop"
    )

  within_edges <- edge_df %>%
    dplyr::inner_join(
      gene_module_df %>% dplyr::rename(from = gene, module_from = module),
      by = "from"
    ) %>%
    dplyr::inner_join(
      gene_module_df %>% dplyr::rename(to = gene, module_to = module),
      by = "to"
    ) %>%
    dplyr::filter(module_from == module_to) %>%
    dplyr::group_by(module = module_from) %>%
    dplyr::summarise(
      E_within = dplyr::n(),
      .groups = "drop"
    )

  mod_sizes %>%
    dplyr::left_join(within_edges, by = "module") %>%
    dplyr::mutate(
      E_within = dplyr::coalesce(E_within, 0L),
      E_possible = dplyr::if_else(
        N_from_sheet1 > 1,
        N_from_sheet1 * (N_from_sheet1 - 1) / 2,
        0
      ),
      density_from_graph = dplyr::if_else(
        E_possible > 0,
        E_within / E_possible,
        NA_real_
      )
    )
}

build_module_metagraph_edges <- function(edge_df, gene_module_df) {
  gene_module_df <- gene_module_df %>%
    dplyr::mutate(
      gene = flatten_chr(gene),
      module = flatten_chr(module)
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(module), module != "") %>%
    dplyr::distinct(gene, module)

  edge_df <- edge_df %>%
    dplyr::mutate(
      from = as.character(from),
      to   = as.character(to)
    )

  edge_df %>%
    dplyr::inner_join(
      gene_module_df %>% dplyr::rename(from = gene, module_from = module),
      by = "from"
    ) %>%
    dplyr::inner_join(
      gene_module_df %>% dplyr::rename(to = gene, module_to = module),
      by = "to"
    ) %>%
    dplyr::filter(module_from != module_to) %>%
    dplyr::mutate(
      m1 = pmin(module_from, module_to),
      m2 = pmax(module_from, module_to)
    ) %>%
    dplyr::group_by(m1, m2) %>%
    dplyr::summarise(
      E = dplyr::n(),
      median_weight = median(weight, na.rm = TRUE),
      mean_weight = mean(weight, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(from = m1, to = m2)
}

# ============================================================
# READ EXCEL
# ============================================================
sheet1_raw <- readxl::read_excel(excel_file, sheet = sheet_gene_modules)
sheet2_raw <- readxl::read_excel(excel_file, sheet = sheet_module_stats)

cat("\n===== Sheet1 original columns =====\n")
print(names(sheet1_raw))

cat("\n===== Sheet2 original columns =====\n")
print(names(sheet2_raw))

# Sheet 1
sheet1_df <- as.data.frame(sheet1_raw, stringsAsFactors = FALSE)
names(sheet1_df) <- clean_names_basic(names(sheet1_df))

gene_col <- find_first_col(
  sheet1_df,
  c("gene", "genes", "symbol", "gene_symbol", "external_gene_name")
)
module_col <- find_first_col(
  sheet1_df,
  c("module", "module_id", "module_name", "cluster", "color")
)

gene_module_df <- tibble::tibble(
  gene   = normalize_text(flatten_chr(sheet1_df[[gene_col]])),
  module = normalize_module_id(flatten_chr(sheet1_df[[module_col]]))
) %>%
  dplyr::filter(!is.na(gene), gene != "", !is.na(module), module != "") %>%
  dplyr::distinct()

# Sheet 2
required_sheet2_cols <- c("Module", "N", "E", "Density", "ModuleType")
missing_sheet2 <- setdiff(required_sheet2_cols, names(sheet2_raw))
if (length(missing_sheet2) > 0) {
  stop("Sheet 2 is missing required columns: ", paste(missing_sheet2, collapse = ", "))
}

sheet2_df <- as.data.frame(sheet2_raw, stringsAsFactors = FALSE)

module_stats_df <- tibble::tibble(
  module = normalize_module_id(flatten_chr(sheet2_df$Module)),
  N_sheet2 = flatten_num(sheet2_df$N),
  E_sheet2 = flatten_num(sheet2_df$E),
  Density_sheet2 = flatten_num(sheet2_df$Density),
  ModuleType = normalize_moduletype(flatten_chr(sheet2_df$ModuleType))
) %>%
  dplyr::filter(!is.na(module), module != "", !is.na(ModuleType), ModuleType != "") %>%
  dplyr::distinct()

cat("\n===== Structure of imported Sheet1 =====\n")
print(str(gene_module_df))

cat("\n===== Structure of imported Sheet2 =====\n")
print(str(module_stats_df))

cat("\n===== Original ModuleType values in Sheet2 =====\n")
print(sort(unique(module_stats_df$ModuleType)))

cat("\n===== Counts per ModuleType in Sheet2 =====\n")
print(module_stats_df %>% dplyr::group_by(ModuleType) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop"))

module_stats_df <- module_stats_df %>%
  dplyr::filter(ModuleType != "Too small / unstable")

cat("\n===== Counts per ModuleType after exclusion =====\n")
print(module_stats_df %>% dplyr::group_by(ModuleType) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop"))

# Keep only modules present in sheet 2
gene_module_df <- gene_module_df %>%
  dplyr::semi_join(module_stats_df %>% dplyr::select(module), by = "module")

if (nrow(gene_module_df) == 0) {
  stop("No gene-module assignments remain after joining Sheet1 to Sheet2.")
}

# ============================================================
# READ GRAPH
# ============================================================
obj <- readRDS(rds_file)
g_gene <- extract_graph_from_rds(obj)

if (!inherits(g_gene, "igraph")) {
  stop("Extracted object is not an igraph.")
}

graph_genes <- get_graph_gene_names(g_gene)

gene_module_df <- gene_module_df %>%
  dplyr::filter(gene %in% graph_genes)

if (nrow(gene_module_df) == 0) {
  stop("No genes from Sheet1 match graph vertex names.")
}

edge_df <- graph_to_edge_df(g_gene) %>%
  dplyr::filter(from %in% gene_module_df$gene, to %in% gene_module_df$gene)

if (nrow(edge_df) == 0) {
  stop("No graph edges remain after filtering to mapped genes.")
}

# ============================================================
# NODE METADATA
# ============================================================
density_df <- compute_module_density_from_graph(edge_df, gene_module_df)

node_df <- gene_module_df %>%
  dplyr::distinct(gene, module) %>%
  dplyr::group_by(module) %>%
  dplyr::summarise(N_from_sheet1 = dplyr::n(), .groups = "drop") %>%
  dplyr::left_join(module_stats_df, by = "module") %>%
  dplyr::left_join(density_df %>% dplyr::select(module, density_from_graph), by = "module") %>%
  dplyr::mutate(
    size_plot_source = N_sheet2,
    alpha_plot_source = Density_sheet2
  ) %>%
  dplyr::filter(!is.na(ModuleType))

cat("\n===== Final modules per ModuleType =====\n")
print(node_df %>% dplyr::group_by(ModuleType) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop"))

# ============================================================
# MODULE-MODULE EDGES
# ============================================================
meta_edges <- build_module_metagraph_edges(edge_df, gene_module_df) %>%
  dplyr::filter(from %in% node_df$module, to %in% node_df$module)

if (nrow(meta_edges) == 0) {
  stop("No inter-module edges found.")
}

# ============================================================
# BUILD GRAPH
# ============================================================
g_meta <- igraph::graph_from_data_frame(
  d = meta_edges,
  vertices = node_df %>% dplyr::rename(name = module),
  directed = FALSE
)

# remove isolated nodes
g_meta <- igraph::delete_vertices(
  g_meta,
  igraph::V(g_meta)[igraph::degree(g_meta) == 0]
)

if (igraph::vcount(g_meta) == 0) {
  stop("All modules are isolated after metagraph construction.")
}

igraph::V(g_meta)$ModuleType <- normalize_moduletype(igraph::V(g_meta)$ModuleType)

# ============================================================
# PLOT ATTRIBUTES
# ============================================================
igraph::V(g_meta)$size_plot  <- rescale01_safe(igraph::V(g_meta)$size_plot_source, to = c(5, 18))
igraph::V(g_meta)$alpha_plot <- rescale01_safe(igraph::V(g_meta)$alpha_plot_source, to = c(0.35, 1.0))
igraph::E(g_meta)$width_plot <- rescale01_safe(igraph::E(g_meta)$E, to = c(0.5, 4))

pal <- c(
  "Dense coherent program" = "#D55E00",
  "Moderate program"       = "#009E73",
  "Sparse, large"          = "#0072B2"
)

cat("\n===== ModuleType values in final plotted graph =====\n")
print(sort(unique(igraph::V(g_meta)$ModuleType)))

cat("\n===== Palette-mismatched ModuleType values =====\n")
print(setdiff(sort(unique(igraph::V(g_meta)$ModuleType)), names(pal)))

# ============================================================
# PLOT
# ============================================================
p <- ggraph::ggraph(g_meta, layout = "star") +
  ggraph::geom_edge_link(
    aes(width = width_plot),
    colour = "grey65",
    alpha = 0.4,
    show.legend = FALSE
  ) +
  ggraph::geom_node_point(
    aes(size = size_plot, colour = ModuleType, alpha = 1)
  ) +
  ggraph::geom_node_text(
    aes(label = name),
    repel = F,
    size = 5,
    family = "Helvetica"
  ) +
  ggplot2::scale_colour_manual(
    values = pal,
    breaks = c("Dense coherent program", "Moderate program", "Sparse, large"),
    drop = FALSE,
    na.value = "magenta",
    name = "ModuleType"
  ) +
  ggplot2::scale_size_identity() +
  ggplot2::scale_alpha_identity() +
  ggraph::theme_graph(base_family = "Helvetica") +
  ggplot2::ggtitle("Module metagraph") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
  )

# ============================================================
# SAVE
# ============================================================
ggplot2::ggsave(out_png, p, width = 9, height = 7, dpi = 300)

try({
  ggplot2::ggsave(
    out_pdf,
    p,
    width = 9,
    height = 7,
    device = cairo_pdf
  )
}, silent = TRUE)

# ============================================================
# EXPORT TABLES
# ============================================================
node_export <- tibble::tibble(
  module = igraph::V(g_meta)$name,
  ModuleType = igraph::V(g_meta)$ModuleType,
  N_sheet2 = igraph::V(g_meta)$N_sheet2,
  E_sheet2 = igraph::V(g_meta)$E_sheet2,
  Density_sheet2 = igraph::V(g_meta)$Density_sheet2,
  N_from_sheet1 = igraph::V(g_meta)$N_from_sheet1,
  density_from_graph = igraph::V(g_meta)$density_from_graph,
  degree = igraph::degree(g_meta)
)

edge_export <- igraph::as_data_frame(g_meta, what = "edges") %>%
  tibble::as_tibble()

write.csv(node_export, sub("\\.png$", "_nodes.csv", out_png), row.names = FALSE)
write.csv(edge_export, sub("\\.png$", "_edges.csv", out_png), row.names = FALSE)

cat("\nDone.\n")
cat("Saved: ", out_png, "\n", sep = "")
cat("Saved: ", out_pdf, " (if Cairo PDF worked)\n", sep = "")
cat("Saved: ", sub("\\.png$", "_nodes.csv", out_png), "\n", sep = "")
cat("Saved: ", sub("\\.png$", "_edges.csv", out_png), "\n", sep = "")
