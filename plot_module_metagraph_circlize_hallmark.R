# =========================================================
# INPUTS EXPECTED
# =========================================================

plot_title <- paste(group, "Module Metagraph", sep = " ") 
min_edge_weight_to_plot <- 1
label_cex <- 0.9
grid_border_col <- "white"
background_col <- "white"
hallmark_sheet <- "HALLMARK_Summary"
hallmark_q_cutoff <- 0.0001

# =========================================================
# HELPERS
# =========================================================
clean_names_basic <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

normalize_moduletype_bio <- function(x) {
  x <- normalize_text(x)

  dplyr::case_when(
    x == "Dense + rewired (highest priority)" ~ "Core rewired program",
    x == "Moderate program" ~ "Intermediate program",
    x == "Sparse, large (usually noise)" ~ "Diffuse / low-coherence program",
    x == "Too small / unstable" ~ "Small / unstable module",
    TRUE ~ x
  )
}

moduletype_palette_bio <- c(
  "Core rewired program" = "#B2182B",
  "Intermediate program" = "#F4A582",
  "Diffuse / low-coherence program" = "#2166AC",
  "Small / unstable module" = "gray"
)

normalize_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  x <- trimws(x)
  x
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

format_module_label <- function(x) {
  x <- normalize_text(x)
  paste0("module ", x)
}

# ---------------------------------------------------------
# Hallmark functional grouping
# ---------------------------------------------------------
hallmark_to_group <- function(pathway) {
  p <- toupper(normalize_text(pathway))

  dplyr::case_when(
    p %in% c(
      "HALLMARK_E2F_TARGETS",
      "HALLMARK_G2M_CHECKPOINT",
      "HALLMARK_MYC_TARGETS_V1",
      "HALLMARK_MYC_TARGETS_V2",
      "HALLMARK_MITOTIC_SPINDLE",
      "HALLMARK_DNA_REPAIR"
    ) ~ "Proliferation",

    p %in% c(
      "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
      "HALLMARK_GLYCOLYSIS",
      "HALLMARK_FATTY_ACID_METABOLISM",
      "HALLMARK_ADIPOGENESIS",
      "HALLMARK_MTORC1_SIGNALING",
      "HALLMARK_PROTEIN_SECRETION",
      "HALLMARK_ESTROGEN_RESPONSE_EARLY",
      "HALLMARK_ESTROGEN_RESPONSE_LATE"
    ) ~ "Metabolism",

    p %in% c(
      "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_IL6_JAK_STAT3_SIGNALING",
      "HALLMARK_IL2_STAT5_SIGNALING",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_COMPLEMENT",
      "HALLMARK_ALLOGRAFT_REJECTION",
      "HALLMARK_KRAS_SIGNALING_UP"
    ) ~ "Immune",

    p %in% c(
      "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
      "HALLMARK_HYPOXIA",
      "HALLMARK_UV_RESPONSE_DN",
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
      "HALLMARK_MYOGENESIS"
    ) ~ "Stress",

    TRUE ~ "Other"
  )
}

hallmark_group_colors <- c(
  "Proliferation" = "#D55E00",
  "Metabolism"    = "#0072B2",
  "Immune"        = "#009E73",
  "Stress"        = "#CC79A7",
  "Other"         = "#999999"
)

# =========================================================
# PREPARE NODES
# =========================================================
node_df2 <- node_df %>%
  dplyr::transmute(
    module = as.character(module),
    ModuleType = normalize_moduletype_bio(as.character(ModuleType)),
    N_sheet2 = as.numeric(N_sheet2)
  ) %>%
  dplyr::mutate(
    module = normalize_text(module),
    ModuleType = normalize_text(ModuleType)
  ) %>%
  dplyr::filter(!is.na(module), module != "", !is.na(ModuleType), ModuleType != "", !is.na(N_sheet2)) %>%
  dplyr::distinct()

cat("\n===== ModuleType labels found in node_df =====\n")
print(sort(unique(node_df2$ModuleType)))

node_df2 <- node_df2 %>%
  dplyr::arrange(ModuleType, dplyr::desc(N_sheet2))

modules <- node_df2$module

# =========================================================
# BUILD PALETTE FROM FOUND LABELS ONLY
# =========================================================
found_types <- unique(node_df2$ModuleType)

pal_use <- moduletype_palette_bio[names(moduletype_palette_bio) %in% found_types]

missing_types <- setdiff(found_types, names(pal_use))
if (length(missing_types) > 0) {
  extra_cols <- setNames(scales::hue_pal()(length(missing_types)), missing_types)
  pal_use <- c(pal_use, extra_cols)
}

pal_use <- pal_use[unique(node_df2$ModuleType)]

cat("\n===== ModuleType color mapping used =====\n")
print(pal_use)

module_colors <- pal_use[node_df2$ModuleType]
names(module_colors) <- node_df2$module

# =========================================================
# PREPARE EDGES
# =========================================================
edge_df2 <- meta_edges %>%
  dplyr::transmute(
    from = as.character(from),
    to   = as.character(to),
    E    = as.numeric(E)
  ) %>%
  dplyr::mutate(
    from = normalize_text(from),
    to   = normalize_text(to)
  ) %>%
  dplyr::filter(
    !is.na(from), !is.na(to), !is.na(E),
    from %in% modules, to %in% modules,
    from != to,
    E >= min_edge_weight_to_plot
  )

if (nrow(edge_df2) > 0) {
  edge_df2$link_lwd <- scales::rescale(edge_df2$E, to = c(1, 5))
} else {
  edge_df2$link_lwd <- numeric(0)
}

edge_df2$link_col <- grDevices::adjustcolor("grey", alpha.f = 0.45)

# =========================================================
# READ HALLMARK SUMMARY
# =========================================================
hallmark_raw <- readxl::read_excel(excel_file, sheet = hallmark_sheet)
hallmark_df <- as.data.frame(hallmark_raw, stringsAsFactors = FALSE)

cat("\n===== HALLMARK_summary original columns =====\n")
print(names(hallmark_df))

module_col_hm <- find_first_col(
  hallmark_df,
  c("module", "module_id", "cluster", "module_name")
)

pathway_col_hm <- find_first_col(
  hallmark_df,
  c("pathway", "term", "hallmark", "geneset", "gene_set", "description")
)

qval_col_hm <- find_first_col(
  hallmark_df,
  c("qval", "q_value", "q.value", "padj", "fdr", "adj_p_val", "adj_pvalue")
)

count_col_hm <- find_first_col(
  hallmark_df,
  c("count", "gene_count", "overlap", "genes_in_term")
)

hallmark_df2 <- tibble::tibble(
  module = normalize_text(as.character(hallmark_df[[module_col_hm]])),
  pathway = normalize_text(as.character(hallmark_df[[pathway_col_hm]])),
  qval = suppressWarnings(as.numeric(hallmark_df[[qval_col_hm]])),
  Count = suppressWarnings(as.numeric(hallmark_df[[count_col_hm]]))
) %>%
  dplyr::filter(
    !is.na(module), module != "",
    !is.na(pathway), pathway != "",
    !is.na(qval), qval < hallmark_q_cutoff,
    !is.na(Count), Count > 0,
    module %in% modules
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    HallmarkGroup = hallmark_to_group(pathway)
  )

cat("\n===== Hallmark pathways retained (qval < ", hallmark_q_cutoff, ") =====\n", sep = "")
print(hallmark_df2 %>% dplyr::group_by(module) %>% dplyr::summarise(n_pathways = dplyr::n(), .groups = "drop"))

cat("\n===== Hallmark functional groups found =====\n")
print(hallmark_df2 %>% dplyr::count(HallmarkGroup, sort = TRUE))

# consistent pathway ordering across modules:
# order first by functional group, then by frequency, then pathway name
pathway_order_df <- hallmark_df2 %>%
  dplyr::group_by(HallmarkGroup, pathway) %>%
  dplyr::summarise(freq = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    HallmarkGroup = factor(HallmarkGroup, levels = c("Proliferation", "Metabolism", "Immune", "Stress", "Other"))
  ) %>%
  dplyr::arrange(HallmarkGroup, dplyr::desc(freq), pathway)

pathway_levels <- pathway_order_df$pathway
hallmark_df2$pathway <- factor(hallmark_df2$pathway, levels = pathway_levels)

# color by functional group
hallmark_df2$pathway_color <- hallmark_group_colors[hallmark_df2$HallmarkGroup]

hallmark_split <- split(hallmark_df2, hallmark_df2$module)

# pretty y-axis breaks for enrichment barplot track
global_max_count <- max(hallmark_df2$Count, na.rm = TRUE)
y_breaks <- pretty(c(0, global_max_count), n = 4)
y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= global_max_count]

# =========================================================
# ORDERING / WIDTHS
# =========================================================
sector_width <- scales::rescale(node_df2$N_sheet2, to = c(6, 18))
names(sector_width) <- node_df2$module

# =========================================================
# EXPORT SETTINGS
# =========================================================
out_png <- paste("network/network_06/corr06/DCN_outputs/", group, "/module_metagraph_circlize.png", sep = "")

png(
  filename = out_png,
  width = 4000,
  height = 3000,
  res = 300
)

# =========================================================
# LAYOUT: left = circos, right = legend
# =========================================================
layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1.5))

# =========================================================
# LEFT PANEL → CIRCOS
# =========================================================
par(mar = c(1, 1, 2, 1))

circos.clear()

# make sure modules are unique and in plotting order
node_df2 <- node_df2 %>%
  dplyr::distinct(module, .keep_all = TRUE)

modules <- node_df2$module

# build one gap value per sector
class_id <- as.character(node_df2$ModuleType)
gap_after <- rep(2, length(modules))
block_end <- c(class_id[-1] != class_id[-length(class_id)], TRUE)
gap_after[block_end] <- 8
stopifnot(length(gap_after) == length(modules))

circos.par(
  start.degree = 90,
  gap.after = gap_after,
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0.002, 0.002),
  points.overflow.warning = FALSE,
  canvas.xlim = c(-1.2, 1.2),
  canvas.ylim = c(-1.2, 1.2)
)

circos.initialize(
  factors = modules,
  xlim = cbind(rep(0, length(modules)), sector_width[modules])
)

# ---------------------------------------------------------
# TRACK 1: module class color ring
# ---------------------------------------------------------
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.055,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")

    circos.rect(
      xleft = xlim[1], ybottom = 0,
      xright = xlim[2], ytop = 1,
      col = module_colors[sector.name],
      border = grid_border_col,
      lwd = 0.8
    )
  }
)

# ---------------------------------------------------------
# TRACK 2: module labels, just inside the class ring
# ---------------------------------------------------------
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.07,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xcenter <- mean(get.cell.meta.data("xlim"))
    lab <- format_module_label(sector.name)

    circos.text(
      x = xcenter,
      y = 0.5,
      labels = lab,
      facing = "bending.inside",
      niceFacing = TRUE,
      adj = c(0.5, 0.5),
      cex = label_cex,
      col = "black",
      font = 2
    )
  }
)

# ---------------------------------------------------------
# TRACK 3: hallmark circular barplots with y-scale
# ---------------------------------------------------------
circos.trackPlotRegion(
  ylim = c(0, global_max_count),
  track.height = 0.22,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    hm <- hallmark_split[[sector.name]]

    circos.rect(
      xleft = xlim[1],
      ybottom = 0,
      xright = xlim[2],
      ytop = global_max_count,
      col = "#F7F7F7",
      border = NA
    )

    for (yb in y_breaks) {
      circos.segments(
        x0 = xlim[1], y0 = yb,
        x1 = xlim[2], y1 = yb,
        col = "#D9D9D9",
        lwd = 0.5
      )
    }

    if (!is.null(hm) && nrow(hm) > 0) {
      hm <- hm %>%
        dplyr::arrange(HallmarkGroup, pathway)

      n_bars <- nrow(hm)
      x_lefts <- seq(xlim[1], xlim[2], length.out = n_bars + 1)[-(n_bars + 1)]
      x_rights <- seq(xlim[1], xlim[2], length.out = n_bars + 1)[-1]
      x_mids <- (x_lefts + x_rights) / 2
      bar_width <- (x_rights - x_lefts) * 0.88

      circos.barplot(
        value = hm$Count,
        pos = x_mids,
        bar_width = bar_width,
        col = hm$pathway_color,
        border = NA
      )
    }

    if (sector.name == modules[1]) {
      circos.yaxis(
        side = "left",
        at = y_breaks,
        labels = y_breaks,
        labels.cex = 0.45,
        tick.length = 0.01,
        lwd = 0.6
      )
    }
  }
)

# ---------------------------------------------------------
# LINKS
# ---------------------------------------------------------
for (i in seq_len(nrow(edge_df2))) {
  x1 <- mean(get.cell.meta.data("xlim", sector.index = edge_df2$from[i], track.index = 1))
  x2 <- mean(get.cell.meta.data("xlim", sector.index = edge_df2$to[i], track.index = 1))

  circos.link(
    sector.index1 = edge_df2$from[i],
    point1 = x1,
    sector.index2 = edge_df2$to[i],
    point2 = x2,
    col = edge_df2$link_col[i],
    lwd = edge_df2$link_lwd[i],
    border = NA
  )
}

title(plot_title, cex.main = 2)

# =========================================================
# RIGHT PANEL → LEGEND
# =========================================================
par(mar = c(2, 2, 2, 2))
plot.new()

legend(
  "top",
  legend = names(pal_use),
  fill = unname(pal_use),
  border = NA,
  bty = "n",
  cex = 1,
  title = "Module class"
)

legend(
  "center",
  legend = c("Proliferation", "Metabolism", "Immune", "Stress", "Other"),
  fill = unname(hallmark_group_colors[c("Proliferation", "Metabolism", "Immune", "Stress", "Other")]),
  border = NA,
  bty = "n",
  cex = 0.9,
  title = "Hallmark functional groups"
)

# =========================================================
# CLOSE DEVICE
# =========================================================
dev.off()

cat("Saved PNG: ", out_png, "\n")
