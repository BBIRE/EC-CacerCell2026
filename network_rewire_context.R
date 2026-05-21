suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# ----------------------------
# Inputs (edit paths)
# ----------------------------
gene_file   <- "top_rewired_genes_modules_k_neighbors.csv"   
module_xlsx <- "results.xlsx" 
module_sheet <- "Module_Ranking"

# ----------------------------
# Load
# ----------------------------
genes <- read_csv(gene_file, show_col_types = FALSE) %>%
  mutate(
    Module_case = as.character(Module_case),
    Module_ctrl = as.character(Module_ctrl)
  )

# ----------------------------
# Load module ranking from Excel (sheet: Module_Ranking)
# ----------------------------
mods_raw <- readxl::read_excel(module_xlsx, sheet = module_sheet) %>%
  as.data.frame()

# Helper: robust numeric parsing (handles 0,123 and 0.123)
to_num <- function(x) readr::parse_number(as.character(x), locale = locale(decimal_mark = ","))

# Standardize expected columns (adjust if your headers differ)
# Requires at least: Module, N, Score, ModuleType, PassSize
mods <- mods_raw %>%
  mutate(
    Module = as.character(.data$Module),
    N      = to_num(.data$N),
    Score  = to_num(.data$Score),
    PassSize = as.logical(.data$PassSize),
    ModuleType = as.character(.data$ModuleType)
  )

# ----------------------------
# Derive rank (1 = highest score)
# ----------------------------
mods_ranked <- mods %>%
  arrange(desc(Score)) %>%
  mutate(
    Rank = row_number(),
    Tier = case_when(
      str_detect(ModuleType, regex("^Dense\\s*\\+\\s*rewired", ignore_case = TRUE)) ~ "Dense+rewired",
      str_detect(ModuleType, regex("Moderate program", ignore_case = TRUE)) ~ "Moderate",
      str_detect(ModuleType, regex("^Sparse", ignore_case = TRUE)) ~ "Sparse/large",
      TRUE ~ "Other"
    )
  )

# ----------------------------
# Join module attributes for case & ctrl assignments
# ----------------------------
genes_ctx <- genes %>%
  left_join(
    mods_ranked %>%
      transmute(Module_case = Module,
                Score_case = Score,
                Rank_case = Rank,
                ModuleType_case = ModuleType,
                Tier_case = Tier,
                PassSize_case = PassSize,
                N_case = N),
    by = "Module_case"
  ) %>%
  left_join(
    mods_ranked %>%
      transmute(Module_ctrl = Module,
                Score_ctrl = Score,
                Rank_ctrl = Rank,
                ModuleType_ctrl = ModuleType,
                Tier_ctrl = Tier,
                PassSize_ctrl = PassSize,
                N_ctrl = N),
    by = "Module_ctrl"
  ) %>%
  mutate(
    ModuleSwitch = Module_case != Module_ctrl,

    # Positive ΔScore means gene moved to a higher-scoring module in CASE relative to CTRL assignment
    DeltaScore = Score_case - Score_ctrl,

    # Rank: 1 is best. Negative ΔRank means improvement (moved to better rank) in CASE
    DeltaRank = Rank_case - Rank_ctrl,

    ShiftDirection = case_when(
      !ModuleSwitch ~ "No switch",
      is.na(DeltaScore) ~ "Unknown (missing module info)",
      DeltaScore > 0 ~ "To higher-scoring module (case)",
      DeltaScore < 0 ~ "To lower-scoring module (case)",
      TRUE ~ "No score change"
    ),

    TierTransition = if_else(
      is.na(Tier_ctrl) | is.na(Tier_case),
      NA_character_,
      paste0(Tier_ctrl, " → ", Tier_case)
    ),

    ModuleTypeTransition = if_else(
      is.na(ModuleType_ctrl) | is.na(ModuleType_case),
      NA_character_,
      paste0(ModuleType_ctrl, " → ", ModuleType_case)
    ),

    # A simple “contextual priority” label you can put in figures/tables
    ContextLabel = case_when(
      ModuleSwitch & Tier_ctrl == "Sparse/large" & Tier_case != "Sparse/large" ~ "Escapes sparse/noise module",
      ModuleSwitch & Tier_ctrl != "Sparse/large" & Tier_case == "Sparse/large" ~ "Falls into sparse/noise module",
      ModuleSwitch & Tier_ctrl != Tier_case & Tier_case == "Dense+rewired" ~ "Moves into dense+rewired program",
      ModuleSwitch & Tier_ctrl != Tier_case & Tier_ctrl == "Dense+rewired" ~ "Moves out of dense+rewired program",
      ModuleSwitch & DeltaScore > 0 ~ "Escalates module priority",
      ModuleSwitch & DeltaScore < 0 ~ "De-escalates module priority",
      ModuleSwitch ~ "Switch (similar priority)",
      TRUE ~ "No switch"
    )
  )

# ----------------------------
# Summaries
# ----------------------------

# 1) How many genes move between tiers / module types
tier_flow <- genes_ctx %>%
  filter(ModuleSwitch) %>%
  count(TierTransition, sort = TRUE)

type_flow <- genes_ctx %>%
  filter(ModuleSwitch) %>%
  count(ModuleTypeTransition, sort = TRUE)

# 2) Which modules attract the highest rewiring genes (case modules)
top_case_modules_by_rewire <- genes_ctx %>%
  group_by(Module_case, ModuleType_case, Tier_case) %>%
  summarise(
    n_genes = n(),
    median_rewire = median(RewireScore, na.rm = TRUE),
    median_degdiff = median(DegreeDiff, na.rm = TRUE),
    median_log2k = median(log2_k_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_rewire), desc(n_genes))

# 3) Top “contextualized shifts”: big rewiring + big module priority change
top_gene_shifts <- genes_ctx %>%
  filter(ModuleSwitch) %>%
  arrange(desc(RewireScore), desc(abs(DeltaScore))) %>%
  select(Gene, Module_ctrl, Module_case,
         ModuleType_ctrl, ModuleType_case,
         Tier_ctrl, Tier_case,
         Score_ctrl, Score_case, DeltaScore,
         Rank_ctrl, Rank_case, DeltaRank,
         k_ctrl, k_case, log2_k_ratio,
         RewireScore, DegreeDiff,
         ContextLabel)

# ----------------------------
# Write outputs
# ----------------------------
write_tsv(genes_ctx, "gene_rewire_with_module_context.tsv")
write_tsv(tier_flow, "tier_transition_counts.tsv")
write_tsv(type_flow, "moduletype_transition_counts.tsv")
write_tsv(top_case_modules_by_rewire, "case_modules_ranked_by_rewiring.tsv")
write_tsv(top_gene_shifts, "top_gene_shifts_contextualized.tsv")

message("Done.
- gene_rewire_with_module_context.tsv
- tier_transition_counts.tsv
- moduletype_transition_counts.tsv
- case_modules_ranked_by_rewiring.tsv
- top_gene_shifts_contextualized.tsv")
