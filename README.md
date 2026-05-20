# EC-CacerCell2026
Code repository for Ferretti et al. 2026 Cancer Cell
EC-CacerCell2026
├── network_rewire_context.R
├── network_analysis_core.R
└── plot_module_metagraph_circlize_hallmark.R

## network_analysis_core.R - Differential Co-Expression Network Analysis Pipeline

This pipeline implements a complete workflow to reconstruct, analyze, and compare gene co-expression networks across three biological conditions:
- Healthy endometrium
- Decidua
- Endometrial cancer (EC)
The objective is to quantify how transcriptional coordination is reorganized during physiological invasion (decidualization) and pathological invasion (cancer), identify rewired genes and modules, and characterize their biological functions.

### 1. Input Data
  1.1 Expression Matrices: Gene-level TPM matrix used for correlation-based network construction. Raw count matrix used for differential expression analyses.
  1.2 Sample Metadata: Clinical and biological annotation table.

### 2. Required Packages
   dplyr, tibble, stringr, purrr, tidyr, Data Import, readxl, readr, Network Analysis, igraph, networktools, Functional Enrichment, clusterProfiler, msigdbr, Gene Ontology Consortium, Visualization, ggplot2, ComplexHeatmap, Export, openxlsx, CRAN, Bioconductor

### 3. Analysis Steps and Functions
  #### Step 1 — Data Import and Harmonization: Load expression matrices and metadata, standardize sample identifiers, and ensure consistent sample matching.
  #### Step 2 — Differential Co-Expression Graph Construction: Build weighted gene co-expression networks and identify edges significantly altered between conditions. For each pairwise comparison:
   - Computes gene–gene correlations.
   - Applies significance and effect-size thresholds.
   - Builds condition-specific networks.
   - Constructs differential rewiring networks.
   - Key Parameters: alpha = 0.01, t_percentile = 0.90, diff_t = 0.3, fdr_method = "BH"

  #### Step 3 — Single Network Topological Analysis: Characterize the internal structure of each graph.
  - Graph Cleaning
  - Removes isolated nodes
  - Extracts largest connected component
  - Community Detection
  - Louvain clustering using absolute edge weights
  - Node-Level Metrics 
  - Global Metrics 

#### Step 4 — Module Prioritization: Identify biologically important co-expression modules.
- Calculates module size, density, clustering, and connectivity.
- Ranks modules using a composite score.
- Assigns labels such as: Dense + rewired, Large dense, Hub-dominated, Moderate program, Sparse large

#### Step 5 — Functional Enrichment: Assign biological meaning to modules.
- Gene Ontology Biological Process enrichment (MSigDB Hallmark enrichment)

#### Step 6 — Network Visualization: Generate publication-quality figures.
- Constructs a network where nodes represent modules.
- Visualizes inter-module connectivity.
- Visualizes module lables and enrichment

#### Step 7 — Differential Rewiring Analysis
- Quantify structural changes between networks: Edge overlap (Jaccard), Edge-rank correlation, Hub preservation, Gain/loss balance, Fraction of positive changes

#### Step 8 — Cross-Network Comparison
Global topology comparison: Degree/strength distributions, Shared-node scatterplots, Module overlap heatmaps, Edge overlap summaries

### 4. Macroscopic Biological Interpretation
The pipeline addresses a central biological question: How does the transcriptional coordination of healthy endometrium reorganize under physiological and pathological invasive conditions?
- Differential graphs identify genes whose interaction patterns change most strongly between conditions, highlighting candidate regulatory switches.
- Co-Expression Modules
- Modules represent coordinated transcriptional programs
- Rewiring Metrics (Preservation of existing hubs, Global decoherence, Emergence of new regulatory programs, Redistribution of connectivity)
- Network Comparisons: Quantify the extent to which decidua and cancer share or diverge in their regulatory architectures.
