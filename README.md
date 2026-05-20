# Ferretti, Betti et. al 2026 
# Code repository

---
# `TIME_analysis.R`

## 1. Input Data

### 1.1 CIBERSORTx Fraction Matrix

* CibersortX tsv files for run in relativ and absolute mode

### 1.2 Sample Metadata

* Clinical and biological annotation table used to define comparison groups.

### 1.3 IHC quantifications
* Infiltration scores from for digital pathology IHC staining

## 2. Required Packages

`dplyr`, `tidyr`, `readr`, `stringr`, `ggplot2`, `patchwork`, `vegan`

## 3. Analysis Steps and Functions

### Step 1 — Input validation and metadata loading

The script checks that all required input files are available and loads sample metadata. Metadata are standardized to ensure that sample IDs and biological groups can be matched to CIBERSORTx outputs.

### Step 2 — CIBERSORTx fraction processing

The CIBERSORTx fraction table is filtered by deconvolution p-value. Technical columns such as RMSE and correlation are removed. LM22 immune cell subsets are then collapsed into broader immune classes:

- B cells
- CD4 T cells
- CD8 T cells
- Macrophages
- Dendritic cells
- Mast cells
- NK cells

### Step 3 — Group-level immune composition

The script computes mean immune cell abundance for each biological group and generates a stacked relative composition barplot.

### Step 4 — Absolute immune score analysis

The optional absolute score file is processed to estimate global immune infiltration per sample. The script then generates a combined figure with:

### Step 5 — PERMANOVA analysis and dispersion control

The aggregated immune composition matrix is tested using Bray-Curtis PERMANOVA. The global test evaluates whether immune composition differs across all groups. Pairwise PERMANOVA then tests each group comparison separately, with Benjamini-Hochberg FDR correction. The script tests whether pairwise PERMANOVA differences could be driven by unequal within-group dispersion rather than shifts in immune composition centroids.

## 4. Output Files

* aggregated metrics per class
* permanova pairwise statistics

## 5. Macroscopic Biological Interpretation

This workflow evaluates how immune infiltrate composition changes across healthy endometrium, decidua, EC-adjacent healthy tissue, and endometrial cancer.
The relative composition plot describes how major immune cell families are redistributed across groups. The absolute score estimates the overall immune infiltration burden. PERMANOVA evaluates whether the global immune landscape differs between biological states, while pairwise PERMANOVA identifies which contrasts drive this separation. 

## `network_analysis_core.R`

### 1. Input Data

#### 1.1 Expression Matrices

* Gene-level TPM matrix used for correlation-based network construction.
* Raw count matrix used for differential expression analyses.

#### 1.2 Sample Metadata

* Clinical and biological annotation table used to define comparison groups.

#### 1.3 Precomputed Graphs (optional)

* Condition-specific and differential `igraph` objects stored as `.rds` files.

---

### 2. Required Packages

**Data Manipulation:** `dplyr`, `tibble`, `stringr`, `purrr`, `tidyr`
**Data Import:** `readxl`, `readr`
**Network Analysis:** `igraph`, `networktools`
**Functional Enrichment:** `clusterProfiler`, `msigdbr`
**Visualization:** `ggplot2`, `ComplexHeatmap`
**Export:** `openxlsx`

---

### 3. Analysis Steps and Functions

#### Step 1 — Data Import and Harmonization

Loads expression matrices, metadata, and network objects; standardizes identifiers and verifies consistency.

#### Step 2 — Differential Co-Expression Graph Construction

Builds correlation-based networks and differential rewiring graphs using predefined statistical thresholds.

#### Step 3 — Single Network Topological Analysis

Computes:

* Largest connected component
* Degree and strength
* Centrality measures
* Louvain communities
* Global topology metrics

#### Step 4 — Module Prioritization

Ranks modules using a composite score based on:

* Size
* Density
* Connectivity
* Rewiring burden

Assigns descriptive labels (e.g., Dense + Rewired, Hub-Dominated).

#### Step 5 — Functional Enrichment

Performs GO and Hallmark enrichment to assign biological meaning to each module.

#### Step 6 — Differential Rewiring Analysis

Quantifies:

* Edge overlap (Jaccard)
* Hub preservation
* Gain/loss balance
* Positive/negative rewiring fractions

#### Step 7 — Cross-Network Comparison

Compares topology and module structure across conditions.

### 4. Output

* Topology tables
* Module ranking tables
* Rewiring metrics
* Enrichment results

### 5. Macroscopic Biological Interpretation

This script reconstructs transcriptional regulatory architectures and quantifies how gene–gene coordination is reorganized across biological conditions. It identifies preserved and disrupted modules, emergent hubs, and condition-specific regulatory programs, enabling a systems-level comparison between healthy endometrium, decidua, and endometrial cancer.

---

# `network_rewire_context.R`

### 1. Input Data

* Top rewired gene tables
* Module ranking tables from case and control networks

### 2. Required Packages

`dplyr`, `readr`, `readxl`, `stringr`

### 3. Analysis Steps and Functions

#### Step 1 — Load Rewired Genes

Reads genes ranked by rewiring magnitude.

#### Step 2 — Join Module Context

Associates each gene with:

* Case module
* Control module
* Module scores and ranks

#### Step 3 — Contextual Annotation

Defines transitions such as:

* Low → High priority
* Stable core
* Module switching

### 4. Output

* Contextualized rewired gene tables
* Tier transition summaries

### 5. Macroscopic Biological Interpretation

Adds systems-level context to rewired genes by identifying whether they migrate between weak and highly prioritized modules, highlighting candidate regulatory switches associated with disease or physiological adaptation.

---

# `plot_module_metagraph_circlize_hallmark.R`

### 1. Input Data

* Module assignments
* Hallmark enrichment results
* Module meta-graph

### 2. Required Packages

`circlize`, `igraph`, `dplyr`, `ComplexHeatmap`

### 3. Analysis Steps and Functions

#### Step 1 — Functional Categorization

Groups Hallmark pathways into broad biological classes.

#### Step 2 — Circular Visualization

Displays:

* Modules as sectors
* Hallmark annotations as outer bars
* Inter-module links

### 4. Output

* Circular module metagraphs with functional annotations

### 5. Macroscopic Biological Interpretation

Integrates network topology and pathway enrichment into a single systems-level visualization, revealing how coordinated transcriptional modules connect to specific biological functions and how these functions are reorganized across conditions.


## `Spatial_ICI.py`

### 1. Purpose

This Python script generates publication-ready spatial plots for two immune-inhibitory axes:

1. **CD47/SIRPA macrophage axis**
2. **NECTIN2/TIGIT lymphoid axis**

It processes one or more `.h5ad` spatial transcriptomics objects and produces per-sample tables, summary plots, spatial zoom-ins, and signature-based comparisons.

### 2. Input Data

#### 2.1 Spatial transcriptomics object

Each input file must be an `.h5ad` object containing:

* gene expression matrix
* spatial coordinates in `adata.obsm["spatial"]` or equivalent `obs` columns
* cell-type annotations in `adata.obs[label_col]`

#### 2.2 Required genes

The script specifically evaluates:

* `EPCAM`
* `CD47`
* `SIRPA`
* `NECTIN2`
* `TIGIT`

#### 2.3 Cell labels

The expected cell types are:

* Endometrial cells
* Macrophages
* NK cells
* T cells

### 3. Required Packages

`scanpy`, `numpy`, `pandas`, `scipy`, `sklearn`, `matplotlib`

### 4. Main Analysis Steps

#### Step 1 — Load sample

Optional filtering removes cells with fewer than a minimum number of transcripts.

#### Step 2 — Extract spatial coordinates

These coordinates are used to build local spatial neighborhoods.

#### Step 3 — Extract gene expression and positivity calls

The function `add_expression_and_signature_calls()` extracts expression for:

* `EPCAM`
* `CD47`
* `SIRPA`
* `NECTIN2`
* `TIGIT`

A cell is considered positive if expression is greater than the chosen transcript threshold.

#### Step 4 — Score macrophage inhibitory signature

This estimates the inhibitory or anti-phagocytic state of macrophages.

#### Step 5 — Build spatial neighborhood graph

The function `build_spatial_graph()` builds a k-nearest-neighbor graph from spatial coordinates.

This graph is used to determine whether immune cells are spatially close to ligand-positive endometrial cells.

#### Step 6 — Define cell masks

The function `build_cell_masks()` identifies:

* endometrial cells
* macrophages
* NK cells
* T cells

These masks drive all downstream cell-type-specific analyses.

#### Step 7 — CD47/SIRPA analysis

The script quantifies:

* EPCAM/CD47 status in endometrial cells
* SIRPA status in macrophages
* macrophages spatially exposed to CD47-positive endometrial cells
* inhibitory macrophage signature across CD47/SIRPA groups

Main output groups:

```text
CD47+ / SIRPA+
CD47+ / SIRPA-
CD47- / SIRPA+
CD47- / SIRPA-
```

#### Step 8 — NECTIN2/TIGIT analysis

The script quantifies:

* EPCAM/NECTIN2 status in endometrial cells
* TIGIT status in NK cells
* TIGIT status in T cells
* NK/T cells spatially exposed to NECTIN2-positive endometrial cells

Main output groups:

```text
TIGIT+ / NECTIN2+
TIGIT+ / NECTIN2-
TIGIT- / NECTIN2+
TIGIT- / NECTIN2-
```

#### Step 9 — Transcript-depth trajectory analysis

The script compares total transcript counts across marker-defined groups.

This is used as a quality-aware readout to verify whether observed positivity patterns are associated with transcript depth or represent biologically meaningful spatial states.

#### Step 10 — Spatial zoom-ins

The function `plot_cd47_sirpa_high_signature_zoom()` identifies macrophages that are:

* `SIRPA+`
* spatially exposed to `CD47+` endometrial cells
* high-scoring for the negative phagocytosis signature

It then generates spatial crop plots centered on these macrophages.

### 5. Output Files

#### Tables

The `files/` folder contains CSV tables such as:

* signature genes present
* cell-type composition
* endometrial EPCAM/CD47 table
* endometrial EPCAM/NECTIN2 table
* macrophage SIRPA table
* macrophage signature table
* NK TIGIT table
* T-cell TIGIT table
* transcript-depth tables

#### Figures

The `plots/` folder contains:

* CD47/SIRPA publication panel
* NECTIN2/TIGIT publication panel
* CD47/SIRPA signature violin plot
* CD47/SIRPA high-signature spatial zoom-ins

### 6. Macroscopic Biological Interpretation

The pipeline tests whether immune-inhibitory ligand–receptor axes are spatially organized in endometrial tissue.

The **CD47/SIRPA axis** evaluates whether endometrial cells expressing CD47 are positioned near SIRPA-positive macrophages and whether those macrophages show increased anti-phagocytic/inhibitory transcriptional programs.

The **NECTIN2/TIGIT axis** evaluates whether NECTIN2-positive endometrial cells are spatially associated with TIGIT-positive NK or T cells, suggesting local lymphoid inhibitory signaling.

Overall, the script provides a spatial framework to identify immune-tolerance niches and compare how these niches may differ between physiological invasion and pathological tumor-associated immune evasion.
