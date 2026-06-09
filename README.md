# LogicComm

**Logic-Based Cell-Cell Communication Analysis for Multi-Sample scRNA-seq**

LogicComm converts single-cell expression matrices into **Relative Expression Ordering (REO)** binary logic states, preserves optional within-cell rank evidence, and scores ligand-receptor communication with **Logic Consensus Scores (LCS)**. It is designed for single-sample exploration and cohort-level **Case vs Control** comparisons.

---

## What LogicComm Does

1. **REO conversion**: each gene is marked active/inactive within each cell by comparing it with that cell's expression anchor.
2. **LCS scoring**: each ligand-receptor pair receives a score based on ligand-active sender cells and receptor-active receiver cells.
3. **Rank-aware REO evidence**: rank dominance, margin, threshold stability, and cell-type-pair specificity complement binary active calls.
4. **Multi-sample comparison**: active L-R pairs are compared across N Case vs N Control samples using frequency voting, continuous rank evidence, and statistical tests.
5. **Per-cell activity scoring**: cells can be ranked by sender, receiver, or combined communication potential.
5. **Cell-type-resolved interpretation**: summarize celltypeA -> celltypeB communication counts, logic strength, pathway drivers, and signaling roles.
7. **Publication-level diagnostics**: spatial neighborhoods, pseudotime dynamics, sample-level GLMs, QC plots, and markdown report helpers.

| Feature | LogicComm behavior |
|---|---|
| Input orientation | genes x cells |
| Input types | base matrix, sparse `Matrix`, Seurat object, optional BPCells IterableMatrix |
| Cell-type annotation | not required for global LCS; optional for cell-type-resolved summaries |
| KNN/SNN graph | optional but recommended; binary or weighted graph scoring supported |
| Multimer complexes | all listed subunits must be present and active |
| Large sparse matrices | all genes used for anchors by default; retained genes processed in chunks |
| Interpretability | cell-type roles, receiver response, bootstrap/null diagnostics, spatial/pseudotime modules, sample-level GLMs, and report helpers |

---

## Installation

### Core dependencies

```r
install.packages(c(
  "Matrix", "dplyr", "ggplot2", "ggrepel", "scales", "Rcpp",
  "testthat", "rmarkdown"
))
```

### Optional dependencies

```r
# Required only for heatmap plotting
install.packages("pheatmap")

# Required only for Seurat object input / Seurat KNN extraction
install.packages(c("Seurat", "SeuratObject"))

# Optional for graph centrality acceleration and publication diagnostics
install.packages("igraph")

# Required only for BPCells on-disk matrix input; install according to BPCells docs
# install.packages("BPCells") or use the official installation method
```

### Install LogicComm from a local source tarball

```r
install.packages("LogicComm_0.10.4.tar.gz", repos = NULL, type = "source")

# Or install a published release directly from GitHub:
# install.packages(
#   "https://github.com/awaragml00029-debug/LogicComm/releases/download/v0.10.4/LogicComm_0.10.4.tar.gz",
#   repos = NULL, type = "source")
```

If installing from the unpacked source directory during development:

```sh
R CMD INSTALL path/to/LogicComm
```

---

## Rank-Aware REO Evidence

LogicComm 0.7 adds an experimental continuous evidence layer inspired by rank-order methods. The original binary REO/LCS functions remain unchanged, but users can now keep the rank margin that is lost when every above-threshold gene becomes a simple active call.

```r
reo <- calc_REO_matrix(count_matrix, lr_genes = all_lr_genes(lr_pairs_human),
                       return_rank = TRUE, verbose = FALSE)
rank_ev <- IdentifyRankLogicConsensus(reo, lr_db = lr_pairs_human,
                                      threshold_grid = c(0.4, 0.5, 0.6, 0.7))
head(rank_ev[order(-rank_ev$specificity_weighted_rank_lcs), ])
```

For multi-sample designs, use `CompareRankLogicGroups()` on one rank-evidence table per biological sample. This is most useful when a permissive binary threshold makes Case and Control both active, but Case has stronger rank dominance and better threshold stability.
## Quick Start: Single Sample

```r
library(LogicComm)
data(lr_pairs_human)

# Expression input must be genes x cells.
# Supported examples:
# - count_matrix: base matrix or Matrix sparse matrix
# - my_seurat: Seurat object

lr_genes <- all_lr_genes(lr_pairs_human)
# By default, anchors are computed from all genes in count_matrix;
# lr_genes only controls which rows are retained in the output REO matrix.

reo <- calc_REO_matrix(
  count_matrix,
  lr_genes        = lr_genes,
  rank_threshold  = 0.5,
  chunk_size      = 5000,
  gene_background = "quantile",   # cell-type-free ambient guard (recommended)
  verbose         = TRUE
)

# Global mode: no KNN graph required.
lcs_global <- IdentifyLogicConsensus(reo, verbose = FALSE)

# View top pairs.
sort(lcs_global, decreasing = TRUE)[1:20]
```

> **Ambient RNA guard (`gene_background`).** The REO anchor is *within-cell*: a
> gene is active if it exceeds that cell's own median gene. A marker such as
> `CD8A` that leaks into non-CD8 cells as ambient RNA (a few counts, still above
> the cell's low median) is therefore wrongly called active, which can surface
> biologically impossible axes (for example an `HLA-B -> CD8A` interaction with a
> Treg receiver). `gene_background = "quantile"` additionally requires each gene
> to exceed its own background, taken over the cells that *detect* it (nonzero
> counts) so the guard stays effective even for markers expressed in a minority
> of cells. This removes ambient activity **without any cell-type annotation**.
> It is `"none"` by default for backward compatibility but recommended for
> biological interpretation; for rigorous decontamination run SoupX, DecontX, or
> CellBender upstream.

---

## What LogicComm scores (v0.12+)

LogicComm scores cell-cell communication at the **cell-type level from REO
co-expression** — it no longer aggregates over a per-cell KNN/SNN neighborhood
graph. For dissociated scRNA-seq that graph lives in expression space
(transcriptomic similarity), not physical space, so it cannot license spatial
juxtacrine/paracrine claims; cell-type co-expression is the honest,
sample-comparable unit. See the cell-type workflow below, and
`discover_celltype_communication()` for the recommended one-call pipeline.

Legacy neighborhood arguments (`knn_mat`, `mode`, `graph_name`,
`remove_self_edges`, `graph_symmetrize`, `edge_weight_mode`) are still accepted
for backward compatibility but are **deprecated and ignored with a warning**.

---

## PBMC3K Four-Sample Validation Demo

A publication-style PBMC validation walkthrough is included in `vignettes/PBMC_multisample_demo.Rmd`. It splits PBMC3K into `pbmc1`, `pbmc2`, `pbmc3`, and `pbmc4`, adds a known Case perturbation to `pbmc1` and `pbmc2`, and verifies that `run_multisample()` recovers the expected Case-enriched MIF, galectin, and MHC-I LR axes at the sample level.

```r
source(system.file("examples", "pbmc_multisample_demo.R", package = "LogicComm"))
demo <- run_pbmc_multisample_demo()
demo$target_rows[, c("lr_pair", "case_freq", "ctrl_freq", "asymmetry", "log2fc_lcs")]
```

## Quick Start: Multi-Sample Case vs Control

### One-stop pipeline

```r
library(LogicComm)
data(lr_pairs_human)

# sample_list: named list of genes x cells matrices or Seurat objects.
# Example names: Case_1, Case_2, Ctrl_1, Ctrl_2

group_info <- c("Case", "Case", "Ctrl", "Ctrl")
names(group_info) <- names(sample_list)

result <- run_multisample(
  sample_list,
  group_info     = group_info,
  lr_db          = lr_pairs_human,
  rank_threshold = 0.5,
  lcs_threshold  = 0.01,
  case_label     = "Case",
  ctrl_label     = "Ctrl",
  chunk_size     = 5000,
  mc_cores       = 1,
  verbose        = TRUE
)

print(result$comparison, n = 10)
```

### Manual multi-sample workflow

```r
lr_genes <- all_lr_genes(lr_pairs_human)

lcs_list <- lapply(sample_list, function(x) {
  reo <- calc_REO_matrix(x, lr_genes = lr_genes, verbose = FALSE)
  IdentifyLogicConsensus(reo, verbose = FALSE)
})

result <- CompareLogicGroups(
  lcs_list,
  group_info     = group_info,
  case_label     = "Case",
  ctrl_label     = "Ctrl",
  lcs_threshold  = 0.01,
  verbose        = TRUE
)
```

Notes:

- `group_info` must contain labels for every sample in `lcs_list`.
- L-R pairs are aligned across all samples by name; missing pairs are stored as `NA`.
- Missing LCS values are treated as unavailable, not as negative calls; per-pair frequencies use non-missing sample counts.
- Pass the same custom `lr_db` to `CompareLogicGroups()` to keep ligand, receptor, pathway, and annotation metadata.
- On Windows, `mc_cores > 1` does not use fork-based `mclapply`; use `mc_cores = 1` unless another parallel backend is implemented.

---

## Visualization

```r
# Bubble plot: Case frequency vs Control frequency.
plot_lcs_bubble(
  result$comparison,
  top_n         = 25,
  fdr_cutoff    = 0.05,
  min_asymmetry = 0.2,
  color_by      = "log2fc"   # "log2fc" or "asymmetry"
)

# Heatmap requires pheatmap.
plot_lcs_heatmap(
  result$comparison,
  group_info = group_info,
  top_n      = 40,
  fdr_cutoff = 0.1
)
```


---

**Publication figure system.** Every `plot_*` function shares one colorblind-safe,
print-safe theme (`theme_logiccomm()`) and brand palette. Restyle custom layers with
the brand scales (`scale_color/fill_logiccomm_c()`, `_diverging()`, `_d()`,
`scale_*_tier()`) and `logiccomm_brand` colours, and export at journal column widths
(single 89 mm / double 183 mm; vector PDF via cairo) with `save_logiccomm_figure()`.
Compose multi-panel figures with `patchwork`. See the Seurat tutorial's
"Publication figures" section for the full reference.

## Biology-Oriented Cell-Type Communication Visualization

Once a cell-type communication result is available, LogicComm can summarize it
around concrete biological questions:

- **How many active interactions exist between celltypeA and celltypeB?** Use
  `n_active_lr` / `active_lr_event_count`, the number of active
  `sender -> receiver -> LR pair` events.
- **How strong is the communication from celltypeA to celltypeB?** Use `sum_lcs`,
  the summed active LCS across L-R pairs. This is a logic-consensus strength, not
  raw ligand/receptor expression.
- **How much cell-cell edge evidence supports a result?** Inspect `n_edges`,
  `edge_weight_sum`, `n_active_edges`, and `sum_active_edges`.
- **What role does each cell type play in the communication network?** Use
  `role_summary`, which reports sender/out-degree, receiver/in-degree,
  mediator/betweenness, and influencer/information-flow scores with low-
  communication filtering and role confidence.

```r
# cell_labels can be a named vector, a factor, Seurat Idents, or a metadata field.
cell_labels <- SeuratObject::Idents(my_seurat)

ct_comm <- summarize_celltype_communication(
  reo_mat       = reo,
  cell_labels   = cell_labels,        # no graph needed: scored from cell-type co-expression
  lr_db         = lr_pairs_human,
  lcs_threshold = 0.01,
  min_edges     = 20                  # min sender x receiver cell-count opportunity
)

ct_comm$pair_summary[order(ct_comm$pair_summary$sum_lcs, decreasing = TRUE), ][1:10, ]
ct_comm$role_summary[order(ct_comm$role_summary$hub_score, decreasing = TRUE), ]
```

### Interaction quantity, communication strength, and evidence support

```r
# Number of active L-R events for each sender -> receiver cell-type pair.
plot_celltype_heatmap(ct_comm, metric = "n_active_lr")

# Summed logic consensus strength for each sender -> receiver cell-type pair.
plot_celltype_heatmap(ct_comm, metric = "sum_lcs")

# How many active cell-cell edge observations support active L-R events?
plot_celltype_heatmap(ct_comm, metric = "sum_active_edges")
```

Interpretation: `n_active_lr` highlights broad communication richness, whereas
`sum_lcs` highlights total logic-supported strength. A cell-type pair with many
weak L-R events can therefore differ from a pair with a few strong, recurrent
events. `sum_active_edges` is an evidence-support metric and helps avoid
interpreting high LCS values from very small edge sets.

### Signaling roles of cell types

```r
plot_celltype_roles(ct_comm)
plot_celltype_role_heatmap(ct_comm)
plot_celltype_role_radar(ct_comm, cell_types = c("DC", "CD8 T", "CD14+ Mono"))
plot_celltype_network(ct_comm, metric = "sum_lcs", top_n_edges = 30)
```

Role metrics:

| Role name | Recommended column | Backward-compatible alias | Biological interpretation |
|---|---|---|---|
| Sender / out-degree | `sender_role_score` | `outdegree_score` | Cell type mainly sends ligand programs to other cell types |
| Receiver / in-degree | `receiver_role_score` | `indegree_score` | Cell type mainly receives receptor stimulation from other cell types |
| Mediator / betweenness | `mediator_role_score` | `betweenness_score` | Cell type bridges communication paths between other cell types |
| Influencer / information flow | `influencer_role_score` | `information_score` | Cell type can reach many others through efficient directed communication paths |

The role table now also distinguishes several counts that are often confused:

| Field | Meaning |
|---|---|
| `outgoing_lr_event_count`, `incoming_lr_event_count` | Number of active `sender -> receiver -> LR pair` events |
| `outgoing_unique_lr_count`, `incoming_unique_lr_count` | Number of unique LR pairs involving the cell type |
| `outgoing_target_type_count`, `incoming_source_type_count` | Number of receiver/source cell types connected by active communication |
| `outgoing_active_edge_support`, `incoming_active_edge_support` | Total active edge observations supporting active L-R events |

`dominant_role` includes a `Low-communication` category when the hub score or
active event count is too small. From v0.6.1 onward, role interpretation is split
into three safer concepts:

| Field | Interpretation |
|---|---|
| `role_separation_label` | Whether one of sender/receiver/mediator/influencer clearly dominates the other role axes |
| `communication_evidence_label` | Whether the cell type has enough hub strength, LR events, active-edge support, and connected source/target types |
| `role_reliability_label` / `dominant_role_strict` | A conservative label that combines topology and evidence |

This distinction matters. A B cell can have a clear topological tendency, for
example `Mediator (betweenness)`, but only moderate communication evidence. In
that case the recommended interpretation is not "B cells are definitive
gatekeepers"; it is:

> B cells are **candidate mediator-like bridges** in the inferred PBMC
> communication graph. The role is topology-supported, but biological confidence
> depends on the total edge support, the specificity of the L-R axes, and whether
> the result reproduces across samples or spatial neighborhoods.

Use the helper below when writing biological results:

```r
interpret_celltype_roles(ct_comm)
plot_celltype_role_dotplot(ct_comm)
```

### Explain a specific cell-type pair or pathway

```r
plot_lr_bubble_by_celltype(
  ct_comm,
  sender   = "Macrophage",
  receiver = "Endothelial",
  top_n    = 30,
  color_by = "pathway"
)

plot_lr_activity_balance(ct_comm, sender = "Macrophage", receiver = "Endothelial")
plot_pathway_heatmap(ct_comm, metric = "sum_lcs", top_n_pathways = 25)
plot_lr_evidence(ct_comm, sender = "Macrophage", receiver = "Endothelial", lr_pair = "VEGFA_KDR")

ex <- explain_celltype_interaction(
  ct_comm,
  sender     = "Macrophage",
  receiver   = "Endothelial",
  lr_pair    = "VEGFA_KDR",
  reo_mat    = reo,
  seurat_obj = my_seurat
)
ex$evidence
ex$plot
```

### Add receiver downstream response evidence

L-R logic alone measures **communication potential**. A stronger biological claim
asks whether the receiver also shows downstream response logic. Provide a simple
response database with `lr_pair` and `response_genes` / `target_genes` / `genes`:

```r
response_db <- data.frame(
  lr_pair = c("IFNG_IFNGR1_IFNGR2", "TGFB1_TGFBR1_TGFBR2"),
  response_genes = c("STAT1;IRF1;CXCL10", "SMAD7;SERPINE1;COL1A1")
)

ct_comm <- add_receiver_response_score(
  ct_comm,
  reo_mat      = reo,
  response_db  = response_db,
  response_mode = "any"
)

ct_comm$lr_table[order(ct_comm$lr_table$response_integrated_score, decreasing = TRUE), ][1:10, ]
```

### Recommended: one-call discovery (`discover_celltype_communication`)

The footgun-free entry point — cell-type co-expression scoring, specificity +
proliferation-confound annotation, an **axis-level** permutation null, evidence
ranking, a confound-filtered view, and an FDR-passing shortlist, in one call:

```r
# REO with the rank matrix (needed only for lcs_weighting = "rank").
reo <- calc_REO_matrix(count_matrix, lr_genes = lr_genes, return_rank = TRUE)

res <- discover_celltype_communication(
  reo_mat       = reo,
  cell_labels   = cell_labels,
  lr_db         = lr_pairs_human,
  expr          = count_matrix,   # enables the proliferation/breadth-hub filter
  lcs_weighting = "rank",         # REO-intensity score (optional; default "binary")
  n_perm        = 1000,
  n_cores       = 8,
  seed          = 1,
  fdr_cutoff    = 0.1
)
res$shortlist   # axes passing permutation_fdr <= 0.1, ranked
res$view        # confound-filtered discovery view
res$ranked      # all active axes with per-axis permutation_empirical_p / permutation_fdr
```

### Uncertainty, null models, and threshold sensitivity

The individual steps that `discover_celltype_communication()` wraps:

```r
# Per-cell bootstrap for confidence intervals (resamples cells, not the
# sender x receiver cell-count product).
boot_lr <- bootstrap_celltype_communication(ct_comm, n_boot = 200, level = "celltype_lr")
boot_pair <- bootstrap_celltype_communication(ct_comm, n_boot = 200, level = "celltype_pair")

# Axis-level cell-label permutation null (default metric = "lcs"): every
# sender -> receiver -> L-R axis gets its OWN empirical p and BH FDR. The
# co-expression null needs no graph, so a publication-grade n_perm is cheap;
# n_cores parallelizes on Unix/macOS.
null_pair <- permute_celltype_communication(
  ct_comm = ct_comm,
  reo_mat = reo,
  n_perm  = 1000,
  n_cores = 8,
  seed    = 1
)

diagnose_permutation_resolution(null_pair)

# If null_sd == 0 and observed > null_mean, LogicComm now reports
# z_score = Inf and degenerate_positive_null = TRUE. Biologically, this means
# the shuffled cell-label null never generated comparable communication. Treat it
# as a structural-null signal, not as an ordinary Gaussian z-score; verify edge
# support, cell abundance, specificity, and sample reproducibility.

# Score whether LR axes are pair-specific or broad/identity-associated.
ct_comm <- score_communication_specificity(ct_comm)
ct_comm$specificity_summary[1:10, ]
plot_pathway_dominance(ct_comm)
# Check whether key discoveries are robust to REO threshold choice.
sens <- sensitivity_REO_threshold(
  my_seurat,
  rank_threshold_grid = c(0.4, 0.5, 0.6, 0.7),
  lr_db       = lr_pairs_human,
  cell_labels = cell_labels
)
sens$stability[1:20, ]

# Stability + specificity landscape: high-stability/high-specificity axes are
# stronger mechanistic candidates; high-stability/low-specificity axes often
# represent broad immune-recognition programs such as MHC/HLA.
plot_specificity_stability(sens, ct_comm = ct_comm)

# Publication-level diagnostics and evidence-aware ranking.
diag <- diagnose_celltype_communication(ct_comm, null_pair = null_pair, sens = sens)
diag

ranked_axes <- rank_communication_axes(ct_comm, null_pair = null_pair, sens = sens)
ranked_axes[1:20, ]
```

For multi-sample studies, compute `ct_comm` **within each biological sample** and
then convert cell-type-level features to named vectors for sample-level group
comparison:

```r
ct_lcs_list <- lapply(names(sample_list), function(sname) {
  reo_i <- calc_REO_matrix(sample_list[[sname]], lr_genes = lr_genes, verbose = FALSE)
  ct_i <- summarize_celltype_communication(
    reo_i,
    cell_labels      = sample_cell_labels[[sname]],
    lr_db            = lr_pairs_human,
    lcs_threshold    = 0.01,
    min_edges        = 20,
    verbose          = FALSE
  )
  celltype_comm_to_lcs(ct_i, level = "celltype_lr", metric = "lcs")
})
names(ct_lcs_list) <- names(sample_list)

ct_comparison <- CompareLogicGroups(
  ct_lcs_list,
  group_info    = group_info,
  case_label    = "Case",
  ctrl_label    = "Ctrl",
  lcs_threshold = 0.01
)

plot_differential_celltype_heatmap(ct_comparison, metric = "asymmetry")
plot_differential_celltype_volcano(ct_comparison, x = "log2fc_lcs")
```

Important interpretation limits:

- LCS is a binary REO logic-consensus proportion, not expression intensity.
- The direction celltypeA -> celltypeB comes from ligand/receptor roles; if a
  graph is symmetrized, directionality still comes from ligand and receptor
  assignment, not from graph storage direction.
- Small cell types can have too few edges; use `min_edges`, inspect `n_edges`,
  and consider bootstrap/permutation results.
- In Case vs Control designs, the statistical unit should be the sample, not the
  cell.

---

## Per-Cell Communication Activity

```r
scores <- score_lr_activity(
  reo,
  lr_db     = lr_pairs_human,
  mode      = "both",     # "sender", "receiver", or "both"
  aggregate = TRUE,
  verbose   = TRUE
)

# Add aggregate score to Seurat metadata if using Seurat.
my_seurat$comm_score <- scores$comm_score

# Rank top communication-active cells.
top_cells <- rank_comm_cells(scores, n = 50)
```

If you call `score_lr_activity(..., mode = "sender")` or `mode = "receiver"`, no combined `comm_score` is produced. `rank_comm_cells()` requires `comm_score`, so use `mode = "both", aggregate = TRUE` before ranking cells.

---

## Rank Shift Analysis

Rank shift detects genes whose within-cell expression rank changes between Case and Control, even when absolute fold-change is modest.

```r
rs <- calc_rank_shift(
  sample_list,
  group_info = group_info,
  genes      = all_lr_genes(lr_pairs_human),
  case_label = "Case",
  ctrl_label = "Ctrl"
)

head(rs)
plot_rank_shift(rs)
```

---

## Ligand-Receptor Database Provenance

`data(lr_pairs_human)` loads a compact human LR database with 109 interactions.
It was compiled as a practical seed set from public ligand-receptor resources,
including CellChatDB, CellPhoneDB, and NicheNet, then normalized into the
LogicComm schema.

Important provenance notes:

- **CellChatDB source**: CellChat describes CellChatDB as a manually curated
  database of literature-supported ligand-receptor interactions. In the
  `jinworks/CellChat` repository, the human and mouse databases are distributed
  under `data/` as `CellChatDB.human.rda` and `CellChatDB.mouse.rda`.
- **CellChatDB structure**: CellChatDB is a list-like object with components such
  as `interaction`, `geneInfo`, `complex`, and `cofactor`. The interaction table
  commonly contains fields such as ligand, receptor, pathway name, annotation,
  interaction name, and cofactor/complex references.
- **LogicComm conversion**: LogicComm does not bundle the full CellChatDB object.
  `lr_pairs_human` is a smaller, package-ready table that keeps one row per
  interaction and adds explicit `ligand_genes` / `receptor_genes` list columns.
  These list columns are what LogicComm uses for multimer logic: all listed
  subunits must be present and active for that ligand or receptor complex to be
  active.
- **Custom database support**: If you want to use the full CellChatDB, CellPhoneDB,
  LIANA, OmniPath, or another source, convert it to the schema shown below.

---

## Custom Ligand-Receptor Database

### Convert CellChatDB to LogicComm schema

If CellChat is installed, you can convert the full CellChat human database into
LogicComm format:

```r
# Depending on your CellChat installation, this object may be available through
# data(CellChatDB.human, package = "CellChat") or after library(CellChat).
data(CellChatDB.human, package = "CellChat")

cellchat_lr <- as_logiccomm_lr_db_from_cellchat(
  CellChatDB.human,
  include_modulators = TRUE
)

lr_genes <- all_lr_genes(cellchat_lr, include_modulators = TRUE)
reo <- calc_REO_matrix(count_matrix, lr_genes = lr_genes, return_rank = TRUE)

# Base LogicComm LCS still works and ignores modulators.
lcs <- IdentifyLogicConsensus(reo$logic, lr_db = cellchat_lr)

# Rank-aware Logic Gate LCS uses CellChat agonist/antagonist/co-receptor fields.
gate_lcs <- IdentifyLogicGateConsensus(
  reo,
  lr_db = cellchat_lr,
  positive_gate = "all",
  negative_gate = "rank_block"
)

# Inspect modulator activity without changing LCS.
mod_summary <- summarize_lr_modulators(reo$logic, cellchat_lr, rank_mat = reo$rank)
```

The converter resolves CellChat complex names through `CellChatDB.human$complex`.
For example, if an interaction uses a receptor complex name, LogicComm expands it
into the corresponding receptor subunits in `receptor_genes`. Those subunits are
then evaluated with LogicComm's all-subunits-active rule.

From version 0.4.1, the converter can preserve CellChatDB modulators:
`agonist`, `antagonist`, `co_A_receptor`, and `co_I_receptor`. LogicComm does not
try to reproduce CellChat's probability model. Instead, the optional
`IdentifyLogicGateConsensus()` function translates these fields into REO logic
gates: positive modulators are configurable gates (`positive_gate = "all"` or `"any"`), and negative modulators can block a
communication event when their within-cell rank overtakes the ligand or receptor.
This gives a rank-order explanation such as: ligand and receptor are still active,
but antagonist rank reversal blocks the effective communication logic.

To check overlap between a converted CellChatDB and LogicComm's built-in compact
seed set:

```r
logic_key <- function(db) {
  vapply(seq_len(nrow(db)), function(i) {
    lig <- paste(sort(unique(db$ligand_genes[[i]])), collapse = "+")
    rec <- paste(sort(unique(db$receptor_genes[[i]])), collapse = "+")
    paste(lig, rec, sep = "__")
  }, character(1))
}

cellchat_lr <- as_logiccomm_lr_db_from_cellchat(CellChatDB.human)
built_in_key <- logic_key(lr_pairs_human)
cellchat_key <- logic_key(cellchat_lr)

mean(built_in_key %in% cellchat_key)
setdiff(lr_pairs_human$lr_pair, lr_pairs_human$lr_pair[built_in_key %in% cellchat_key])
```

### Manual custom table

A custom LR database should contain at least:

- `lr_pair`: unique pair name
- `ligand_genes`: list column of ligand subunit gene vectors
- `receptor_genes`: list column of receptor subunit gene vectors

Recommended metadata columns:

- `ligand`
- `receptor`
- `pathway`
- `annotation`

Example:

```r
custom_lr <- data.frame(
  lr_pair    = "TGFB1_TGFBR1_TGFBR2",
  ligand     = "TGFB1",
  receptor   = "TGFBR1+TGFBR2",
  pathway    = "TGFb",
  annotation = "custom",
  stringsAsFactors = FALSE
)
custom_lr$ligand_genes <- list("TGFB1")
custom_lr$receptor_genes <- list(c("TGFBR1", "TGFBR2"))

lr_genes <- all_lr_genes(custom_lr)
reo <- calc_REO_matrix(count_matrix, lr_genes = lr_genes)
lcs <- IdentifyLogicConsensus(reo, lr_db = custom_lr)
```

Multimer rule: all subunits in `ligand_genes` or `receptor_genes` must be present in `reo` and active in a cell for the complex to be active.

---

## Publication-Level Extensions

LogicComm 0.7 extends several modules intended for manuscript-grade analyses.
They do not change the core REO/LCS definition; instead they add stronger
biological context, statistical modeling, and figure/report scaffolding.

### Spatial communication (graph scoring pending in v0.12)

> **Status.** With genuine spatial coordinates a physical neighborhood graph *is*
> a legitimate, non-expression-space basis for juxtacrine/paracrine scoring. But
> the v0.12 cell-type-co-expression rewrite removed the per-cell graph scorer that
> `summarize_spatial_communication()` relied on, so spatial **graph scoring is
> temporarily unavailable**: the function now warns and returns cell-type
> co-expression (the graph is *not* used). A dedicated spatial edge scorer over
> physical coordinates is planned. `build_spatial_graph()` and
> `plot_spatial_logic()` (which visualizes REO states on coordinates) still work.

```r
# Graph construction and per-spot visualization still work:
sp_graph <- build_spatial_graph(coords, mode = "knn", k = 6, distance_weight = "gaussian")
plot_spatial_logic(reo, coords, lr_pair = "VEGFA_KDR",
                   lr_db = lr_pairs_human, cell_labels = spatial_cell_labels)
```

### Pseudotime or time-course communication dynamics

For differentiation or treatment time-course studies, bin cells along a numeric
pseudotime or use categorical time labels. Each bin is analyzed as a local
LogicComm summary.

```r
dyn <- summarize_communication_dynamics(
  reo_mat      = reo,
  pseudotime   = pseudotime_vector,
  cell_labels  = cell_labels,
  lr_db        = lr_pairs_human,
  n_bins       = 6,
  bin_method   = "quantile",
  min_cells_per_bin = 50,
  lcs_threshold = 0.01,
  min_edges     = 20
)

plot_communication_dynamics(dyn, level = "pair", metric = "sum_lcs", top_n = 8)
plot_communication_dynamics(dyn, level = "role", metric = "hub_score", top_n = 8)
```

### Sample-level GLM for differential communication

For cohorts, the safest statistical unit is the biological sample. The GLM helper
models active edge counts with available edge opportunities as trials:

```text
active_edges ~ group + covariates, family = quasibinomial
```

```r
# sample_ct_list is a named list of LogicCommCellTypeComm objects, one per sample.
sample_meta <- data.frame(
  group = c("Control", "Control", "Disease", "Disease"),
  batch = c("A", "B", "A", "B"),
  row.names = names(sample_ct_list)
)

glm_res <- fit_celltype_comm_glm(
  sample_ct_list,
  sample_metadata = sample_meta,
  design = ~ group + batch,
  coef = "groupDisease",
  level = "celltype_lr",
  min_samples = 4,
  min_total_trials = 50
)

plot_celltype_glm_volcano(glm_res)
```

Use this together with `CompareLogicGroups()`: the frequency-based comparison is
simple and robust, whereas the GLM is better when covariates, unequal edge
opportunities, or sample-level count modeling matter.

### Publication-ready reports and QC

```r
findings <- summarize_communication_findings(ct_comm, top_n = 15)
findings

plot_communication_qc(ct_comm)
plot_role_confidence(ct_comm)

write_communication_report(
  ct_comm,
  file = "LogicComm_PBMC_report.md",
  title = "PBMC LogicComm communication report",
  top_n = 15
)
```

A typical publication figure set is:

1. **Overview:** `plot_celltype_heatmap(metric = "n_active_lr")` and
   `plot_celltype_heatmap(metric = "sum_lcs")`.
2. **Network structure:** `plot_celltype_network_publication()` plus `plot_communication_qc()`.
3. **Cellular roles:** `plot_celltype_roles()`, `plot_celltype_role_dotplot()`,
   `plot_celltype_role_heatmap()`, and `plot_role_confidence()`.
4. **Specificity:** `score_communication_specificity()`, `plot_pathway_dominance()`,
   and `plot_specificity_stability()` to separate broad/identity-associated axes
   from pair-specific candidates.
5. **Mechanism:** `plot_lr_bubble_by_celltype()`, `plot_lr_activity_balance()`,
   and `plot_lr_evidence()` for a selected sender -> receiver axis.
6. **Robustness:** `bootstrap_celltype_communication()`,
   `permute_celltype_communication()`, `diagnose_permutation_resolution()`, and
   `sensitivity_REO_threshold()`.
7. **Cohort differences:** `CompareLogicGroups()` or `fit_celltype_comm_glm()`
   with heatmap/volcano visualizations.

### Common interpretation pitfalls

- LCS is a **logic-consensus proportion**, not ligand or receptor expression
  magnitude.
- KNN/SNN edges are neighborhood opportunities. They are not physical contacts
  unless you use a spatial graph.
- High LCS does not by itself prove causal signaling. Receiver response genes,
  perturbation data, spatial colocalization, or external validation strengthen
  the claim.
- Small cell types can produce unstable high scores if `n_edges` is small. Always
  inspect edge support and consider bootstrap/permutation diagnostics.
- For group comparisons, do not treat cells as independent replicates. Run
  sample-level summaries first and compare across biological samples.

- In PBMC-like data, stable MHC/HLA axes are often biologically real but broad;
  combine specificity, pathway dominance, and receiver-response evidence before
  describing them as pair-specific signaling mechanisms.
- A zero-variance positive permutation null (`degenerate_positive_null = TRUE`)
  means shuffled labels never reproduced the observed communication score. This
  is a structural-null flag, not an ordinary Gaussian z-score. Report edge
  support, cell abundance, specificity, and replication when interpreting it.


## Troubleshooting

### Error: dependencies are not available

Install missing dependencies first:

```r
install.packages(c("dplyr", "ggplot2", "ggrepel", "scales", "Rcpp"))
```

### Warning: a `knn_mat` / `mode` / `graph_*` argument is deprecated and ignored

LogicComm scores at the cell-type level from co-expression and no longer uses a
per-cell graph. Drop these arguments; they are accepted only for backward
compatibility and have no effect.

### Heatmap fails

`plot_lcs_heatmap()` requires `pheatmap`:

```r
install.packages("pheatmap")
```

### Large matrix uses too much memory

Use `lr_genes = all_lr_genes(lr_pairs_human)` and a smaller `chunk_size`. LogicComm computes anchors from all genes by default, then converts retained output genes by chunk. For very large matrices, keep the LR output set small and reduce `chunk_size`. Only use `anchor_genes` when you intentionally want anchors based on a defined gene universe.

```r
reo <- calc_REO_matrix(count_matrix, lr_genes = lr_genes, chunk_size = 1000)
```

For very large on-disk data, use BPCells input if available.

---

## End-to-End Discovery Workflow

Once a `summarize_celltype_communication()` result (`ct_comm`) is available,
these steps take you from "the package ran" to "here are the prioritized,
evidence-backed interactions, who drives them, and what changes between
conditions". (`discover_celltype_communication()` wraps steps 2 below.)

```r
# 1. Which cell types communicate, how broadly, and through what?
part <- summarize_celltype_participation(ct_comm, reo_mat = reo)
part                                            # major hubs + communicating-cell fractions
plot_celltype_participation(part)               # fraction of communicating cells per type
plot_celltype_pathway_composition(part, "outgoing")
plot_celltype_communication_profile(part, "CD14+ Mono")   # one subgroup's profile card

# 2. Prioritize effective interactions by integrated evidence.
ct_comm <- score_communication_specificity(ct_comm)
null_pair <- permute_celltype_communication(ct_comm, reo_mat = reo, n_perm = 1000, n_cores = 8)
sens      <- sensitivity_REO_threshold(my_seurat, lr_db = lr_pairs_human,
                                       cell_labels = cell_labels)
ranked <- rank_communication_axes(ct_comm, null_pair = null_pair, sens = sens)
ranked[, c("sender_type","receiver_type","lr_pair","discovery_score","evidence_tier")][1:20, ]
plot_communication_discovery(ranked)            # strength vs specificity, coloured by tier

# 3. For a cohort, where do the differences occur (which sender -> receiver pair)?
diff <- differential_celltype_communication(sample_ct_list, group_info,
                                            case_label = "Case", ctrl_label = "Ctrl")
diff                                            # top differential subgroups + FDR caveat
plot_differential_communication_summary(diff)   # differential L-R counts per subgroup
plot_differential_celltype_heatmap(diff, metric = "asymmetry")
```

`discovery_score` integrates communication strength, cell-type-pair specificity,
REO threshold stability, cell-label permutation support, and (when available)
receiver response; `evidence_tier` summarizes Tier 1 (strong, specific,
supported) to Tier 4. For cohorts, see the FDR-and-sample-size note in
`?differential_celltype_communication`: with few replicates per group Fisher FDR
saturates near 1, so prioritize by effect size and corroborate with
`fit_celltype_comm_glm()`.

---

## Core Functions

| Function | Description |
|---|---|
| `calc_REO_matrix()` | Convert expression matrix to sparse REO logic matrix |
| `IdentifyLogicConsensus()` | Compute LCS for all L-R pairs in one sample |
| `CompareLogicGroups()` | Compare LCS across Case vs Control samples |
| `run_multisample()` | One-stop REO + LCS + group comparison pipeline |
| `score_lr_activity()` | Per-cell sender/receiver/communication scores |
| `rank_comm_cells()` | Rank cells by communication potential |
| `AutoLabelLogicClusters()` | Label clusters by sender/receiver L-R activity |
| `calc_rank_shift()` | Detect genes with within-cell rank shifts |
| `plot_lcs_bubble()` | Bubble plot for group comparison results |
| `plot_lcs_heatmap()` | Heatmap of LCS across samples |
| `plot_umap_logic()` | UMAP colored by a specific L-R logic state |
| `filter_lcs()` | Filter comparison results by FDR/asymmetry/direction |
| `build_spatial_graph()` | Construct spatial kNN/radius graphs from coordinates |
| `summarize_spatial_communication()` | Run cell-type communication on spatial neighborhoods |
| `summarize_communication_dynamics()` | Analyze communication along pseudotime or time bins |
| `fit_celltype_comm_glm()` | Sample-level quasibinomial differential communication model |
| `summarize_communication_findings()` | Create publication-ready top findings and QC tables |
| `write_communication_report()` | Write a markdown report from a cell-type communication object |
| `score_communication_specificity()` | Distinguish pair-specific LR axes from broad/identity-associated axes |
| `interpret_celltype_roles()` | Biology-oriented cell-type role interpretation with evidence labels |
| `rank_communication_axes()` | Rank LR axes by strength, specificity, stability, and null evidence |
| `diagnose_celltype_communication()` | Publication-level diagnostic warnings and interpretation notes |
| `diagnose_permutation_resolution()` | Report p-value resolution and degenerate-null status |
| `plot_pathway_dominance()` | Show whether a few pathways dominate total communication strength |
| `plot_celltype_role_dotplot()` | Dotplot of sender/receiver/mediator/influencer scores with evidence alpha |
| `plot_specificity_stability()` | Compare REO-threshold stability and cell-type-pair specificity |
| `all_lr_genes()` | Extract all unique genes from an LR database |
| `summarize_celltype_participation()` | Per-cell-type communicating-cell fraction, pathway/L-R composition, and major-hub label |
| `plot_celltype_participation()` | Fraction of communicating cells per cell type |
| `plot_celltype_pathway_composition()` | Pathway/L-R composition of a cell type's outgoing/incoming communication |
| `plot_celltype_communication_profile()` | Single cell type's communication profile card |
| `differential_celltype_communication()` | Subgroup-resolved multi-sample differential communication |
| `plot_differential_communication_summary()` | Differential L-R pairs per sender -> receiver subgroup |
| `rank_communication_axes()` | Integrated discovery ranking (strength + specificity + stability + permutation + response) |
| `plot_communication_discovery()` | Strength-vs-specificity discovery landscape by evidence tier |

---

## License

MIT © 2026 LogicComm Developer
