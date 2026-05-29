# LogicComm 0.7.1

## Cell-type communication consistency and visualization audit

* `permute_celltype_communication()` now inherits graph and scoring parameters from `ct_comm$params` so observed summaries and permutation nulls use matching definitions by default.
* Added warnings when permutation calls intentionally override observed scoring parameters with different values.
* Added local-support annotations to cell-type pair and pathway summaries, including `communication_support_label` and `local_support_fraction_active`, to distinguish local graph-supported, mixed, and global-only distal candidates.
* Updated cell-type network and LR bubble plots to visually fade global-only distal candidates and document their interpretation in captions.
* Split report output into primary local/mixed cell-type pairs and global-only distal candidate pairs so candidate distal signals are not overinterpreted as direct neighborhood interactions.
* Updated the Seurat demo and regression tests to audit observed/null scoring consistency and global-only distal candidate reporting.

# LogicComm 0.7

## Rank-aware REO evidence for publication-grade interpretation

* Added `calc_REO_rank_score_matrix()` to expose within-cell rank percentiles for ligand-receptor genes without replacing the original binary REO workflow.
* Added `IdentifyRankLogicConsensus()` with binary LCS, rank dominance, rank margin, threshold stability, cell-label specificity, and evidence-tier annotations.
* Added `score_lr_rank_activity()` for continuous per-cell sender and receiver LR rank activity.
* Added `CompareRankLogicGroups()` for sample-level Case/Ctrl comparisons when binary active frequency is saturated but rank evidence differs.
* Added a 0.7 design note at `inst/design/LogicComm_0.7_rank_aware_reo.md` documenting publication-level package targets, evidence tiers, and RankComp/DRM-inspired rationale.
* Added regression tests showing that rank-aware evidence recovers simulated Case/Ctrl strength differences even when binary LCS is identical.
* Extended the PBMC multi-sample demo with 0.7 rank-aware Case/Ctrl evidence output.

# LogicComm 0.6.7

## PBMC multi-sample validation and publication cleanup

* Added an executable PBMC3K four-sample validation demo at `inst/examples/pbmc_multisample_demo.R`. The script splits PBMC3K into `pbmc1`, `pbmc2`, `pbmc3`, and `pbmc4`, applies a controlled Case perturbation to `pbmc1` and `pbmc2`, and verifies that LogicComm recovers the expected sample-level MIF, galectin, and MHC-I LR signals.
* Added `vignettes/PBMC_multisample_demo.Rmd` with detailed PBMC1-4 Case/Ctrl simulation design, expected output, plotting code, and interpretation cautions for low sample counts.
* Added a regression test for four named PBMC-style samples with a known ligand-receptor perturbation.
* Quieted known role-summary warnings from directed acyclic igraph eigen-centrality calls and all-NA role-score rescaling.
* Updated package citation and installation examples to the 0.6.7 source tarball.

# LogicComm 0.6.1

## Specificity, diagnostics, and safer biological interpretation

* Added `score_communication_specificity()` to distinguish pair-specific LR axes from broad or identity-associated programs such as MHC/HLA-dominated PBMC communication. The function annotates `lr_table`, `pair_summary`, and adds `specificity_summary`.
* Split cell-type role interpretation into separate concepts: `role_separation_label`, `communication_evidence_label`, `role_reliability_label`, and `dominant_role_strict`. This avoids interpreting a high topological role margin as high biological evidence when total communication support is modest.
* Added biology-oriented role text columns: `role_biological_interpretation`, `role_interpretation_caution`, and `role_recommended_followup`.
* Improved `permute_celltype_communication()` with optional adaptive top-candidate refinement, `degenerate_null`, `degenerate_positive_null`, infinite z-scores for observed values above zero-variance nulls, p-value resolution fields, and biological interpretation text for each sender-receiver null result.
* Added `diagnose_permutation_resolution()` and `diagnose_celltype_communication()` for publication-oriented QC warnings covering pathway dominance, low-evidence role labels, small cell types, broad/identity-associated LR axes, insufficient permutations, and degenerate nulls.
* Added `interpret_celltype_roles()` and `rank_communication_axes()` to produce interpretable role summaries and evidence-aware candidate LR rankings.
* Added publication-oriented plots: `plot_pathway_dominance()`, `plot_celltype_role_dotplot()`, `plot_specificity_stability()`, `plot_celltype_network_publication()`, and `plot_pathway_network()`.
* Updated README and Seurat demo with PBMC-specific interpretation guidance: stable-but-broad MHC/MIF axes should be separated from pair-specific mechanistic candidates, and zero-variance nulls should be reported as structural-null cases rather than ordinary z-score estimates.

# LogicComm 0.6.0

## Publication-level analysis extensions

* Added spatial communication support: `extract_spatial_coordinates()`,
  `build_spatial_graph()`, `summarize_spatial_communication()`, and
  `plot_spatial_logic()`. Spatial graphs can be used directly as LogicComm
  neighborhood graphs.
* Added pseudotime/time-bin dynamics: `summarize_communication_dynamics()` and
  `plot_communication_dynamics()`.
* Added sample-level quasibinomial differential communication modeling with
  `fit_celltype_comm_glm()` and `plot_celltype_glm_volcano()`. This keeps the
  biological sample as the statistical unit while modeling active-edge counts.
* Added publication-oriented reporting helpers: `summarize_communication_findings()`,
  `write_communication_report()`, `plot_communication_qc()`, and
  `plot_role_confidence()`.
* Cell-type summaries now include `celltype_sizes`, `sender_n_cells`,
  `receiver_n_cells`, and the total cell count in object parameters for clearer
  edge-opportunity interpretation.
* README and Seurat vignette were expanded with story-driven interpretation,
  spatial analysis, pseudotime dynamics, GLM-based multi-sample modeling,
  publication figure recipes, and common interpretation pitfalls.

# LogicComm 0.5.0

## Publication-oriented cell-type communication interpretation

- Added weighted KNN/SNN scoring via `edge_weight_mode = "weighted"` and graph
  symmetrization options (`graph_symmetrize = "none"`, `"or"`, or `"max"`) to
  `IdentifyLogicConsensus()`, `IdentifyLogicGateConsensus()`,
  `summarize_celltype_communication()`, and `run_multisample()`.
- Expanded `role_summary` with clearer count fields: active L-R event counts,
  unique L-R counts, source/target cell-type degree counts, and active edge
  support. Backward-compatible aliases such as `outgoing_count` and
  `out_degree_n` are retained.
- Added low-communication role filtering, `dominant_role_unfiltered`,
  `secondary_role`, `role_confidence`, and `role_confidence_label` so weak or
  mixed cell types are not overinterpreted as definitive mediators/influencers.
- Revised autocrine fractions to avoid double-counting self-communication in
  total hub strength and added outgoing/incoming autocrine fractions.
- Added receiver-response scoring with `score_receiver_response()` and
  `add_receiver_response_score()` for downstream response gene evidence.
- Added uncertainty and robustness utilities: `bootstrap_celltype_communication()`,
  `permute_celltype_communication()`, and `sensitivity_REO_threshold()`.
- Added interpretability visualizations: role heatmap/radar, ligand/receptor
  activity balance, single-event evidence plots, and differential cell-type
  heatmap/volcano plots.
- Updated README and Seurat demo with weighted SNN analysis, role-count
  interpretation, receiver-response examples, bootstrap/null diagnostics, and
  differential cell-type communication visualizations.

# LogicComm 0.4.2

## Cell-type-resolved biology and role analysis

- Added `summarize_celltype_communication()` to compute sender-cell-type to
  receiver-cell-type LCS tables, interaction counts, communication strength,
  pathway summaries, and cell-type signaling roles.
- Added network role metrics inspired by signaling-role centrality analysis:
  `outdegree_score` for sender roles, `indegree_score` for receiver roles,
  `betweenness_score` for mediator roles, and `information_score` for
  influencer roles.
- Added `celltype_comm_to_lcs()` so cell-type-level LCS features can be compared
  across biological samples with `CompareLogicGroups()`.
- Added visualization helpers: `plot_celltype_heatmap()`,
  `plot_celltype_network()`, `plot_celltype_roles()`,
  `plot_lr_bubble_by_celltype()`, `plot_pathway_heatmap()`, and
  `explain_celltype_interaction()`.
- `run_multisample()` now exposes `graph_name`, `layer`, and
  `remove_self_edges`, making Seurat multi-sample workflows easier and more
  reproducible.
- Updated README and the Seurat tutorial with a biology-oriented visualization
  workflow: interaction quantity, communication strength, pathway drivers,
  sender/receiver/mediator/influencer role positioning, and sample-level
  cell-type communication comparison.

# LogicComm 0.4.1

## Core fixes

- `calc_REO_matrix()` now computes REO anchors from all genes by default while
  using `lr_genes` only to choose retained output rows. This keeps binary logic
  states stable when switching or subsetting ligand-receptor databases.
- Added `anchor_genes` for deliberate sensitivity analyses with a custom anchor
  universe.
- Fixed `positive_gate = "any"` in `IdentifyLogicGateConsensus()` so any active
  positive modulator is sufficient.
- Negative modulator active-block summaries now use any-active logic for
  multi-gene modulator sets.
- `CompareLogicGroups()` now treats missing LCS values as unavailable rather than
  negative calls and reports per-pair non-missing group counts.
- `CompareLogicGroups()` accepts a custom `lr_db` so metadata are preserved for
  non-built-in LR databases.
- `filter_lcs()` no longer copies incompatible data-frame attributes.
- KNN scoring now removes self-loop edges by default and accepts dense adjacency
  matrices after validation.

## Package quality

- Added `tests/testthat.R` so testthat tests run under `R CMD check`.
- Added regression tests for anchor scope, positive modulator `any` gates,
  missing-LCS handling, KNN self-loop removal, and `filter_lcs()` attributes.
- Added manual Rd help files for exported functions and datasets.
- Updated tutorials and README for the corrected anchor semantics and v0.4.1.
- Fixed MIT license file format for `MIT + file LICENSE`.
