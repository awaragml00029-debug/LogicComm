# LogicComm 0.12.3

## Major change: cell-type co-expression scoring (neighborhood removal, stage 1)

LogicComm is pivoting to a pure cell-type-level cell-cell communication method,
positioned as an interpretable, sample-comparable alternative to CellChat. The
per-cell KNN/SNN *neighborhood* scoring is being removed because, for
dissociated scRNA-seq, that graph lives in expression space (transcriptomic
similarity), not physical space, and therefore cannot license spatial
juxtacrine/paracrine distance claims.

* `summarize_celltype_communication()` now scores communication **only** at the
  cell-type level: for each sender -> receiver pair, the LCS is the fraction of
  the pair's opportunity universe (sender x receiver cell-count product) in which
  the ligand is active in the sender and the receptor is active in the receiver.
  The KNN/SNN neighborhood branch, edge construction, and graph helpers have been
  removed.
* The legacy neighborhood arguments (`knn_mat`, `graph_name`, `mode`,
  `remove_self_edges`, `graph_symmetrize`, `edge_weight_mode`) are accepted via
  `...` for backward compatibility but are **deprecated and ignored with a
  warning**. They will be removed in a later release.
* The former neighborhood/local/distal/range output fields
  (`lcs_neighborhood`, `local_active`, `distal_candidate`, `communication_range`,
  ...) are retained as inert transitional stubs and will be removed once readers
  are migrated.

### Stage 3: per-cell scoring engines de-neighborhooded

The standalone scoring engines turned out to be dual-mode (neighborhood +
global), and the flagship multi-sample pipeline `run_multisample()` is built on
`IdentifyLogicConsensus()`. Rather than hard-deleting and breaking the flagship,
their KNN/neighborhood paths were stripped and their **global co-expression
cores retained**:

* `IdentifyLogicConsensus()` is now global co-expression only (LCS = fraction of
  cells co-expressing the complete ligand and receptor logic). KNN graph
  resolution, symmetrization, edge construction, and weighted-edge scoring
  removed.
* `logic_score_lr()` drops its `mode`/`seurat_obj`/`knn_mat`/graph arguments and
  scores globally (gate-aware path unchanged).
* `run_multisample()` computes per-sample LCS by global co-expression; it no
  longer passes a per-sample KNN graph.
* Legacy neighborhood arguments remain accepted via `...` and ignored with a
  warning.

### Stage 5: trustworthy statistics for selecting real communications

* **Axis-level permutation null.** `permute_celltype_communication()` now
  defaults to an axis-level null (`metric = "lcs"`): every sender -> receiver ->
  L-R axis gets its **own** empirical p-value and BH FDR, instead of one
  cell-type-pair-level p-value being broadcast onto all of its L-R pairs. The
  output gains an `lr_pair` column; `rank_communication_axes()` joins the
  permutation evidence per axis and exposes `permutation_fdr`. Pair-level nulls
  (`metric = "sum_lcs"`) remain available for backward compatibility. The
  cell-type co-expression pivot makes per-axis nulls cheap (no graph to rescan).
* **Per-cell bootstrap.** `bootstrap_celltype_communication()` now resamples on
  the independent unit (cells) rather than the sender x receiver cell-count
  *product*. Co-expression LCS = (ligand-active fraction of sender cells) x
  (receptor-active fraction of receiver cells), so each fraction is bootstrapped
  as a binomial proportion over its own cell count. The previous approach used
  `n_edges` (the product) as the trial size and produced anticonservative,
  far-too-narrow intervals.

### Standard discovery pipeline (footgun-free)

* New `discover_celltype_communication()` wires the whole workflow in one call
  with the correct defaults: cell-type co-expression scoring, specificity +
  proliferation-confound annotation, an **axis-level** permutation null
  (`metric = "lcs"`, so it is impossible to accidentally request a pair-level
  null via `metric = "sum_lcs"`), evidence ranking, a confound-filtered
  discovery view, and an FDR-passing `shortlist`. Pass `expr =` the counts matrix
  to enable the proliferation/breadth filter (cycling clusters otherwise dominate
  the ranking). A copy-paste template lives in
  `inst/workflows/standard_discovery.R`.

### Stage 4: REO rank-weighted LCS (opt-in)

* `summarize_celltype_communication(lcs_weighting = "rank")` scores by REO
  intensity instead of the binary co-expression fraction: the prevalence-weighted
  within-cell rank of the ligand in the sender type (fraction expressing times
  mean rank among expressers) times the same for the receptor. This restores
  dynamic range that binarizing discards (a ligand at the 99th within-cell
  percentile separates from one at the 51st) while keeping the prevalence signal
  that separates specific from ubiquitous axes, addressing the compressed,
  near-floor LCS values seen on real data. On the benchmark it edges out the
  binary score (AUPRC/sens@k); averaging rank among expressers only -- an earlier
  formulation the benchmark flagged -- discards prevalence and hurts precision.
  Requires the rank matrix from `calc_REO_matrix(..., return_rank = TRUE)`. The
  default remains `"binary"`, so existing results are unchanged.
* The `active` call still uses the binary co-expression fraction, so the active
  axis set is identical under either weighting -- only the reported `lcs` differs.
* `permute_celltype_communication()` inherits `lcs_weighting` from the scored
  object and carries the rank matrix through, so a rank-weighted observed score is
  tested against a rank-weighted null (coherent significance). It also flows
  through `discover_celltype_communication(..., lcs_weighting = "rank")` when the
  input is built with `return_rank = TRUE`.

### Stage 6 (started): benchmark harness vs. baselines

* `inst/benchmark/benchmark_vs_baselines.R`: a sandbox-runnable harness that
  simulates data with a known ground truth and the confounds that break naive
  scores -- ubiquitous/housekeeping pairs, a broad moderate-abundance pair, a
  CYCLING / transcriptional-breadth hub cell type that co-expresses many L-R
  genes, uneven cell-type sizes incl. a rare type, and true axes at three signal
  strengths. It scores every sender -> receiver -> L-R axis with LogicComm
  (binary, rank, and +proliferation-filter) and with re-implemented baselines
  (CellPhoneDB/CellChat-style mean-expression product; naive co-detection), and
  reports AUROC / AUPRC / sensitivity-at-k (random tie-breaking, averaged over
  sims). Real, runnable `score_cellchat()` and `score_liana()` adapters (and a
  note on CellPhoneDB via LIANA) are included for use where those packages exist.
* Result (mean over sims; 1188 axes, 12 true positives): LogicComm with the
  proliferation filter reaches AUROC 0.995 / AUPRC 0.80, LogicComm (rank) 0.97 /
  0.63 and (binary) 0.99 / 0.58, while the mean-expression-product and naive
  baselines collapse to AUPRC ~0.04 -- ubiquitous, abundance and cycling-breadth
  confounds saturate their scores. The proliferation filter is the single biggest
  contributor to precision.

Still pending (non-blocking; the package is functional and green without them):
removal of the inert transitional stub columns (`lcs_neighborhood`,
`communication_range`, `local_active`, `distal_candidate`, ...) that the
cell-type path still emits as dead weight -- a cross-cutting refactor of the
downstream readers; and a possible re-architecture of `IdentifyRankLogicConsensus()`
to cell-type level (its graph-free mode is autocrine-only, so it keeps an optional
neighborhood mode for now). The gate-aware consensus and the spatial module retain
a cell-level graph by design. Also pending: real-package (CellChat/LIANA) and
real-data benchmark runs.

# LogicComm 0.11.1

## Documentation

* Seurat tutorial: the permutation-null step now demonstrates `n_cores`
  (fork-parallel) and explains the `1 / (n_perm + 1)` empirical p-value floor and
  how to pick `n_perm`; the discovery-view step notes that
  `communication_discovery_view()` accepts a `ct_comm` object directly.

# LogicComm 0.11.0

## Faster permutation null + discovery-view ergonomics

* `permute_celltype_communication()` gains an `n_cores` argument: the permutation
  loop now runs on forked workers (Unix/macOS) for a near-linear speedup, so a
  publication-grade `n_perm` (e.g. 200-1000) is affordable. Results are
  reproducible across core counts when `seed` is set; the serial path
  (`n_cores = 1`, default) is unchanged. Note: with too few permutations the
  empirical p-value floor is `1 / (n_perm + 1)`, which caps `null_support` and
  prevents any Tier 1 -- use enough permutations to resolve p < 0.05.
* `permute_celltype_communication()` now fails early with a clear message when
  `knn_mat` is missing for neighborhood mode (previously a cryptic
  `.extract_knn` error).
* `communication_discovery_view()` now accepts a `LogicCommCellTypeComm` object
  directly (it ranks the axes on the fly), in addition to a ranked data.frame
  from `rank_communication_axes()`.

# LogicComm 0.10.4

## Documentation

* Expanded the Seurat tutorial's publication-figures section (9.12) with the
  brand scale reference and the per-figure publication controls added across
  0.10.x (volcano `fdr_cutoff` / `lfc_threshold`, network `layout`, bubble
  `top_n_pathways`, discovery/roles `subtitle`, heatmap clustering).
* README documents the figure system and `save_logiccomm_figure()`; refreshed
  stale install-version strings in the README and intro vignette.

# LogicComm 0.10.3

## Figure polish (round 4 -- remaining figures)

* `plot_celltype_roles()` fills points by the brand role palette, adds a
  diagonal-orientation subtitle and a `subtitle` argument, and uses
  white-stroked points (it previously used default ggplot hues).
* `plot_role_confidence()`, `plot_celltype_pathway_composition()`,
  `plot_celltype_communication_profile()`, and `plot_communication_dynamics()`
  now use the brand role / qualitative palettes instead of default ggplot hues;
  `plot_celltype_participation()` uses the brand sender/receiver/communicating
  colours.
* `plot_communication_qc()` uses a white-stroked, brand-filled point style.

# LogicComm 0.10.2

## Figure polish (round 3)

* `plot_celltype_glm_volcano()` now matches `plot_differential_celltype_volcano()`:
  points coloured by significance (up / down / n.s.), FDR and effect-size guide
  lines, readable axis labels, and `fdr_cutoff` / `effect_threshold` arguments
  (it previously coloured by sender type with default hues and drew no thresholds).
* `plot_lr_bubble_advanced()` caps the pathway facets via a new `top_n_pathways`
  argument (default 12, ranked by total LCS) so the figure stays legible instead
  of fragmenting into many tiny panels.

# LogicComm 0.10.1

## Figure polish (round 2)

* `plot_celltype_network()` now uses the same force-directed igraph layout as
  `plot_celltype_network_publication()` (with a circle fallback) and gains a
  `layout` argument; edge pathways use the brand qualitative palette (the basic
  network previously fell back to default ggplot hues).
* `plot_celltype_role_radar()` uses the brand role palette
  (Sender / Receiver / Mediator / Influencer) instead of default hues.
* `plot_pathway_heatmap()` clusters rows and columns by default (with a safe
  fallback to the original order) so related pathways and cell-type pairs form
  readable blocks.
* `plot_lcs_heatmap()` uses brand colours for the heatmap ramp and the
  Case/Ctrl annotation.

# LogicComm 0.10.0

## Publication-grade figure system

All plotting functions now share one visual language designed for journal
figures: colorblind-safe, perceptually ordered, and consistent across panels.

* New exported visual system: `theme_logiccomm()` (export-ready typography, with
  `grid`/`legend` controls), the LogicComm brand palette (`logiccomm_brand`), and
  brand colour/fill scales -- `scale_color/fill_logiccomm_d` (qualitative),
  `_logiccomm_c` (sequential mako), `_logiccomm_diverging` (teal-amber), and an
  ordered `scale_color/fill_tier` for `evidence_tier`. Every fixed-meaning palette
  avoids red/green confusion and differs in luminance for greyscale printing.
  (`theme_logiccomm()` was previously defined but never actually exported.)
* `save_logiccomm_figure()` exports at journal column widths (single 89 mm,
  double 183 mm) to vector (PDF/SVG) or high-DPI raster.
* Every `plot_*` function was converted to the shared theme and brand scales,
  replacing ad-hoc per-figure gradients, ColorBrewer Set1 (red/green), and the
  `turbo` palette.
* `plot_communication_discovery()` now colours by an *ordered* tier scale
  (strong = dark, broad/non-specific = amber, weak = grey), adds a strong+specific
  quadrant cue, and gains a `subtitle` argument.
* `plot_differential_celltype_volcano()` colours points by significance
  (up / down / n.s.), draws the FDR and effect-size guide lines, and gains
  `fdr_cutoff`, `lfc_threshold`, and `x_lab` arguments.
* `plot_celltype_network_publication()` uses a real force-directed layout
  (igraph, with a circle fallback) instead of a fixed circle, and gains a
  `layout` argument; edge pathways use the brand qualitative palette.
* New Suggests: `igraph` (graph layouts) and `patchwork` (multi-panel figures).

# LogicComm 0.9.4

## One-call confound-filtered discovery view

`rank_communication_axes()` scales strength/specificity across the whole active
set, so after the 0.9.3 broad-axis demotion the genuine but lower-LCS candidates
still inherited Tier 3/4 -- the high-LCS broad/cycling axes set the ceiling.

* New `communication_discovery_view()` filters a ranked table down to the
  confound-free candidates -- dropping broad/ubiquitous axes (`broad_axis_flag`)
  and proliferation-hub axes (`proliferation_confound_flag`), and optionally
  identity-associated axes (`drop_identity`) -- and, with `rescale = TRUE`
  (default), recomputes `strength_score`, `specificity_score`,
  `discovery_score`, and `evidence_tier` *within the kept set* using the same
  weights and tier rules as the ranker. Surviving candidates are therefore tiered
  on their own merits instead of against the noise that was removed.
* The discovery scoring/tiering/demotion logic is refactored into a single shared
  internal core used by both `rank_communication_axes()` and
  `communication_discovery_view()`, so the two always agree.
  `rank_communication_axes()` output is unchanged.
* Added regression tests for the discovery view.

# LogicComm 0.9.3

## Broad/ubiquitous axes are down-ranked, not just annotated

`score_communication_specificity()` already detected broad axes (MHC-I -> CD8A,
CD99 - CD99, and other identity/ubiquitous interactions), but the flag was
annotation only: a strong, null-supported broad axis still reached Tier 1/2 in
`rank_communication_axes()` because the tier rule never consulted it. This let
non-specific interactions dominate the discovery ranking -- including
identity-associated broad axes where a marker survives the `min_expr_frac` gate
(e.g. CD8A logic-active in ~26% of an activated Treg cluster from
ambient/doublet contamination, well above the 10% gate).

* `rank_communication_axes()` gains `demote_broad` (default `TRUE`) and
  `broad_penalty` (default `0.5`). Axes with `ubiquitous_interaction_flag == TRUE`
  have their `discovery_score` multiplied by `broad_penalty` and their
  `evidence_tier` capped at Tier 3 ("broad / non-specific"), so a broad axis can
  no longer outrank a genuinely cell-type-pair-specific candidate. This is a
  ranking-layer change only -- the `active` flag and LCS are untouched -- and
  `demote_broad = FALSE` restores the previous behaviour. A `broad_axis_flag`
  column is added to the ranked output.

## Proliferation / transcriptional-breadth confound diagnostic

Cycling cells have unusually broad transcriptomes, so they call many genes
logic-active and form active L-R edges with almost every partner: a proliferation
"hub" (e.g. a Cycling_Lymphocytes cluster) can appear in the majority of
top-ranked axes without representing specific biology.

* New `diagnose_proliferation_confound(ct_comm, expr)` scores each cell's S/G2M
  signatures (Tirosh et al. 2016 sets by default) and transcriptional breadth,
  aggregates them per cell type, flags proliferation/breadth-outlier cell types,
  and annotates `lr_table` with `sender_proliferation_hub`,
  `receiver_proliferation_hub`, `proliferation_confound_flag`, and
  `proliferation_hub_role`. It is a diagnostic only: it does not change the
  `active` flag, LCS, or `discovery_score`. Use the flags to filter or
  down-weight axes, or pass the object straight to `rank_communication_axes()`
  (the flags ride along). A per-cell-type `proliferation_summary` is stored on
  the object.

* Added regression tests for the broad-axis demotion and the proliferation
  diagnostic.

# LogicComm 0.9.2

## Cell-type expressing-fraction gate removes hub-inflated phantom axes

Neighborhood LCS is an *edge* fraction, so a handful of high-degree (hub)
receiver cells carrying residual ambient receptor signal could make an entire
sender -> receiver L-R axis look active even when the receptor is expressed in a
negligible fraction of that cell type. This is why, even after the 0.9.1 ambient
guard cut `CD8A` activity in Treg to ~0.4% of cells, a `... -> Treg` `CD8A` axis
could still surface with a high `lcs_neighborhood` (~0.5) and a Tier 2
`discovery_score` -- the few contaminated Treg cells happened to be KNN hubs.

* `summarize_celltype_communication()` gains `min_expr_frac` (default `0.1`): an
  L-R axis is called active only if the ligand is logic-active in at least this
  fraction of **sender** cells *and* the receptor in at least this fraction of
  **receiver** cells. This is a graph-independent, cell-type-free gate (the same
  expressing-fraction guard used by CellChat/CellPhoneDB), so it is immune to the
  hub effect. Set `min_expr_frac = 0` to restore the previous behaviour.
* Because the gate flips the `active` flag, the phantom axes now disappear from
  every downstream consumer that filters on active rows -- `rank_communication_axes()`
  (no `discovery_score`/`evidence_tier`), `score_communication_specificity()`,
  the differential workflow, and the bubble/heatmap plots -- rather than merely
  having a lower LCS. In the hub reproduction the spurious `CD8T -> Treg` and
  `Treg -> Treg` `CD8A` axes drop out of the discovery ranking entirely while the
  genuine CD8 T signal (receptor frac ~0.45) is preserved.
* `permute_celltype_communication()` inherits `min_expr_frac` from the observed
  object so the null permutations score with the same definition.
* Added regression tests reproducing the KNN-hub artifact end to end.

## Actionable errors for cell-type-filtered plots

`plot_lr_bubble_advanced()`, `plot_lr_bubble_by_celltype()`, and
`plot_lr_activity_balance()` previously failed with a blunt
`"No active events to plot."` when a `senders`/`receivers` name was mistyped or
carried no signal. They now distinguish the two cases: a name that matches no
cell type lists the available types, while a valid type with no active events
explains which gates removed them and what to adjust.

# LogicComm 0.9.1

## Make the ambient guard effective for minority-expressed genes

The 0.9.0 `gene_background = "quantile"` guard took each gene's background
quantile over *all* cells. For a marker detected in fewer than half of cells
(the common case, e.g. `CD8A` in a mixed immune dataset) that quantile is 0, so
the guard had no effect and ambient receptors such as `CD8A` in Treg still
appeared.

* The background is now taken over the cells that **detect** the gene (nonzero
  counts), so it stays effective regardless of detection rate. In a low-detection
  simulation the spurious sender -> Treg `CD8A` axis drops from LCS ~0.16 to 0
  while genuine CD8 T expression is preserved; in high-detection data behaviour
  is essentially unchanged.
* Added a regression test for the minority-detection case.

# LogicComm 0.9.0

## Cell-type-free ambient guard for REO activity calls

The within-cell REO anchor (a gene is active if it exceeds that cell's own
median gene) cannot tell genuine expression from ambient RNA: a marker such as
`CD8A` that contaminates non-CD8 cells with a few ambient counts still clears the
cell's low median and is called active. This produced biologically impossible
results such as a `HLA-B -> CD8A` axis with a Treg receiver.

* `calc_REO_matrix()` (and `logic_prepare()`) gain `gene_background` and
  `gene_background_quantile`. With `gene_background = "quantile"`, a gene is
  active in a cell only if it also exceeds its own across-cell background (the
  given quantile of that gene over all cells), which removes ambient activity
  **without requiring any cell-type annotation**. The default is `"none"` for
  backward compatibility; the Seurat demo and README now recommend the guard.
* In a controlled ambient simulation this drops the spurious sender -> Treg CD8A
  axis from LCS ~0.79 to ~0.02 while preserving genuine CD8 T self-signal.
* Documentation recommends upstream decontamination (SoupX / DecontX /
  CellBender) for rigorous ambient removal, and explains that REO is a
  within-cell relative measure.
* Added regression tests for the ambient guard.

# LogicComm 0.8.1

## Publication-figure readability pass

A visualization audit fixed labels and palettes that became unreadable on real
datasets with long cell-type names and many pathways:

* `plot_specificity_stability()` no longer piles overlapping full
  `sender|receiver|L-R` keys at the top of the panel. It labels with short L-R
  names by default, gains a `label_field` argument (`"lr_pair"`, `"feature"`,
  `"none"`), breaks ties in the crowded high-stability/high-specificity corner by
  total LCS, and uses a finite `max.overlaps`.
* `plot_celltype_network_publication()` switched the pathway colour scale from
  `Set2` (only 8 colours, which failed when there were more pathways) to a
  viridis scale, and shortens node labels.
* `plot_communication_discovery()`, `plot_differential_celltype_volcano()`,
  `plot_celltype_glm_volcano()`, `plot_celltype_network()`, and
  `plot_celltype_roles()` now use short/compact point and node labels and a
  finite `max.overlaps` so labels repel legibly instead of stacking.
* `plot_celltype_glm_volcano()` guards the optional `sender_type` colour column.
* Added internal `.short_label()` / `.compact_feature_label()` helpers (ASCII
  source, `…` / `→` escapes) and regression tests.

# LogicComm 0.8.0

## A clearer discovery workflow: who communicates, what changes, and what to prioritize

This release adds three layers that turn LogicComm's computed evidence into
results you can read off directly. The REO/LCS core is unchanged.

### Per-cell-type communication participation (which subgroup, how broadly, through what)

* `summarize_celltype_participation()`: for every cell type, the fraction of
  cells that actually express the ligand/receptor machinery
  (`frac_communicating` / `frac_sender_active` / `frac_receiver_active`, using
  graph-independent within-cell REO activity), the pathway and ligand-receptor
  composition of its outgoing and incoming communication, and a transparent
  major-hub importance label.
* `plot_celltype_participation()`, `plot_celltype_pathway_composition()`, and
  `plot_celltype_communication_profile()` visualize the communicating-cell
  fraction, the pathway/L-R composition, and a single subgroup's profile card.

### Subgroup-resolved multi-sample differential communication (which cell-type pair carries the change)

* `differential_celltype_communication()`: compares cell-type-resolved
  communication between sample groups and splits results into explicit
  `sender_type` / `receiver_type` / `lr_pair` / `pathway` / `direction` columns
  plus a per sender -> receiver subgroup summary, so the differential pairs and
  the subgroups that carry them are directly visible.
* `plot_differential_communication_summary()` shows differential L-R counts per
  subgroup; `plot_differential_celltype_heatmap()` now plots
  sender -> receiver x L-R pair (it previously dropped the receiver).
* The output documents the small-sample FDR caveat: with few replicates per
  group Fisher FDR saturates near 1, so rank by effect size and corroborate with
  `fit_celltype_comm_glm()` and more replicates.

### Integrated discovery ranking (a single prioritized candidate list)

* `rank_communication_axes()` now integrates communication strength,
  cell-type-pair specificity, REO threshold stability (`sens`), cell-label
  permutation support (`null_pair`), and receiver response into a single
  `discovery_score` with an `evidence_tier` (previously it ignored `null_pair`
  and `sens`).
* `plot_communication_discovery()` draws the strength-vs-specificity discovery
  landscape coloured by evidence tier.

# LogicComm 0.7.2

## Correctness audit, documentation pipeline, and check hygiene

Correctness fixes (with new regression tests in
`tests/testthat/test-correctness-audit.R`):

* Fixed an inverted Mediator role: `summarize_celltype_communication()` role
  centrality passed communication strength to `igraph::betweenness()`, which
  treats edge weights as distances, so betweenness now uses inverse strength.
* Non-finite betweenness/centrality (for example a single cell type, where
  normalized betweenness is `NaN`) is sanitized to `0` instead of propagating to
  `NA` roles, and self-loops are excluded from the centrality graph.
* The Influencer / information-flow score now uses PageRank instead of
  `eigen_centrality(directed = TRUE)`, which returned all zeros (or collapsed
  onto sink nodes) on the acyclic / weakly connected communication graphs that
  are common in practice. PageRank is the robust analogue of incoming
  eigenvector centrality and keeps the same `influencer_role_score` /
  `information_score` columns.
* Resolved a duplicate `lr_pair` in `lr_pairs_human` (`POSTN_ITGAV` is now
  `POSTN_ITGAV_ITGB1` and `POSTN_ITGAV_ITGB3`) that caused silent row-name
  collisions in per-cell scoring.
* `rank_comm_cells()` now aligns unnamed `cell_labels` to the ranked cells
  instead of the original order.
* `fit_celltype_comm_glm()` no longer clamps distal-candidate features to 100%
  active; binomial successes are taken from the same opportunity universe as the
  edge-count denominator.
* `run_multisample()` errors when `group_info` does not label every sample, and
  minor `na.rm`/non-finite hardening was added to modulator and cluster helpers.

Documentation and packaging:

* All `man/*.Rd` pages are regenerated from roxygen source so documentation is
  idempotent; previously 26 hand-written pages had drifted from the source
  comments. Fixed malformed `\code{}` tags and de-duplicated the internal
  `%||%` operator.
* Made the intro-vignette heatmap chunk conditional on the suggested `pheatmap`
  package so vignettes build without it, dropped unused `Suggests`
  (`ggnetwork`, `viridis`), and added a `.gitignore` for build artifacts.

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
