# LogicComm 0.7 Rank-Aware REO Design

## Publication-level target

LogicComm 0.7 should move from a "runs on PBMC" package to a package whose
core signal model is defensible in a methods paper. The main target is not only
statistical significance. The package should separate four claims that are often
mixed in single-cell communication tools:

1. The ligand and receptor are active by a within-cell relative-order rule.
2. The ligand and receptor are high-ranking, not merely above a permissive
   binary threshold.
3. The signal is stable across nearby REO thresholds.
4. The signal is cell-type-pair specific, or explicitly labeled as broad.

## RankComp/DRM-inspired upgrade

The previous core REO implementation used one binary rule:

`gene active = expression above the within-cell anchor quantile`.

That is robust, simple, and backwards compatible, but it loses rank margin. A
gene at the 51st percentile and a gene at the 95th percentile both become
active. This is why a PBMC simulation can show identical low-threshold active
frequency in Case and Control while still having clear Case rank strengthening.

The 0.7 upgrade keeps the binary REO/LCS API intact and adds an experimental
continuous evidence layer:

* `calc_REO_rank_score_matrix()` returns within-cell rank percentiles.
* `IdentifyRankLogicConsensus()` reports binary LCS plus rank dominance, rank
  margin, threshold stability, and optional cell-label specificity.
* `CompareRankLogicGroups()` compares rank-aware evidence across biological
  samples when binary active frequency is saturated.

## Evidence fields

`binary_lcs` is the original active-event frequency. `rank_dominance_lcs` is the
average weakest-component rank for ligand-receptor events. `rank_margin_lcs`
measures how far that weakest component sits above the REO threshold.
`threshold_stability_lcs` averages active-event frequency across a grid of rank
thresholds. `specificity_score` is an entropy-based concentration score over
active cell-type pairs when labels and a graph are supplied.

`specificity_weighted_rank_lcs` is the main 0.7 continuous score. It is designed
to keep broad PBMC axes such as MHC or MIF visible, but prevent them from being
mistaken for highly pair-specific mechanisms.

## Evidence tiers

Tier 1 means strong, specific, and stable rank evidence. Tier 2 means active and
rank-supported, but possibly broad. Tier 3 means rank-supported weak candidate
evidence that may deserve follow-up. Tier 4 means weak or threshold-sensitive.

These tiers are interpretation aids, not hypothesis-test p-values. For cohorts,
sample-level comparison should still use biological samples as the unit.

## Remaining publication goals

Before submission, the package should provide:

1. A stable API with no undocumented exported functions.
2. A PBMC single-sample vignette that reports broad immune sanity-check axes
   separately from pair-specific candidates.
3. A PBMC multi-sample vignette where simulated Case/Ctrl effects are recovered
   by both binary and rank-aware evidence.
4. Regression tests for NA labels, graph alignment, self/same-cell-type edges,
   plot exports, clustered heatmaps, permutation diagnostics, and rank-aware
   evidence.
5. R CMD check OK on the release tarball with vignettes ignored for CI and
   vignette rendering verified separately on the analysis server.
6. A methods manuscript that states the distinction between algorithmic
   validation on PBMC simulations and biological discovery on independent data.

