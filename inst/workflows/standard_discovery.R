# inst/workflows/standard_discovery.R
#
# Standard LogicComm cell-type communication discovery flow.
#
# This is the recommended, footgun-free pipeline. The single call to
# discover_celltype_communication() folds in everything that is easy to miss
# when calling the steps by hand:
#   * specificity scoring (down-ranks ubiquitous/housekeeping axes),
#   * the proliferation / transcriptional-breadth confound (cycling clusters
#     otherwise dominate the ranking) -- requires passing `expr`,
#   * an AXIS-LEVEL permutation null (metric = "lcs") so every L-R axis gets its
#     own empirical p and BH FDR. (Passing metric = "sum_lcs" by hand silently
#     gives a cell-type-pair p broadcast across all of a pair's L-R axes.)
#
# Inputs you provide:
#   counts : gene x cell expression matrix (raw counts), with gene symbols as
#            rownames and cell IDs as colnames.
#   labels : named character vector of cell-type labels (names = cell IDs).
#   lr_db  : ligand-receptor database (defaults to lr_pairs_human).

library(LogicComm)

## ---- 1. inputs (replace with your data) -------------------------------------
# counts <- <gene x cell matrix>
# labels <- <named cell-type vector>
data(lr_pairs_human)
lr_db <- lr_pairs_human

## ---- 2. REO binarization ----------------------------------------------------
lr_genes <- all_lr_genes(lr_db)
reo <- calc_REO_matrix(counts, lr_genes = lr_genes)   # add anchor_genes/rank_threshold as needed

## ---- 3. one-call discovery --------------------------------------------------
# Pass `expr = counts` so cycling clusters (e.g. Cycling_Lymphocytes) are
# diagnosed and dropped. Co-expression nulls are cheap, so use a publication-grade
# n_perm with n_cores to parallelize on Unix/macOS.
res <- discover_celltype_communication(
  reo_mat     = reo,
  cell_labels = labels,
  lr_db       = lr_db,
  expr        = counts,     # enables the proliferation-confound filter
  n_perm      = 1000,
  n_cores     = 8,
  seed        = 1,
  fdr_cutoff  = 0.1
)

## ---- 4. outputs -------------------------------------------------------------
# res$ranked    : every active axis, with per-axis permutation_empirical_p /
#                 permutation_fdr, discovery_score, evidence_tier.
# res$view      : confound-filtered discovery view (broad + proliferation dropped).
# res$shortlist : axes passing permutation_fdr <= fdr_cutoff, ordered by FDR.
print(utils::head(res$shortlist, 30))

# Sanity checks (what the v0.12 statistics fixes give you):
stopifnot("lr_pair" %in% names(res$null))                       # axis-level null
length(unique(res$ranked$permutation_empirical_p))              # >> 1 (per-axis, not broadcast)

# Optional figure:
# library(ggplot2)
# p <- plot_communication_discovery(res$ranked)
# save_logiccomm_figure(p, "discovery.pdf", width = "double")
