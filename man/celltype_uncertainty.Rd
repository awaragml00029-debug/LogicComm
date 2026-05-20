\name{celltype_uncertainty}
\alias{bootstrap_celltype_communication}
\alias{permute_celltype_communication}
\alias{sensitivity_REO_threshold}
\title{Uncertainty and sensitivity analysis for LogicComm cell-type communication}
\description{
Bootstrap, permutation-null, and REO-threshold sensitivity utilities for
cell-type-resolved LogicComm results.
}
\usage{
bootstrap_celltype_communication(ct_comm, n_boot = 100,
  level = c("celltype_lr", "celltype_pair"), seed = NULL)

permute_celltype_communication(ct_comm = NULL, reo_mat, cell_labels = NULL,
  knn_mat = NULL, lr_db = NULL, n_perm = 100, metric = "sum_lcs",
  seed = NULL, adaptive = FALSE, adaptive_top_n = 25,
  adaptive_n_perm = 1000, min_n_perm_publication = 1000,
  verbose = TRUE, ...)

sensitivity_REO_threshold(expr_mat,
  rank_threshold_grid = c(0.4, 0.5, 0.6, 0.7),
  lr_db = lr_pairs_human, cell_labels = NULL, knn_mat = NULL,
  lcs_threshold = 0.01, ...)
}
\arguments{
\item{ct_comm}{A \code{LogicCommCellTypeComm} object.}
\item{n_boot}{Number of bootstrap replicates.}
\item{level}{Whether to summarize bootstrap uncertainty at celltype-LR or celltype-pair level.}
\item{seed}{Optional random seed.}
\item{reo_mat}{Binary REO matrix or \code{LogicCommREOResult}.}
\item{cell_labels}{Named cell-type labels.}
\item{knn_mat}{Optional KNN/SNN graph.}
\item{lr_db}{Ligand-receptor database.}
\item{n_perm}{Number of label permutations.}
\item{metric}{Pair-summary metric tested against the permutation null.}
\item{adaptive}{Whether to refine top preliminary candidates with additional permutations.}
\item{adaptive_top_n}{Number of preliminary candidates to refine when adaptive is TRUE.}
\item{adaptive_n_perm}{Target number of permutations for adaptive candidates.}
\item{min_n_perm_publication}{Recommended permutation count for publication-level empirical p-value resolution.}
\item{verbose}{Whether to print progress messages.}
\item{...}{Additional arguments passed to \code{summarize_celltype_communication()} or \code{calc_REO_matrix()}.}
\item{expr_mat}{Expression matrix, Seurat object, or other input accepted by \code{calc_REO_matrix()}.}
\item{rank_threshold_grid}{REO threshold values to evaluate.}
\item{lcs_threshold}{LCS threshold used to call active features.}
}
\details{
\code{bootstrap_celltype_communication()} uses a lightweight binomial edge
bootstrap based on observed LCS and edge counts. \code{permute_celltype_communication()}
shuffles cell-type labels while preserving the graph and REO states. It also reports p-value resolution, degenerate zero-variance nulls, and biological interpretation text for structural-null cases. 
\code{sensitivity_REO_threshold()} recomputes REO and LCS across multiple rank
thresholds to identify stable findings.
}
\value{
Data frames for bootstrap and permutation summaries; \code{sensitivity_REO_threshold()}
returns a list with per-threshold results, a stability table, and an active matrix.
}
