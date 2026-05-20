\name{run_multisample}
\alias{run_multisample}
\alias{get_comparison}
\alias{export_results}
\title{Run multi-sample LogicComm analysis}
\description{One-stop pipeline for REO conversion, LCS scoring, and group comparison.}
\usage{run_multisample(sample_list, group_info, lr_db = lr_pairs_human,
  rank_threshold = 0.5, anchor_genes = NULL, lcs_threshold = 0.01,
  case_label = "Case", ctrl_label = "Ctrl", knn_list = NULL,
  graph_name = NULL, layer = "counts", remove_self_edges = TRUE,
  graph_symmetrize = c("none", "or", "max"),
  edge_weight_mode = c("binary", "weighted"),
  chunk_size = 5000, mc_cores = 1, min_samples_per_group = 1,
  verbose = TRUE)
get_comparison(x)
export_results(result, file, sig_only = FALSE, fdr_cutoff = 0.05)}
\arguments{
\item{sample_list}{Named list of expression inputs.}
\item{group_info}{Group labels for samples.}
\item{lr_db}{Ligand-receptor database.}
\item{rank_threshold}{REO anchor quantile.}
\item{anchor_genes}{Optional anchor gene universe.}
\item{lcs_threshold}{LCS threshold for positive calls.}
\item{case_label}{Case group label.}
\item{ctrl_label}{Control group label.}
\item{knn_list}{Optional named list of KNN/SNN matrices.}
\item{graph_name}{Optional Seurat graph name, e.g. \code{"RNA_nn"} or \code{"RNA_snn"}.}
\item{layer}{Seurat assay layer/slot used for REO conversion.}
\item{remove_self_edges}{Whether to remove KNN self-loops before scoring.}
\item{graph_symmetrize}{How to symmetrize graph entries before scoring.}
\item{edge_weight_mode}{\code{"binary"} treats non-zero graph entries equally; \code{"weighted"} uses graph weights.}
\item{chunk_size}{Cells per REO chunk.}
\item{mc_cores}{Number of forked parallel cores on Unix-like systems.}
\item{min_samples_per_group}{Minimum non-missing samples per group for comparisons.}
\item{verbose}{Whether to print progress messages.}
\item{x}{A \code{LogicCommMulti} object.}
\item{result}{A \code{LogicCommResult} or \code{LogicCommMulti}.}
\item{file}{Output CSV file path.}
\item{sig_only}{Whether to export only significant rows.}
\item{fdr_cutoff}{FDR cutoff for \code{sig_only}.}
}
\value{\code{run_multisample()} returns a \code{LogicCommMulti}; \code{get_comparison()} returns a \code{LogicCommResult}; \code{export_results()} invisibly returns the written data frame.}
\examples{\dontrun{multi <- run_multisample(sample_list, group_info,
  graph_name = "RNA_snn", edge_weight_mode = "weighted", graph_symmetrize = "max")}}
