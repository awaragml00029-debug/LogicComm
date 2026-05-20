\name{logic_api}
\alias{logic_api}
\alias{logic_prepare}
\alias{logic_get_lr_genes}
\alias{logic_import_lrdb}
\alias{logic_check_lrdb}
\alias{logic_score_lr}
\alias{logic_summarize_celltypes}
\alias{logic_compare_groups}
\alias{logic_run}
\alias{print.LogicCommLRDBCheck}
\title{LogicComm main pipeline API}
\description{
A clear pipeline-oriented API for preparing REO logic states, scoring ligand-receptor
communication, summarizing cell-type communication, comparing sample groups, and
running multi-sample LogicComm analyses.
}
\usage{
logic_prepare(object, genes = NULL, lr_db = lr_pairs_human,
  include_modulators = FALSE, rank_threshold = 0.5, layer = "counts",
  chunk_size = 5000, return_rank = FALSE,
  rank_output = c("percentile", "none"), anchor_genes = NULL,
  verbose = TRUE)

logic_get_lr_genes(lr_db = lr_pairs_human, include_modulators = FALSE)

logic_import_lrdb(x, source = c("cellchat", "logiccomm"), ...)

logic_check_lrdb(lr_db, require_metadata = FALSE, stop_on_error = FALSE)

logic_score_lr(reo, seurat_obj = NULL, knn_mat = NULL,
  lr_db = lr_pairs_human, graph_name = NULL,
  mode = c("auto", "neighborhood", "global"), gates = FALSE,
  lcs_threshold = 0.01, remove_self_edges = TRUE,
  graph_symmetrize = NULL, edge_weight_mode = c("binary", "weighted"),
  verbose = TRUE, ...)

logic_summarize_celltypes(reo, cell_labels = NULL, seurat_obj = NULL,
  label_col = NULL, knn_mat = NULL, lr_db = lr_pairs_human,
  graph_name = NULL, mode = c("auto", "neighborhood", "global"),
  include_self = TRUE, remove_self_edges = TRUE, ...)

logic_compare_groups(x, group_info = NULL, ...)

logic_run(samples, group_info, ...)

\method{print}{LogicCommLRDBCheck}(x, ...)
}
\arguments{
\item{object}{Expression matrix, sparse matrix, Seurat object, or BPCells input.}
\item{genes}{Optional genes to retain in the REO output. If \code{NULL}, genes are extracted from \code{lr_db}.}
\item{lr_db}{LogicComm ligand-receptor database.}
\item{include_modulators}{Include optional LR modulator genes when extracting genes from \code{lr_db}.}
\item{rank_threshold}{Per-cell REO anchor quantile.}
\item{layer}{Seurat assay layer/slot.}
\item{chunk_size}{Cells per processing chunk.}
\item{return_rank}{Return rank-percentile matrix with the REO logic matrix.}
\item{rank_output}{Rank output mode.}
\item{anchor_genes}{Optional gene universe used for REO anchors.}
\item{verbose}{Whether to print progress messages.}
\item{x}{For \code{logic_import_lrdb()}, a source database object; for \code{logic_compare_groups()}, an LCS list or LogicComm analysis object; for the print method, a \code{LogicCommLRDBCheck}.}
\item{source}{Input LR database source format. Currently \code{"cellchat"} and \code{"logiccomm"} are supported.}
\item{...}{Additional arguments passed to the underlying implementation.}
\item{require_metadata}{Require recommended metadata columns during LR database validation.}
\item{stop_on_error}{Stop when the LR database has required-field or schema problems.}
\item{reo}{REO matrix or \code{LogicCommREOResult}.}
\item{seurat_obj}{Optional Seurat object for graph extraction.}
\item{knn_mat}{Optional cells x cells graph.}
\item{graph_name}{Optional Seurat graph name.}
\item{mode}{\code{"auto"}, \code{"neighborhood"}, or \code{"global"}.}
\item{gates}{Use gate-aware LR scoring via \code{IdentifyLogicGateConsensus()}.}
\item{lcs_threshold}{LCS activity threshold.}
\item{remove_self_edges}{Remove cell-level graph self-loops. This does not remove same-cell-type communication.}
\item{graph_symmetrize}{Optional graph symmetrization mode. If \code{NULL}, the scorer default is used.}
\item{edge_weight_mode}{\code{"binary"} treats each non-zero graph entry as one edge; \code{"weighted"} uses graph weights.}
\item{cell_labels}{Named vector/factor mapping cells to cell types.}
\item{label_col}{Optional Seurat metadata column used as cell type labels.}
\item{include_self}{Include same-cell-type sender/receiver pairs, interpreted as autocrine-like or homotypic signaling potential.}
\item{group_info}{Named group labels used for sample-level comparison.}
\item{samples}{Named list of expression matrices or Seurat objects.}
}
\details{
The \code{logic_*} API is the recommended user-facing workflow:
\enumerate{
\item \code{logic_prepare()} converts expression to REO logic states.
\item \code{logic_score_lr()} scores ligand-receptor logic as an LCS vector or gate-aware table.
\item \code{logic_summarize_celltypes()} summarizes sender-cell-type to receiver-cell-type communication.
\item \code{logic_compare_groups()} compares sample-level LCS features between groups.
\item \code{logic_run()} runs the multi-sample convenience pipeline.
}

Cell-level self-loops and same-cell-type communication are distinct. The default
\code{remove_self_edges = TRUE} removes diagonal graph entries such as cell A ->
cell A from local graph scoring, while \code{include_self = TRUE} keeps
celltypeA -> celltypeA results for autocrine-like or homotypic interpretation.
}
\value{
\code{logic_prepare()} returns a REO matrix or \code{LogicCommREOResult}.
\code{logic_get_lr_genes()} returns a character vector.
\code{logic_import_lrdb()} returns a LogicComm-format LR database.
\code{logic_check_lrdb()} returns a \code{LogicCommLRDBCheck} object.
\code{logic_score_lr()} returns an \code{LCSVector} or \code{LogicGateResult}.
\code{logic_summarize_celltypes()} returns a \code{LogicCommCellTypeComm}.
\code{logic_compare_groups()} returns a \code{LogicCommResult}.
\code{logic_run()} returns a \code{LogicCommAnalysis} / \code{LogicCommMulti} object.
}
\examples{
\dontrun{
genes <- logic_get_lr_genes(lr_pairs_human)
reo <- logic_prepare(pbmc, genes = genes, layer = "counts")

lcs <- logic_score_lr(
  reo,
  seurat_obj = pbmc,
  graph_name = "RNA_nn",
  lr_db = lr_pairs_human
)

ct_comm <- logic_summarize_celltypes(
  reo,
  seurat_obj = pbmc,
  cell_labels = celltype_labels,
  graph_name = "RNA_nn",
  lr_db = lr_pairs_human
)

comparison <- logic_compare_groups(lcs_list, group_info = group_info)
}
}
