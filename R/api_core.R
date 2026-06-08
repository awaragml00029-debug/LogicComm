# R/api_core.R

#' Prepare REO logic input for LogicComm
#'
#' @param object Expression matrix, sparse matrix, Seurat object, or BPCells input.
#' @param genes Optional genes to retain. If \code{NULL}, genes are extracted from
#'   \code{lr_db}; set both \code{genes = NULL} and \code{lr_db = NULL} to retain
#'   all genes.
#' @param lr_db Optional LogicComm ligand-receptor database used to choose retained
#'   genes when \code{genes} is \code{NULL}.
#' @param include_modulators Include LR modulator genes when extracting genes from
#'   \code{lr_db}.
#' @param rank_threshold Per-cell REO anchor quantile.
#' @param layer Seurat assay layer/slot.
#' @param chunk_size Cells per processing chunk.
#' @param return_rank Return rank-percentile matrix with REO logic.
#' @param rank_output Rank output mode.
#' @param anchor_genes Optional gene universe for REO anchors.
#' @param gene_background Cell-type-free ambient guard passed to
#'   \code{\link{calc_REO_matrix}}: \code{"none"} or \code{"quantile"}.
#' @param gene_background_quantile Across-cell background quantile used when
#'   \code{gene_background = "quantile"}.
#' @param verbose Print progress messages.
#' @return A REO sparse matrix, or a \code{LogicCommREOResult} when
#'   \code{return_rank = TRUE}.
#' @export
logic_prepare <- function(object,
                          genes = NULL,
                          lr_db = lr_pairs_human,
                          include_modulators = FALSE,
                          rank_threshold = 0.5,
                          layer = "counts",
                          chunk_size = 5000,
                          return_rank = FALSE,
                          rank_output = c("percentile", "none"),
                          anchor_genes = NULL,
                          gene_background = c("none", "quantile"),
                          gene_background_quantile = 0.5,
                          verbose = TRUE) {
  rank_output <- match.arg(rank_output)
  gene_background <- match.arg(gene_background)
  if (is.null(genes) && !is.null(lr_db)) {
    genes <- logic_get_lr_genes(lr_db, include_modulators = include_modulators)
  }

  calc_REO_matrix(
    object,
    lr_genes = genes,
    rank_threshold = rank_threshold,
    layer = layer,
    chunk_size = chunk_size,
    return_rank = return_rank,
    rank_output = rank_output,
    anchor_genes = anchor_genes,
    gene_background = gene_background,
    gene_background_quantile = gene_background_quantile,
    verbose = verbose
  )
}

#' Extract LR genes for LogicComm scoring
#'
#' @param lr_db LogicComm ligand-receptor database.
#' @param include_modulators Include optional modulator genes when present.
#' @return Character vector of unique gene symbols.
#' @export
logic_get_lr_genes <- function(lr_db = lr_pairs_human, include_modulators = FALSE) {
  logic_check_lrdb(lr_db, stop_on_error = TRUE)
  all_lr_genes(lr_db, include_modulators = include_modulators)
}

#' Import an LR database into LogicComm format
#'
#' @param x Source database object.
#' @param source Source format. Currently \code{"cellchat"} converts a CellChatDB
#'   object; \code{"logiccomm"} validates and returns an existing LogicComm-format
#'   database.
#' @param ... Arguments passed to the source-specific importer.
#' @return LogicComm-format ligand-receptor database.
#' @export
logic_import_lrdb <- function(x, source = c("cellchat", "logiccomm"), ...) {
  source <- match.arg(source)
  if (source == "cellchat") {
    out <- as_logiccomm_lr_db_from_cellchat(x, ...)
    logic_check_lrdb(out, stop_on_error = TRUE)
    return(out)
  }
  logic_check_lrdb(x, stop_on_error = TRUE)
  x
}

#' Check a LogicComm ligand-receptor database
#'
#' @param lr_db Candidate ligand-receptor database.
#' @param require_metadata Require recommended metadata columns \code{ligand},
#'   \code{receptor}, \code{pathway}, and \code{annotation}.
#' @param stop_on_error Stop when required fields are missing or malformed.
#' @return A \code{LogicCommLRDBCheck} list with validation details.
#' @export
logic_check_lrdb <- function(lr_db, require_metadata = FALSE, stop_on_error = FALSE) {
  required <- c("lr_pair", "ligand_genes", "receptor_genes")
  metadata <- c("ligand", "receptor", "pathway", "annotation")

  is_df <- is.data.frame(lr_db)
  names_lr <- if (is_df) names(lr_db) else character(0)
  missing_required <- setdiff(required, names_lr)
  missing_metadata <- setdiff(metadata, names_lr)

  list_columns <- character(0)
  empty_ligand <- NA_integer_
  empty_receptor <- NA_integer_
  duplicated_pairs <- NA_integer_

  if (is_df && length(missing_required) == 0) {
    gene_cols <- c("ligand_genes", "receptor_genes")
    list_columns <- gene_cols[!vapply(lr_db[gene_cols], is.list, logical(1))]
    ligand_lengths <- lengths(lr_db$ligand_genes)
    receptor_lengths <- lengths(lr_db$receptor_genes)
    empty_ligand <- sum(ligand_lengths == 0 | is.na(ligand_lengths))
    empty_receptor <- sum(receptor_lengths == 0 | is.na(receptor_lengths))
    duplicated_pairs <- sum(duplicated(lr_db$lr_pair))
  }

  ok <- is_df &&
    length(missing_required) == 0 &&
    length(list_columns) == 0 &&
    empty_ligand == 0 &&
    empty_receptor == 0 &&
    (!isTRUE(require_metadata) || length(missing_metadata) == 0)

  problems <- character(0)
  if (!is_df) problems <- c(problems, "lr_db is not a data.frame")
  if (length(missing_required) > 0) {
    problems <- c(problems, paste0("missing required columns: ", paste(missing_required, collapse = ", ")))
  }
  if (length(list_columns) > 0) {
    problems <- c(problems, paste0("columns must be list columns: ", paste(list_columns, collapse = ", ")))
  }
  if (is.finite(empty_ligand) && empty_ligand > 0) {
    problems <- c(problems, paste0(empty_ligand, " rows have empty ligand_genes"))
  }
  if (is.finite(empty_receptor) && empty_receptor > 0) {
    problems <- c(problems, paste0(empty_receptor, " rows have empty receptor_genes"))
  }
  if (isTRUE(require_metadata) && length(missing_metadata) > 0) {
    problems <- c(problems, paste0("missing metadata columns: ", paste(missing_metadata, collapse = ", ")))
  }

  out <- structure(
    list(
      ok = ok,
      n_pairs = if (is_df) nrow(lr_db) else NA_integer_,
      missing_required = missing_required,
      missing_metadata = missing_metadata,
      malformed_list_columns = list_columns,
      empty_ligand_rows = empty_ligand,
      empty_receptor_rows = empty_receptor,
      duplicated_lr_pairs = duplicated_pairs,
      problems = problems
    ),
    class = "LogicCommLRDBCheck"
  )

  if (isTRUE(stop_on_error) && !ok) {
    stop("Invalid LogicComm LR database: ", paste(problems, collapse = "; "))
  }
  out
}

#' Score ligand-receptor communication logic
#'
#' Computes global REO co-expression Logic Consensus Scores. LogicComm scores
#' communication at the cell-type / co-expression level and no longer aggregates
#' over a per-cell neighborhood graph.
#'
#' @param reo REO matrix or \code{LogicCommREOResult}.
#' @param lr_db LogicComm ligand-receptor database.
#' @param gates Use gate-aware LR scoring via \code{IdentifyLogicGateConsensus()}.
#' @param lcs_threshold LCS activity threshold.
#' @param verbose Print progress messages.
#' @param ... For the ordinary scorer, deprecated neighborhood arguments accepted
#'   for backward compatibility and ignored with a warning; for \code{gates =
#'   TRUE}, additional arguments passed to \code{IdentifyLogicGateConsensus()}.
#' @return \code{LCSVector} for ordinary scoring, or \code{LogicGateResult} for
#'   gate-aware scoring.
#' @export
logic_score_lr <- function(reo,
                           lr_db = lr_pairs_human,
                           gates = FALSE,
                           lcs_threshold = 0.01,
                           verbose = TRUE,
                           ...) {
  logic_check_lrdb(lr_db, stop_on_error = TRUE)

  reo_mat <- if (isTRUE(gates)) reo else if (is.list(reo) && !is.null(reo$logic)) reo$logic else reo
  args <- c(list(reo_mat = reo_mat, lr_db = lr_db,
                 lcs_threshold = lcs_threshold, verbose = verbose), list(...))

  if (isTRUE(gates)) {
    return(do.call(IdentifyLogicGateConsensus, args))
  }
  .deprecate_neighborhood_args(list(...), "logic_score_lr")
  do.call(IdentifyLogicConsensus, args)
}

#' Summarize cell-type-resolved LogicComm communication
#'
#' @param reo REO matrix or \code{LogicCommREOResult}.
#' @param cell_labels Named cell type labels.
#' @param seurat_obj Optional Seurat object.
#' @param label_col Optional Seurat metadata column.
#' @param knn_mat Optional cells x cells graph.
#' @param lr_db LogicComm ligand-receptor database.
#' @param graph_name Optional Seurat graph name.
#' @param mode \code{"auto"}, \code{"neighborhood"}, or \code{"global"}.
#' @param include_self Include same-cell-type sender/receiver pairs.
#' @param remove_self_edges Remove cell-level graph self-loops.
#' @param ... Additional arguments passed to \code{summarize_celltype_communication()}.
#' @return \code{LogicCommCellTypeComm} object.
#' @export
logic_summarize_celltypes <- function(reo,
                                      cell_labels = NULL,
                                      seurat_obj = NULL,
                                      label_col = NULL,
                                      knn_mat = NULL,
                                      lr_db = lr_pairs_human,
                                      graph_name = NULL,
                                      mode = c("auto", "neighborhood", "global"),
                                      include_self = TRUE,
                                      remove_self_edges = TRUE,
                                      ...) {
  mode <- match.arg(mode)
  logic_check_lrdb(lr_db, stop_on_error = TRUE)
  summarize_celltype_communication(
    reo_mat = reo,
    cell_labels = cell_labels,
    seurat_obj = seurat_obj,
    label_col = label_col,
    knn_mat = knn_mat,
    lr_db = lr_db,
    graph_name = graph_name,
    mode = mode,
    include_self = include_self,
    remove_self_edges = remove_self_edges,
    ...
  )
}

#' Compare LogicComm scores between sample groups
#'
#' @param x Named list of LCS vectors, \code{LogicCommMulti}, or
#'   \code{LogicCommAnalysis}.
#' @param group_info Named group labels. Optional when \code{x} already contains a
#'   comparison result.
#' @param ... Additional arguments passed to \code{CompareLogicGroups()}.
#' @return \code{LogicCommResult}.
#' @export
logic_compare_groups <- function(x, group_info = NULL, ...) {
  if (inherits(x, "LogicCommAnalysis") || inherits(x, "LogicCommMulti")) {
    if (is.null(group_info) && !is.null(x$comparison)) return(x$comparison)
    if (is.null(x$lcs_list)) stop("analysis object does not contain lcs_list.")
    x <- x$lcs_list
  }
  if (is.null(group_info)) stop("group_info is required when x is an LCS list.")
  CompareLogicGroups(x, group_info = group_info, ...)
}

#' Run the LogicComm multi-sample pipeline
#'
#' @param samples Named list of expression matrices or Seurat objects.
#' @param group_info Named group labels.
#' @param ... Additional arguments passed to \code{run_multisample()}.
#' @return A \code{LogicCommAnalysis} object.
#' @export
logic_run <- function(samples, group_info, ...) {
  out <- run_multisample(sample_list = samples, group_info = group_info, ...)
  class(out) <- unique(c("LogicCommAnalysis", class(out)))
  out
}

#' @export
print.LogicCommLRDBCheck <- function(x, ...) {
  cat(sprintf("LogicCommLRDBCheck: %s", if (isTRUE(x$ok)) "OK" else "FAILED"), "\n", sep = "")
  cat(sprintf("L-R pairs: %s\n", ifelse(is.na(x$n_pairs), "NA", as.character(x$n_pairs))))
  if (length(x$missing_metadata) > 0) {
    cat("Missing optional metadata: ", paste(x$missing_metadata, collapse = ", "), "\n", sep = "")
  }
  if (length(x$problems) > 0) {
    cat("Problems:\n")
    for (p in x$problems) cat("- ", p, "\n", sep = "")
  }
  invisible(x)
}
