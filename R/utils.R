# R/utils.R — Internal helpers for LogicComm

#' Extract all unique genes from an LR database
#'
#' @param lr_db Data frame with columns \code{ligand_genes} and
#'   \code{receptor_genes} (each a list of character vectors).
#' @param include_modulators Logical; if \code{TRUE}, include optional CellChatDB
#'   modulator list columns such as \code{agonist_genes},
#'   \code{antagonist_genes}, \code{co_A_receptor_genes}, and
#'   \code{co_I_receptor_genes} when present.
#' @return Character vector of unique gene symbols.
#' @export
#' @examples
#' data(lr_pairs_human)
#' genes <- all_lr_genes(lr_pairs_human)
all_lr_genes <- function(lr_db, include_modulators = FALSE) {
  cols <- c("ligand_genes", "receptor_genes")
  if (isTRUE(include_modulators)) {
    cols <- c(cols, "agonist_genes", "antagonist_genes",
              "co_A_receptor_genes", "co_I_receptor_genes")
  }
  cols <- intersect(cols, names(lr_db))
  genes <- unlist(lr_db[cols], use.names = FALSE)
  sort(unique(genes[nzchar(genes) & !is.na(genes)]))
}

#' Convert a CellChatDB object to a LogicComm LR database
#'
#' Converts the list-like CellChat database object (for example
#' \code{CellChatDB.human}) into the compact LogicComm ligand-receptor schema.
#' The CellChat interaction table is expected at \code{cellchat_db$interaction};
#' complex names referenced by ligand or receptor are resolved against
#' \code{cellchat_db$complex} when present.
#'
#' @param cellchat_db A CellChatDB list with an \code{interaction} data frame and,
#'   optionally, a \code{complex} data frame.
#' @param pathway_col Column name in \code{interaction} to use as pathway. The
#'   CellChat default is \code{"pathway_name"}.
#' @param annotation_col Column name in \code{interaction} to use as annotation.
#'   If absent, annotation is filled with \code{"CellChatDB"}.
#' @param interaction_name_col Column name in \code{interaction} to use as pair
#'   identifier. If absent, identifiers are generated from ligand and receptor.
#' @param source_label Text stored in \code{annotation} when \code{annotation_col}
#'   is absent.
#' @param make_unique Logical; if \code{TRUE}, duplicate \code{lr_pair} names are
#'   made unique with \code{make.unique()}.
#' @param include_modulators Logical; if \code{TRUE}, preserve CellChatDB
#'   agonist/antagonist/co-receptor fields and resolve them to list columns.
#' @return A data frame with LogicComm columns \code{lr_pair}, \code{ligand},
#'   \code{receptor}, \code{ligand_genes}, \code{receptor_genes}, \code{pathway},
#'   and \code{annotation}.
#' @export
#' @examples
#' \dontrun{
#' data(CellChatDB.human, package = "CellChat")
#' lr_db <- as_logiccomm_lr_db_from_cellchat(CellChatDB.human)
#' head(lr_db)
#' }
#' @seealso \code{\link{all_lr_genes}}
as_logiccomm_lr_db_from_cellchat <- function(cellchat_db,
                                             pathway_col = "pathway_name",
                                             annotation_col = "annotation",
                                             interaction_name_col = "interaction_name",
                                             source_label = "CellChatDB",
                                             make_unique = TRUE,
                                             include_modulators = TRUE) {
  if (!is.list(cellchat_db) || is.null(cellchat_db$interaction)) {
    stop("cellchat_db must be a CellChatDB-like list with an 'interaction' data frame.")
  }
  interaction <- as.data.frame(cellchat_db$interaction, stringsAsFactors = FALSE)
  if (!all(c("ligand", "receptor") %in% names(interaction))) {
    stop("cellchat_db$interaction must contain 'ligand' and 'receptor' columns.")
  }

  complex_tbl <- NULL
  if (!is.null(cellchat_db$complex)) {
    complex_tbl <- as.data.frame(cellchat_db$complex, stringsAsFactors = FALSE)
  }

  split_symbol_string <- function(x) {
    x <- as.character(x)
    x <- gsub("\\s+", "", x)
    x <- unlist(strsplit(x, "[+|,;/]", perl = TRUE), use.names = FALSE)
    x <- x[nzchar(x) & !is.na(x)]
    unique(x)
  }

  resolve_partner <- function(x) {
    x <- as.character(x)
    if (is.na(x) || !nzchar(x)) return(character(0))

    if (!is.null(complex_tbl)) {
      rn <- rownames(complex_tbl)
      if (!is.null(rn) && x %in% rn) {
        row <- complex_tbl[x, , drop = FALSE]
        vals <- unlist(row, use.names = FALSE)
        vals <- as.character(vals)
        vals <- vals[nzchar(vals) & !is.na(vals)]
        return(unique(vals))
      }
      if ("complex_name" %in% names(complex_tbl) && x %in% complex_tbl$complex_name) {
        row <- complex_tbl[match(x, complex_tbl$complex_name), , drop = FALSE]
        row$complex_name <- NULL
        vals <- unlist(row, use.names = FALSE)
        vals <- as.character(vals)
        vals <- vals[nzchar(vals) & !is.na(vals)]
        return(unique(vals))
      }
    }

    split_symbol_string(x)
  }

  ligand_genes <- lapply(interaction$ligand, resolve_partner)
  receptor_genes <- lapply(interaction$receptor, resolve_partner)

  pair_name <- if (interaction_name_col %in% names(interaction)) {
    as.character(interaction[[interaction_name_col]])
  } else {
    paste(interaction$ligand, interaction$receptor, sep = "_")
  }
  missing_pair <- is.na(pair_name) | !nzchar(pair_name)
  pair_name[missing_pair] <- paste(interaction$ligand[missing_pair], interaction$receptor[missing_pair], sep = "_")
  pair_name <- gsub("[^A-Za-z0-9_.+-]+", "_", pair_name)
  if (make_unique) pair_name <- make.unique(pair_name, sep = "_")

  pathway <- if (pathway_col %in% names(interaction)) {
    as.character(interaction[[pathway_col]])
  } else {
    rep(NA_character_, nrow(interaction))
  }

  annotation <- if (annotation_col %in% names(interaction)) {
    as.character(interaction[[annotation_col]])
  } else {
    rep(source_label, nrow(interaction))
  }

  out <- data.frame(
    lr_pair = pair_name,
    ligand = as.character(interaction$ligand),
    receptor = as.character(interaction$receptor),
    pathway = pathway,
    annotation = annotation,
    stringsAsFactors = FALSE
  )
  out$ligand_genes <- ligand_genes
  out$receptor_genes <- receptor_genes

  if (isTRUE(include_modulators)) {
    raw_modulator_cols <- c("agonist", "antagonist", "co_A_receptor", "co_I_receptor",
                            "evidence", "interaction_name_2")
    for (col in raw_modulator_cols) {
      out[[col]] <- if (col %in% names(interaction)) {
        as.character(interaction[[col]])
      } else {
        rep(NA_character_, nrow(interaction))
      }
    }

    modulator_cols <- c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")
    for (col in modulator_cols) {
      out[[paste0(col, "_genes")]] <- lapply(out[[col]], resolve_partner)
    }
  }

  core_cols <- c("lr_pair", "ligand", "receptor", "ligand_genes", "receptor_genes", "pathway", "annotation")
  if (isTRUE(include_modulators)) {
    mod_cols <- c("agonist", "antagonist", "co_A_receptor", "co_I_receptor",
                  "agonist_genes", "antagonist_genes", "co_A_receptor_genes",
                  "co_I_receptor_genes", "evidence", "interaction_name_2")
    out <- out[, c(core_cols, mod_cols), drop = FALSE]
  } else {
    out <- out[, core_cols, drop = FALSE]
  }

  empty <- lengths(out$ligand_genes) == 0 | lengths(out$receptor_genes) == 0
  if (any(empty)) {
    warning(sum(empty), " interactions have empty ligand or receptor gene sets after conversion.")
  }
  out
}

#' Safely extract an expression matrix from a Seurat object or matrix-like input
#'
#' Handles both Seurat v4 (slots) and v5 (layers) APIs, and also accepts a plain
#' sparse or dense matrix directly. Sparse matrices are returned without being
#' densified so callers can subset before chunk-wise dense conversion.
#'
#' @param x Seurat object, sparse matrix, or numeric matrix.
#' @param layer Layer/slot to use: \code{"counts"} or \code{"data"}.
#' @return A matrix-like object with genes x cells.
#' @keywords internal
.extract_matrix <- function(x, layer = "counts") {
  if (is.matrix(x)) return(x)
  if (inherits(x, "dgCMatrix") || inherits(x, "sparseMatrix")) return(x)
  if (inherits(x, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE) &&
        !requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat / SeuratObject package required when passing a Seurat object.")
    }
    m <- tryCatch(
      SeuratObject::LayerData(x, layer = layer),
      error = function(e) Seurat::GetAssayData(x, slot = layer)
    )
    return(m)
  }
  stop("Unsupported input type: ", class(x)[1])
}

#' Safely extract KNN graph from a Seurat object
#'
#' @param seurat_obj Seurat object with pre-computed nearest-neighbor graphs.
#' @param graph_name Graph name. If NULL, tries \code{RNA_nn}, \code{SCT_nn},
#'   then the first graph found.
#' @return A sparse \code{dgCMatrix} (cells x cells) adjacency matrix.
#' @keywords internal
.extract_knn <- function(seurat_obj, graph_name = NULL) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  graphs <- names(seurat_obj@graphs)
  if (length(graphs) == 0) {
    stop("No graph found in Seurat object. Run FindNeighbors() first.")
  }
  if (!is.null(graph_name)) {
    if (!graph_name %in% graphs) {
      stop("Graph '", graph_name, "' not found. Available: ", paste(graphs, collapse=", "))
    }
    return(seurat_obj@graphs[[graph_name]])
  }
  # Auto-select: prefer *_nn over *_snn
  nn_graphs <- grep("_nn$", graphs, value = TRUE)
  chosen <- if (length(nn_graphs) > 0) nn_graphs[1] else graphs[1]
  message("Using KNN graph: ", chosen)
  seurat_obj@graphs[[chosen]]
}

#' Compute geometric mean of a matrix column-wise (for multimer subunits)
#'
#' @param mat Numeric matrix (n_subunits x n_cells).
#' @return Numeric vector (n_cells) with geometric means.
#' @keywords internal
.geom_mean_cols <- function(mat) {
  if (is.null(dim(mat)) || nrow(mat) == 1) return(as.vector(mat))
  exp(colMeans(log(mat + 1e-9))) - 1e-9
}

#' Resolve a composite gene set to a single logic vector
#'
#' For ligand and receptor complexes, \code{mode = "all"} is used so every
#' subunit must be active. For modulator sets such as multiple antagonists,
#' \code{mode = "any"} can be used so any available active gene is sufficient.
#'
#' @param genes Character vector of gene names or subunits.
#' @param logic_mat Logical matrix (genes x cells).
#' @param require_all_subunits Logical; if TRUE, missing genes cause all cells to
#'   be returned FALSE.
#' @param mode \code{"all"} requires all available genes to be active;
#'   \code{"any"} requires at least one available gene to be active.
#' @return Logical vector of length \code{ncol(logic_mat)}.
#' @keywords internal
.resolve_complex_logic <- function(genes, logic_mat,
                                   require_all_subunits = TRUE,
                                   mode = c("all", "any")) {
  mode <- match.arg(mode)
  genes <- unique(as.character(genes))
  genes <- genes[nzchar(genes) & !is.na(genes)]
  if (length(genes) == 0) return(rep(FALSE, ncol(logic_mat)))

  available <- intersect(genes, rownames(logic_mat))
  if (length(available) == 0) return(rep(FALSE, ncol(logic_mat)))
  if (isTRUE(require_all_subunits) && length(available) < length(genes)) {
    return(rep(FALSE, ncol(logic_mat)))
  }

  mat_sub <- logic_mat[available, , drop = FALSE]
  if (inherits(mat_sub, "sparseMatrix")) {
    sums <- Matrix::colSums(mat_sub)
  } else {
    sums <- colSums(as.matrix(mat_sub), na.rm = TRUE)
  }

  if (mode == "all") {
    as.logical(sums == nrow(mat_sub))
  } else {
    as.logical(sums > 0)
  }
}

#' BH-adjusted p-values with safe handling
#' @keywords internal
.bh_adjust <- function(pvals) {
  p.adjust(pvals, method = "BH")
}

#' Format a p-value as significance stars
#' @keywords internal
.pval_stars <- function(p) {
  ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**",
      ifelse(p < 0.05, "*", "ns")))
}

#' Null-coalescing operator: return \code{a} unless it is \code{NULL}
#'
#' Package-internal helper. Defined once here and shared across LogicComm so the
#' operator is not duplicated in multiple source files.
#' @param a Value to return when not \code{NULL}.
#' @param b Fallback value used when \code{a} is \code{NULL}.
#' @return \code{a} when it is not \code{NULL}, otherwise \code{b}.
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
