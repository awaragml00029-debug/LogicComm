# R/celltype_communication.R

#' Summarize Cell-Type-Resolved Logic Communication
#'
#' Aggregates ligand-receptor Logic Consensus Scores (LCS) by directed
#' sender-cell-type to receiver-cell-type relationships.
#'
#' @param reo_mat Binary REO matrix (genes x cells).
#' @param cell_labels Named character/factor vector mapping cells to cell types.
#' @param seurat_obj Optional Seurat object.
#' @param label_col Optional Seurat metadata column for labels.
#' @param knn_mat Optional KNN adjacency matrix.
#' @param lr_db LogicComm LR database.
#' @param graph_name Optional Seurat graph name.
#' @param lcs_threshold Minimum cell-type LCS.
#' @param min_edges Minimum cell-cell edges.
#' @param min_active_edges Minimum active edges.
#' @param mode Scoring mode.
#' @param remove_self_edges Remove self-loops.
#' @param graph_symmetrize Graph symmetrization.
#' @param edge_weight_mode Edge weighting.
#' @param min_role_hub_quantile Hub score quantile.
#' @param min_role_event_count Minimum event count.
#' @param include_self Include same-type communication.
#' @param verbose Print progress.
#'
#' @return A list of class \code{LogicCommCellTypeComm}.
#' @export
summarize_celltype_communication <- function(reo_mat,
                                             cell_labels = NULL,
                                             seurat_obj = NULL,
                                             label_col = NULL,
                                             knn_mat = NULL,
                                             lr_db = lr_pairs_human,
                                             graph_name = NULL,
                                             lcs_threshold = 0.01,
                                             min_edges = 20,
                                             min_active_edges = 1,
                                             mode = c("auto", "neighborhood", "global"),
                                             remove_self_edges = TRUE,
                                             graph_symmetrize = c("none", "or", "max"),
                                             edge_weight_mode = c("binary", "weighted"),
                                             min_role_hub_quantile = 0.2,
                                             min_role_event_count = 5,
                                             include_self = TRUE,
                                             verbose = TRUE) {
  
  # ── 1. Resolve and Validate Inputs ──────────────────────────────────────────
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_mat <- reo_mat$logic
  mode <- match.arg(mode)
  graph_symmetrize <- match.arg(graph_symmetrize)
  edge_weight_mode <- match.arg(edge_weight_mode)

  labels <- .resolve_celltype_labels(cell_labels, seurat_obj, label_col, colnames(reo_mat))
  cell_types <- sort(unique(labels))
  
  if (identical(mode, "auto")) {
    mode <- if (!is.null(knn_mat) || !is.null(seurat_obj)) "neighborhood" else "global"
  }

  # ── 2. Handle Graph/Edges ──────────────────────────────────────────────────
  if (mode == "neighborhood") {
    if (is.null(knn_mat)) knn_mat <- .extract_knn(seurat_obj, graph_name)
    knn_mat <- .validate_and_align_knn(knn_mat, colnames(reo_mat))
    knn_mat <- .symmetrize_sparse_graph(knn_mat, graph_symmetrize)
    knn_edges <- .sparse_to_edges(knn_mat, remove_self_edges, edge_weight_mode)
    
    edge_i <- knn_edges$i
    edge_j <- knn_edges$j
    edge_w <- knn_edges$w
    
    if (!include_self) {
      keep <- labels[edge_i] != labels[edge_j]
      edge_i <- edge_i[keep]; edge_j <- edge_j[keep]; edge_w <- edge_w[keep]
    }
  } else {
    edge_i <- edge_j <- integer(0)
    edge_w <- numeric(0)
  }

  # ── 3. Pre-calculate complex logic ────────────────────────────────────────
  if (verbose) message("[CellTypeComm] Resolving unique complex logic...")
  unique_complexes <- unique(c(lr_db$ligand_genes, lr_db$receptor_genes))
  complex_keys <- vapply(unique_complexes, function(gs) paste(sort(gs), collapse = "|"), character(1))
  
  complex_logic_map <- lapply(unique_complexes, function(gs) .resolve_complex_logic(gs, reo_mat))
  names(complex_logic_map) <- complex_keys
  
  # Pre-calculate cell-type activity sums for all unique complexes
  complex_type_sums <- lapply(complex_logic_map, function(log_vec) {
    tapply(log_vec, labels, sum, na.rm = TRUE)
  })

  # ── 4. Main Scoring Loop ───────────────────────────────────────────────────
  group_defs <- expand.grid(sender_type = cell_types, receiver_type = cell_types, stringsAsFactors = FALSE)
  if (!include_self) group_defs <- group_defs[group_defs$sender_type != group_defs$receiver_type, ]
  group_defs$key <- paste(group_defs$sender_type, group_defs$receiver_type, sep = "|||")
  n_groups <- nrow(group_defs)
  
  # Pre-calculate edge denominators
  if (mode == "neighborhood") {
    edge_keys <- paste(labels[edge_i], labels[edge_j], sep = "|||")
    edge_group_id <- match(edge_keys, group_defs$key)
    valid_edge <- !is.na(edge_group_id)
    edge_group_id <- edge_group_id[valid_edge]
    edge_i <- edge_i[valid_edge]; edge_j <- edge_j[valid_edge]; edge_w <- edge_w[valid_edge]
    
    n_edges_by_group <- tabulate(edge_group_id, nbins = n_groups)
    edge_weight_sum_by_group <- .sum_by_group(edge_w, edge_group_id, n_groups)
  } else {
    n_by_type <- table(labels)
    n_edges_by_group <- vapply(seq_len(n_groups), function(k) {
      s <- group_defs$sender_type[k]; r <- group_defs$receiver_type[k]
      n <- n_by_type[[s]] * n_by_type[[r]]
      if (!include_self && s == r) n <- n * (n - 1)
      as.integer(n)
    }, integer(1))
    edge_weight_sum_by_group <- as.numeric(n_edges_by_group)
  }

  if (verbose) message(sprintf("[CellTypeComm] Scoring %d L-R pairs for %d CT pairs...", nrow(lr_db), n_groups))

  lr_table_list <- lapply(seq_len(nrow(lr_db)), function(idx) {
    lig_genes <- lr_db$ligand_genes[[idx]]
    rec_genes <- lr_db$receptor_genes[[idx]]
    
    if (!any(lig_genes %in% rownames(reo_mat)) || !any(rec_genes %in% rownames(reo_mat))) return(NULL)

    lig_key <- paste(sort(lig_genes), collapse = "|")
    rec_key <- paste(sort(rec_genes), collapse = "|")
    
    lig_logic <- complex_logic_map[[lig_key]]
    rec_logic <- complex_logic_map[[rec_key]]
    
    # 1. Neighborhood scoring (if applicable)
    lcs_neighborhood <- rep(NA_real_, n_groups)
    n_active_neighborhood <- rep(0, n_groups)
    
    if (mode == "neighborhood") {
      active_events <- lig_logic[edge_i] & rec_logic[edge_j]
      n_active_neighborhood <- tabulate(edge_group_id[active_events], nbins = n_groups)
      w_active <- .sum_by_group(edge_w[active_events], edge_group_id[active_events], n_groups)
      
      denom_nb <- if (edge_weight_mode == "weighted") edge_weight_sum_by_group else n_edges_by_group
      num_nb <- if (edge_weight_mode == "weighted") w_active else n_active_neighborhood
      lcs_neighborhood <- ifelse(denom_nb > 0, num_nb / denom_nb, NA_real_)
    }

    # 2. Global/Potential scoring (always calculate to detect distal signals)
    lig_sum <- complex_type_sums[[lig_key]]
    rec_sum <- complex_type_sums[[rec_key]]
    
    # Calculate global potential (product of active fractions)
    # n_edges_by_group for global mode is |Sender| * |Receiver|
    n_by_type <- table(labels)
    n_global_possible <- vapply(seq_len(n_groups), function(k) {
      as.numeric(n_by_type[[group_defs$sender_type[k]]]) * as.numeric(n_by_type[[group_defs$receiver_type[k]]])
    }, numeric(1))
    
    n_active_global <- vapply(seq_len(n_groups), function(k) {
      as.numeric(lig_sum[[group_defs$sender_type[k]]]) * as.numeric(rec_sum[[group_defs$receiver_type[k]]])
    }, numeric(1))
    
    lcs_global <- n_active_global / pmax(n_global_possible, 1)

    # 3. Decision logic for final LCS
    # If neighborhood mode is ON, we use neighborhood LCS as primary but keep global as reference
    lcs_final <- if (mode == "neighborhood") lcs_neighborhood else lcs_global
    
    # Classification of communication range
    # paracrine: high local, high global
    # juxtacrine: high local, low global (very specific to neighbors)
    # endocrine/distal: low local, high global (the "USA-China" case)
    comm_range <- ifelse(is.na(lcs_neighborhood), "unknown",
                  ifelse(lcs_neighborhood >= lcs_threshold & lcs_global >= lcs_threshold, "paracrine",
                  ifelse(lcs_neighborhood >= lcs_threshold & lcs_global < lcs_threshold, "juxtacrine",
                  ifelse(lcs_neighborhood < lcs_threshold & lcs_global >= lcs_threshold, "distal/endocrine", "inactive"))))

    data.frame(
      sender_type = group_defs$sender_type,
      receiver_type = group_defs$receiver_type,
      lr_pair = lr_db$lr_pair[idx],
      pathway = lr_db$pathway[idx],
      lcs = lcs_final,
      lcs_neighborhood = lcs_neighborhood,
      lcs_global = lcs_global,
      communication_range = comm_range,
      n_edges = n_edges_by_group,
      n_active_edges = if (mode == "neighborhood") n_active_neighborhood else n_active_global,
      active = !is.na(lcs_final) & lcs_final >= lcs_threshold & n_active_neighborhood >= min_active_edges,
      stringsAsFactors = FALSE
    )
  })

  lr_table <- do.call(rbind, lr_table_list)
  
  # ── 5. Summarize and Return ────────────────────────────────────────────────
  pair_summary <- .summarize_celltype_pairs(lr_table, group_defs)
  pathway_summary <- .summarize_celltype_pathways(lr_table)
  
  adj_strength <- .pair_summary_to_matrix(pair_summary, cell_types, "sum_lcs")
  adj_count <- .pair_summary_to_matrix(pair_summary, cell_types, "n_active_lr")
  
  role_summary <- .communication_role_summary(adj_strength, adj_count, 
                                              min_role_hub_quantile = min_role_hub_quantile,
                                              min_role_event_count = min_role_event_count)

  res <- list(lr_table = lr_table, pair_summary = pair_summary, 
              pathway_summary = pathway_summary, role_summary = role_summary,
              cell_labels = labels, lr_db = lr_db,
              params = list(mode = mode, lcs_threshold = lcs_threshold, min_edges = min_edges))
  class(res) <- "LogicCommCellTypeComm"
  
  if (verbose) message(sprintf("[CellTypeComm] Done. %d active L-R events.", sum(lr_table$active, na.rm = TRUE)))
  res
}

# Internal helpers -------------------------------------------------------------

.resolve_celltype_labels <- function(cell_labels, seurat_obj, label_col, cell_names) {
  if (inherits(cell_labels, "Seurat")) {
    seurat_obj <- cell_labels
    cell_labels <- NULL
  }
  if (is.null(cell_labels)) {
    if (is.null(seurat_obj)) {
      stop("cell_labels or seurat_obj must be provided for cell-type-resolved communication.")
    }
    if (!is.null(label_col)) {
      if (!label_col %in% colnames(seurat_obj@meta.data)) {
        stop("label_col '", label_col, "' not found in Seurat metadata.")
      }
      cell_labels <- seurat_obj@meta.data[[label_col]]
      names(cell_labels) <- rownames(seurat_obj@meta.data)
    } else {
      cell_labels <- SeuratObject::Idents(seurat_obj)
      names(cell_labels) <- colnames(seurat_obj)
    }
  }
  label_names <- names(cell_labels)
  cell_labels <- as.character(cell_labels)
  if (!is.null(label_names)) names(cell_labels) <- label_names
  if (is.null(names(cell_labels))) {
    if (length(cell_labels) != length(cell_names)) {
      stop("Unnamed cell_labels must have length ncol(reo_mat).")
    }
    names(cell_labels) <- cell_names
  }
  out <- cell_labels[cell_names]
  out
}

.sum_by_group <- function(values, group, n_groups) {
  out <- numeric(n_groups)
  if (length(values) == 0 || length(group) == 0) return(out)
  ss <- tapply(values, group, sum, na.rm = TRUE)
  out[as.integer(names(ss))] <- as.numeric(ss)
  out
}

.summarize_celltype_pairs <- function(lr_table, group_defs) {
  df <- lr_table[lr_table$active %in% TRUE, , drop = FALSE]
  
  # Aggregate by CT pair
  keys <- paste(lr_table$sender_type, lr_table$receiver_type, sep = "|||")
  sp <- split(seq_len(nrow(lr_table)), keys)
  
  out <- lapply(sp, function(ii) {
    sub <- lr_table[ii, , drop = FALSE]
    active_sub <- sub[sub$active %in% TRUE, , drop = FALSE]
    
    top_lr <- if (nrow(active_sub) > 0) active_sub$lr_pair[which.max(active_sub$lcs)] else NA_character_
    top_pathway <- if (nrow(active_sub) > 0) {
      pw_sum <- tapply(active_sub$lcs, active_sub$pathway, sum, na.rm = TRUE)
      names(pw_sum)[which.max(pw_sum)]
    } else NA_character_

    data.frame(
      sender_type = sub$sender_type[1],
      receiver_type = sub$receiver_type[1],
      n_edges = sub$n_edges[1],
      n_active_lr = sum(sub$active),
      n_juxtacrine = sum(sub$active & sub$communication_range == "juxtacrine"),
      n_paracrine = sum(sub$active & sub$communication_range == "paracrine"),
      n_distal = sum(sub$active & sub$communication_range == "distal/endocrine"),
      sum_lcs = sum(active_sub$lcs, na.rm = TRUE),
      top_lr_pair = top_lr,
      top_pathway = top_pathway,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}

.summarize_celltype_pathways <- function(lr_table) {
  df <- lr_table[lr_table$active %in% TRUE, , drop = FALSE]
  if (nrow(df) == 0) return(data.frame())
  
  keys <- paste(df$sender_type, df$receiver_type, df$pathway, sep = "|||")
  sp <- split(seq_len(nrow(df)), keys)
  
  out <- lapply(sp, function(ii) {
    sub <- df[ii, , drop = FALSE]
    data.frame(
      sender_type = sub$sender_type[1],
      receiver_type = sub$receiver_type[1],
      pathway = sub$pathway[1],
      n_active_lr = nrow(sub),
      sum_lcs = sum(sub$lcs, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}

.pair_summary_to_matrix <- function(pair_summary, cell_types, metric) {
  mat <- matrix(0, nrow = length(cell_types), ncol = length(cell_types),
                dimnames = list(cell_types, cell_types))
  if (nrow(pair_summary) == 0) return(mat)
  
  for (i in seq_len(nrow(pair_summary))) {
    s <- pair_summary$sender_type[i]
    r <- pair_summary$receiver_type[i]
    val <- pair_summary[[metric]][i]
    if (!is.na(val)) mat[s, r] <- val
  }
  mat
}

.extract_knn <- function(seurat_obj, graph_name) {
  if (is.null(seurat_obj)) return(NULL)
  if (is.null(graph_name)) {
    graph_name <- grep("snn$", names(seurat_obj@graphs), value = TRUE)[1]
    if (is.na(graph_name)) graph_name <- names(seurat_obj@graphs)[1]
  }
  if (is.null(graph_name) || !graph_name %in% names(seurat_obj@graphs)) {
    stop("Could not find graph in Seurat object. Please specify graph_name.")
  }
  seurat_obj@graphs[[graph_name]]
}

.validate_and_align_knn <- function(knn_mat, cell_names) {
  if (is.null(knn_mat)) return(NULL)
  if (!all(cell_names %in% colnames(knn_mat))) {
    stop("KNN matrix is missing cells from REO matrix.")
  }
  knn_mat[cell_names, cell_names]
}

.symmetrize_sparse_graph <- function(m, mode) {
  if (mode == "none") return(m)
  m_t <- Matrix::t(m)
  if (mode == "or") return((m + m_t) > 0)
  if (mode == "max") return(Matrix::pmax(m, m_t))
  m
}

.sparse_to_edges <- function(m, remove_self, weight_mode) {
  summary_m <- Matrix::summary(m)
  i <- summary_m$i
  j <- summary_m$j
  w <- summary_m$x
  
  if (remove_self) {
    keep <- i != j
    i <- i[keep]; j <- j[keep]; w <- w[keep]
  }
  
  if (weight_mode == "binary") w <- rep(1, length(w))
  
  list(i = i, j = j, w = w)
}

#' @export
print.LogicCommCellTypeComm <- function(x, ...) {
  cat(sprintf("LogicCommCellTypeComm | %d cells | %d cell types | %d active L-R events\n",
              length(x$cell_labels),
              nrow(x$role_summary),
              sum(x$lr_table$active, na.rm = TRUE)))
  invisible(x)
}

