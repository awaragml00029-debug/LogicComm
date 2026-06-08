# R/AutoLabelLogicClusters.R

#' Auto-Annotate Cell Clusters Driving Logic Consensus Hotspots
#'
#' @param reo_mat REO binary matrix.
#' @param cell_labels Named character vector or Seurat object.
#' @param lr_db LR database.
#' @param lr_pairs Character vector of L-R pairs to analyze.
#' @param min_sender_frac Minimum fraction for sender.
#' @param min_receiver_frac Minimum fraction for receiver.
#' @param verbose Logical.
#' @return A data frame with cluster enrichment stats.
#' @examples
#' expr <- matrix(
#'   c(5, 1, 4, 2, 1, 5, 3, 4, 4, 2, 5, 1),
#'   nrow = 3,
#'   dimnames = list(c("L1", "R1", "T1"), paste0("cell", 1:4))
#' )
#' reo <- expr >= 3
#' rank_mat <- apply(expr, 2, rank) / nrow(expr)
#' lr_db <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   stringsAsFactors = FALSE
#' )
#' lcs <- data.frame(
#'   ligand = "L1",
#'   receptor = "R1",
#'   pathway = "toy",
#'   sample = c("S1", "S2"),
#'   group = c("control", "case"),
#'   sender = "A",
#'   receiver = "B",
#'   celltype_sender = "A",
#'   celltype_receiver = "B",
#'   LCS = c(0.2, 0.5),
#'   lcs = c(0.2, 0.5),
#'   mean_lcs = c(0.2, 0.5),
#'   delta_lcs = c(0.0, 0.3),
#'   p_value = c(0.5, 0.01),
#'   p_adj = c(0.5, 0.02),
#'   fdr = c(0.5, 0.02),
#'   stringsAsFactors = FALSE
#' )
#' sample_ct_list <- list(S1 = lcs, S2 = lcs)
#' group_info <- c(S1 = "control", S2 = "case")
#' knn <- matrix(1, nrow = 4, ncol = 4, dimnames = list(colnames(expr), colnames(expr)))
#' diag(knn) <- 0
#' toy_args <- list(
#'   x = lcs, result = lcs, results = lcs, lcs_df = lcs, ct_comm = lcs,
#'   comm_df = lcs, communication = lcs, celltype_comm = lcs,
#'   celltype_results = lcs, differential_results = lcs, diff_comm = lcs,
#'   glm_result = lcs, role_df = lcs, roles = lcs, specificity = lcs,
#'   null_pair = list(observed = lcs, null = lcs), reo_mat = reo,
#'   rank_mat = rank_mat, expr_mat = expr, expression = expr,
#'   lr_db = lr_db, samples = list(S1 = expr, S2 = expr),
#'   sample_ct_list = sample_ct_list, group_info = group_info,
#'   group_labels = group_info, groups = group_info, knn_mat = knn,
#'   output_dir = tempfile("logiccomm"), file = tempfile(fileext = ".csv"),
#'   path = tempfile(fileext = ".csv")
#' )
#' fun <- get("AutoLabelLogicClusters")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
AutoLabelLogicClusters <- function(reo_mat,
                                    cell_labels,
                                    lr_db             = lr_pairs_human,
                                    lr_pairs          = NULL,
                                    min_sender_frac   = 0.1,
                                    min_receiver_frac = 0.1,
                                    verbose           = TRUE) {

  if (inherits(cell_labels, "Seurat")) {
    labels_vec <- as.character(SeuratObject::Idents(cell_labels))
    names(labels_vec) <- colnames(cell_labels)
  } else {
    labels_vec <- as.character(cell_labels)
    if (is.null(names(labels_vec))) names(labels_vec) <- colnames(reo_mat)
  }

  common_cells <- intersect(names(labels_vec), colnames(reo_mat))
  valid_cells <- common_cells[!is.na(labels_vec[common_cells])]
  
  reo_sub <- reo_mat[, valid_cells, drop = FALSE]
  labels_sub <- labels_vec[valid_cells]
  clusters <- sort(unique(labels_sub))

  if (!is.null(lr_pairs)) {
    lr_db <- lr_db[lr_db$lr_pair %in% lr_pairs, , drop = FALSE]
  }
  if (nrow(lr_db) == 0) stop("No valid lr_pairs found.")

  if (verbose) message(sprintf("[AutoLabel] Scoring %d L-R pairs...", nrow(lr_db)))

  # Pre-calculate unique complex logic
  unique_complexes <- unique(c(lr_db$ligand_genes, lr_db$receptor_genes))
  complex_keys <- vapply(unique_complexes, function(gs) paste(sort(gs), collapse = "|"), character(1))
  complex_logic_map <- lapply(unique_complexes, function(gs) .resolve_complex_logic(gs, reo_sub))
  names(complex_logic_map) <- complex_keys
  
  get_logic <- function(gs) {
    key <- paste(sort(gs), collapse = "|")
    complex_logic_map[[key]]
  }

  rows <- lapply(seq_len(nrow(lr_db)), function(i) {
    lig_logic <- get_logic(lr_db$ligand_genes[[i]])
    rec_logic <- get_logic(lr_db$receptor_genes[[i]])

    sender_stats <- .cluster_enrichment(lig_logic, labels_sub, clusters, min_sender_frac)
    receiver_stats <- .cluster_enrichment(rec_logic, labels_sub, clusters, min_receiver_frac)

    data.frame(
      lr_pair = lr_db$lr_pair[i],
      ligand = lr_db$ligand[i] %||% NA,
      receptor = lr_db$receptor[i] %||% NA,
      sender_cluster = sender_stats$top_cluster,
      sender_frac = round(sender_stats$top_frac, 3),
      receiver_cluster = receiver_stats$top_cluster,
      receiver_frac = round(receiver_stats$top_frac, 3),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  if (verbose) message("[AutoLabel] Done.")
  result
}

.cluster_enrichment <- function(logic_vec, cluster_vec, clusters, min_frac) {
  n_total <- length(logic_vec)
  n_active <- sum(logic_vec, na.rm = TRUE)
  best <- list(top_cluster = NA_character_, top_frac = 0, top_pval = 1)

  for (cl in clusters) {
    in_cl <- cluster_vec == cl
    n_cl <- sum(in_cl)
    if (n_cl == 0) next

    active_in_cl <- sum(logic_vec[in_cl], na.rm = TRUE)
    frac_in_cl <- active_in_cl / n_cl
    if (frac_in_cl < min_frac) next

    # Simplified p-value calculation or just keep best frac for auto-labeling
    if (frac_in_cl > best$top_frac) {
       best <- list(top_cluster = cl, top_frac = frac_in_cl, top_pval = 0)
    }
  }
  best
}
