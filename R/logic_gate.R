# R/logic_gate.R — Rank-aware Logic Gate consensus for LogicComm

#' Validate and align rank matrix to REO matrix
#' @keywords internal
.validate_rank_mat <- function(rank_mat, reo_mat) {
  if (is.null(rank_mat)) stop("rank_mat is required for rank-aware logic-gate scoring.")
  if (is.list(rank_mat) && !is.null(rank_mat$rank)) rank_mat <- rank_mat$rank
  if (is.null(dim(rank_mat))) stop("rank_mat must be a genes x cells matrix.")
  if (is.null(rownames(rank_mat)) || is.null(colnames(rank_mat))) {
    stop("rank_mat must have rownames and colnames.")
  }
  missing_genes <- setdiff(rownames(reo_mat), rownames(rank_mat))
  missing_cells <- setdiff(colnames(reo_mat), colnames(rank_mat))
  if (length(missing_genes) > 0 || length(missing_cells) > 0) {
    stop("rank_mat must contain all reo_mat genes and cells.")
  }
  rank_mat[rownames(reo_mat), colnames(reo_mat), drop = FALSE]
}

#' Compute active logic for optional positive modulators
#' @keywords internal
.positive_modulator_gate <- function(genes, reo_mat, mode = c("all", "any", "ignore")) {
  mode <- match.arg(mode)
  if (mode == "ignore" || length(genes) == 0) return(rep(TRUE, ncol(reo_mat)))
  if (mode == "all") {
    .resolve_complex_logic(genes, reo_mat, require_all_subunits = TRUE, mode = "all")
  } else {
    .resolve_complex_logic(genes, reo_mat, require_all_subunits = FALSE, mode = "any")
  }
}

#' Rank-aware negative gate for antagonist or inhibitory co-receptor genes
#' @keywords internal
.rank_block_gate_v2 <- function(reference_rank, inhibitor_genes, get_rank, get_logic, reo_mat,
                              fallback = c("active_block", "ignore")) {
  fallback <- match.arg(fallback)
  if (length(inhibitor_genes) == 0) return(rep(TRUE, ncol(reo_mat)))
  inhibitor_rank <- get_rank(inhibitor_genes)
  
  if (is.null(inhibitor_rank)) {
     if (fallback == "ignore") return(rep(TRUE, ncol(reo_mat)))
     inhibitor_active <- get_logic(inhibitor_genes, mode = "any") %||% rep(FALSE, ncol(reo_mat))
     return(!inhibitor_active)
  }
  
  gate <- is.na(inhibitor_rank) | is.na(reference_rank) | reference_rank > inhibitor_rank
  as.logical(gate)
}

#' Compute active logic for optional positive modulators (v2)
#' @keywords internal
.positive_modulator_gate_v2 <- function(genes, get_logic, mode = c("all", "any", "ignore")) {
  mode <- match.arg(mode)
  if (mode == "ignore" || length(genes) == 0) return(NULL)
  get_logic(genes, mode = mode)
}

#' Identify rank-aware logic-gate consensus scores
#'
#' Extends the base ligand-receptor Logic Consensus Score with optional modulators.
#'
#' @param reo_mat Binary REO matrix.
#' @param rank_mat Rank-percentile matrix.
#' @param seurat_obj Optional Seurat object.
#' @param knn_mat Optional KNN adjacency matrix.
#' @param lr_db LogicComm LR database.
#' @param graph_name Name of Seurat graph.
#' @param positive_gate How to handle positive modulators: \code{"all"}, \code{"any"}, \code{"ignore"}.
#' @param negative_gate How to handle negative modulators: \code{"rank_block"}, \code{"active_block"}, \code{"ignore"}.
#' @param agonist_side Side used for agonist gates: sender, receiver, either, or ignore.
#' @param antagonist_side Side used for antagonist blocking: sender, receiver, both, either, or ignore.
#' @param lcs_threshold Threshold for verbose summary.
#' @param remove_self_edges Logical.
#' @param graph_symmetrize \code{"none"}, \code{"or"}, or \code{"max"}.
#' @param edge_weight_mode \code{"binary"} or \code{"weighted"}.
#' @param verbose Print progress messages.
#' @return A data frame with base and regulated LCS plus diagnostic rates.
#' @export
IdentifyLogicGateConsensus <- function(reo_mat,
                                       rank_mat = NULL,
                                       seurat_obj = NULL,
                                       knn_mat = NULL,
                                       lr_db = lr_pairs_human,
                                       graph_name = NULL,
                                       positive_gate = c("all", "any", "ignore"),
                                       negative_gate = c("rank_block", "active_block", "ignore"),
                                       agonist_side = c("receiver", "sender", "either", "ignore"),
                                       antagonist_side = c("both", "sender", "receiver", "either", "ignore"),
                                       lcs_threshold = 0.01,
                                       remove_self_edges = TRUE,
                                       graph_symmetrize = c("none", "or", "max"),
                                       edge_weight_mode = c("binary", "weighted"),
                                       verbose = TRUE) {
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) {
    if (is.null(rank_mat) && !is.null(reo_mat$rank)) rank_mat <- reo_mat$rank
    reo_mat <- reo_mat$logic
  }
  stopifnot(inherits(reo_mat, c("dgCMatrix", "sparseMatrix", "matrix")))
  positive_gate <- match.arg(positive_gate)
  negative_gate <- match.arg(negative_gate)
  agonist_side <- match.arg(agonist_side)
  antagonist_side <- match.arg(antagonist_side)
  graph_symmetrize <- match.arg(graph_symmetrize)
  edge_weight_mode <- match.arg(edge_weight_mode)

  if (negative_gate == "rank_block" || !is.null(rank_mat)) {
    rank_mat <- .validate_rank_mat(rank_mat, reo_mat)
  }

  use_knn <- FALSE
  knn_edges <- NULL
  if (!is.null(knn_mat)) {
    use_knn <- TRUE
    knn_mat <- .validate_and_align_knn(knn_mat, colnames(reo_mat))
    knn_mat <- .symmetrize_sparse_graph(knn_mat, graph_symmetrize)
    knn_edges <- .sparse_to_edges(knn_mat, remove_self_edges = remove_self_edges,
                                   edge_weight_mode = edge_weight_mode)
  } else if (!is.null(seurat_obj)) {
    knn_mat <- .extract_knn(seurat_obj, graph_name)
    knn_mat <- .validate_and_align_knn(knn_mat, colnames(reo_mat))
    knn_mat <- .symmetrize_sparse_graph(knn_mat, graph_symmetrize)
    knn_edges <- .sparse_to_edges(knn_mat, remove_self_edges = remove_self_edges,
                                   edge_weight_mode = edge_weight_mode)
    use_knn <- TRUE
  }
  
  if (verbose) {
    message(if (use_knn) sprintf("[LogicGate] KNN mode with %d edges.", length(knn_edges$i)) else "[LogicGate] Global mode.")
    message(sprintf("[LogicGate] Scoring %d L-R pairs with positive_gate='%s', negative_gate='%s'.",
                    nrow(lr_db), positive_gate, negative_gate))
  }

  # ── Pre-calculate unique complex logic ──────────────────────────────────────
  if (verbose) message("[LogicGate] Pre-calculating unique complex logic...")
  all_genes_sets <- unique(c(
    lr_db$ligand_genes, lr_db$receptor_genes,
    lr_db$agonist_genes, lr_db$antagonist_genes,
    lr_db$co_A_receptor_genes, lr_db$co_I_receptor_genes
  ))
  all_genes_sets <- all_genes_sets[lengths(all_genes_sets) > 0]
  
  complex_keys <- vapply(all_genes_sets, function(gs) paste(sort(gs), collapse = "|"), character(1))
  
  # Logic map handles both ALL and ANY modes
  logic_map_all <- lapply(all_genes_sets, function(gs) .resolve_complex_logic(gs, reo_mat, mode = "all"))
  logic_map_any <- lapply(all_genes_sets, function(gs) .resolve_complex_logic(gs, reo_mat, mode = "any"))
  names(logic_map_all) <- names(logic_map_any) <- complex_keys
  
  get_logic <- function(gs, mode = c("all", "any")) {
    if (length(gs) == 0) return(NULL)
    mode <- match.arg(mode)
    key <- paste(sort(gs), collapse = "|")
    if (mode == "all") logic_map_all[[key]] else logic_map_any[[key]]
  }
  
  rank_map <- list()
  if (!is.null(rank_mat)) {
    if (verbose) message("[LogicGate] Pre-calculating unique complex ranks...")
    rank_map <- lapply(all_genes_sets, function(gs) .resolve_complex_rank(gs, rank_mat))
    names(rank_map) <- complex_keys
  }
  
  get_rank <- function(gs) {
    if (length(gs) == 0) return(NULL)
    key <- paste(sort(gs), collapse = "|")
    rank_map[[key]]
  }

  out <- vector("list", nrow(lr_db))
  for (idx in seq_len(nrow(lr_db))) {
    lig_genes <- .lr_list_entry(lr_db, "ligand_genes", idx)
    rec_genes <- .lr_list_entry(lr_db, "receptor_genes", idx)
    
    lig_logic <- get_logic(lig_genes, mode = "all")
    rec_logic <- get_logic(rec_genes, mode = "all")
    if (is.null(lig_logic) || is.null(rec_logic)) {
      out[[idx]] <- .empty_gate_row(lr_db$lr_pair[idx], lr_db, idx)
      next
    }
    
    lig_rank <- get_rank(lig_genes) %||% rep(NA_real_, ncol(reo_mat))
    rec_rank <- get_rank(rec_genes) %||% rep(NA_real_, ncol(reo_mat))

    agonist_genes <- .lr_list_entry(lr_db, "agonist_genes", idx)
    antagonist_genes <- .lr_list_entry(lr_db, "antagonist_genes", idx)
    coa_genes <- .lr_list_entry(lr_db, "co_A_receptor_genes", idx)
    coi_genes <- .lr_list_entry(lr_db, "co_I_receptor_genes", idx)

    sender_gate <- rep(TRUE, ncol(reo_mat))
    receiver_gate <- rep(TRUE, ncol(reo_mat))

    if (length(coa_genes) > 0) {
      receiver_gate <- receiver_gate & (.positive_modulator_gate_v2(coa_genes, get_logic, positive_gate) %||% TRUE)
    }
    if (length(agonist_genes) > 0 && agonist_side != "ignore") {
      ag_gate <- .positive_modulator_gate_v2(agonist_genes, get_logic, positive_gate) %||% TRUE
      if (agonist_side == "sender") sender_gate <- sender_gate & ag_gate
      if (agonist_side == "receiver") receiver_gate <- receiver_gate & ag_gate
      if (agonist_side == "either") {
        sender_gate <- sender_gate & ag_gate
        receiver_gate <- receiver_gate & ag_gate
      }
    }

    if (negative_gate != "ignore" && length(coi_genes) > 0) {
      if (negative_gate == "rank_block") {
        receiver_gate <- receiver_gate & .rank_block_gate_v2(rec_rank, coi_genes, get_rank, get_logic, reo_mat)
      } else {
        receiver_gate <- receiver_gate & ! (get_logic(coi_genes, mode = "any") %||% rep(FALSE, ncol(reo_mat)))
      }
    }
    if (negative_gate != "ignore" && length(antagonist_genes) > 0 && antagonist_side != "ignore") {
      if (negative_gate == "rank_block") {
        ant_sender <- .rank_block_gate_v2(lig_rank, antagonist_genes, get_rank, get_logic, reo_mat)
        ant_receiver <- .rank_block_gate_v2(rec_rank, antagonist_genes, get_rank, get_logic, reo_mat)
      } else {
        ant_active <- get_logic(antagonist_genes, mode = "any") %||% rep(FALSE, ncol(reo_mat))
        ant_sender <- !ant_active
        ant_receiver <- !ant_active
      }
      if (antagonist_side == "sender") sender_gate <- sender_gate & ant_sender
      if (antagonist_side == "receiver") receiver_gate <- receiver_gate & ant_receiver
      if (antagonist_side == "both") {
        sender_gate <- sender_gate & ant_sender
        receiver_gate <- receiver_gate & ant_receiver
      }
      if (antagonist_side == "either") {
        either_gate <- ant_sender | ant_receiver
        sender_gate <- sender_gate & either_gate
        receiver_gate <- receiver_gate & either_gate
      }
    }

    lig_reg <- lig_logic & sender_gate
    rec_reg <- rec_logic & receiver_gate

    if (use_knn) {
      base_events <- lig_logic[knn_edges$i] & rec_logic[knn_edges$j]
      reg_events <- lig_reg[knn_edges$i] & rec_reg[knn_edges$j]
      if (edge_weight_mode == "weighted") {
        w <- knn_edges$w
        denom <- sum(w, na.rm = TRUE)
        base_lcs <- if (denom > 0) sum(w[base_events], na.rm = TRUE) / denom else NA_real_
        regulated_lcs <- if (denom > 0) sum(w[reg_events], na.rm = TRUE) / denom else NA_real_
        candidate_rate <- base_lcs
        base_weight <- sum(w[base_events], na.rm = TRUE)
        block_rate <- if (base_weight == 0) 0 else sum(w[base_events & !reg_events], na.rm = TRUE) / base_weight
      } else {
        base_lcs <- mean(base_events)
        regulated_lcs <- mean(reg_events)
        candidate_rate <- mean(base_events)
        block_rate <- if (sum(base_events) == 0) 0 else sum(base_events & !reg_events) / sum(base_events)
      }
    } else {
      base_events <- lig_logic & rec_logic
      reg_events <- lig_reg & rec_reg
      base_lcs <- mean(base_events)
      regulated_lcs <- mean(reg_events)
      candidate_rate <- mean(base_events)
      block_rate <- if (sum(base_events) == 0) 0 else sum(base_events & !reg_events) / sum(base_events)
    }

    out[[idx]] <- data.frame(
      lr_pair = lr_db$lr_pair[idx],
      ligand = if ("ligand" %in% names(lr_db)) as.character(lr_db$ligand[idx]) else NA_character_,
      receptor = if ("receptor" %in% names(lr_db)) as.character(lr_db$receptor[idx]) else NA_character_,
      base_lcs = base_lcs,
      regulated_lcs = regulated_lcs,
      delta_lcs = regulated_lcs - base_lcs,
      block_rate = block_rate,
      candidate_event_rate = candidate_rate,
      has_agonist = length(agonist_genes) > 0,
      has_antagonist = length(antagonist_genes) > 0,
      has_co_A_receptor = length(coa_genes) > 0,
      has_co_I_receptor = length(coi_genes) > 0,
      pathway = if ("pathway" %in% names(lr_db)) as.character(lr_db$pathway[idx]) else NA_character_,
      annotation = if ("annotation" %in% names(lr_db)) as.character(lr_db$annotation[idx]) else NA_character_,
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, out)
  class(res) <- c("LogicGateResult", class(res))
  if (verbose) {
    n_active <- sum(!is.na(res$regulated_lcs) & res$regulated_lcs >= lcs_threshold)
    message(sprintf("[LogicGate] Done. %d / %d pairs regulated-active.", n_active, nrow(res)))
  }
  res
}

#' Helper for empty result rows
#' @keywords internal
.empty_gate_row <- function(lr_pair, lr_db, idx) {
  data.frame(
    lr_pair = lr_pair,
    ligand = if ("ligand" %in% names(lr_db)) as.character(lr_db$ligand[idx]) else NA_character_,
    receptor = if ("receptor" %in% names(lr_db)) as.character(lr_db$receptor[idx]) else NA_character_,
    base_lcs = NA_real_, regulated_lcs = NA_real_, delta_lcs = NA_real_, block_rate = NA_real_,
    candidate_event_rate = NA_real_, has_agonist = FALSE, has_antagonist = FALSE,
    has_co_A_receptor = FALSE, has_co_I_receptor = FALSE,
    pathway = if ("pathway" %in% names(lr_db)) as.character(lr_db$pathway[idx]) else NA_character_,
    annotation = if ("annotation" %in% names(lr_db)) as.character(lr_db$annotation[idx]) else NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Print method for LogicGateResult
#' @param x LogicGateResult object.
#' @param ... ignored.
#' @export
print.LogicGateResult <- function(x, ...) {
  n_total <- nrow(x)
  n_scored <- sum(is.finite(x$regulated_lcs), na.rm = TRUE)
  n_blocked <- sum(is.finite(x$block_rate) & x$block_rate > 0, na.rm = TRUE)
  cat(sprintf("LogicGateResult: %d pairs | scored: %d | modulator-blocked: %d\n",
              n_total, n_scored, n_blocked))
  show_cols <- intersect(c("lr_pair", "base_lcs", "regulated_lcs", "delta_lcs", "block_rate"), names(x))
  if (length(show_cols) > 0) print(utils::head(as.data.frame(x)[, show_cols, drop = FALSE], 5), row.names = FALSE)
  invisible(x)
}
