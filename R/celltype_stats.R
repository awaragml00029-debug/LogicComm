# R/celltype_stats.R

#' Bootstrap LogicComm Cell-Type Communication Scores
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param n_boot Number of bootstrap replicates.
#' @param level \code{"celltype_lr"} or \code{"celltype_pair"}.
#' @param seed Optional random seed.
#' @return Data frame with observed score and bootstrap interval.
#' @export
bootstrap_celltype_communication <- function(ct_comm,
                                             n_boot = 100,
                                             level = c("celltype_lr", "celltype_pair"),
                                             seed = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  level <- match.arg(level)
  if (!is.null(seed)) set.seed(seed)
  df <- ct_comm$lr_table
  df <- df[!is.na(df$lcs) & df$n_edges > 0, , drop = FALSE]
  if (nrow(df) == 0) stop("No scored L-R rows with positive edge counts.")
  probs <- pmin(pmax(df$lcs, 0), 1)
  trials <- pmax(1, as.integer(round(df$n_edges)))
  boot_mat <- matrix(NA_real_, nrow = nrow(df), ncol = n_boot)
  for (b in seq_len(n_boot)) boot_mat[, b] <- stats::rbinom(nrow(df), size = trials, prob = probs) / trials
  if (level == "celltype_lr") {
    out <- data.frame(
      sender_type = df$sender_type,
      receiver_type = df$receiver_type,
      lr_pair = df$lr_pair,
      observed_lcs = df$lcs,
      boot_mean = rowMeans(boot_mat, na.rm = TRUE),
      boot_sd = apply(boot_mat, 1, stats::sd, na.rm = TRUE),
      ci_low = apply(boot_mat, 1, stats::quantile, probs = 0.025, na.rm = TRUE),
      ci_high = apply(boot_mat, 1, stats::quantile, probs = 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    return(out[order(out$observed_lcs, decreasing = TRUE), , drop = FALSE])
  }
  keys <- paste(df$sender_type, df$receiver_type, sep = "|||")
  sp <- split(seq_len(nrow(df)), keys)
  out <- lapply(sp, function(ii) {
    bm <- colSums(boot_mat[ii, , drop = FALSE], na.rm = TRUE)
    data.frame(sender_type = df$sender_type[ii[1]], receiver_type = df$receiver_type[ii[1]],
               observed_sum_lcs = sum(df$lcs[ii], na.rm = TRUE),
               boot_mean = mean(bm), boot_sd = stats::sd(bm),
               ci_low = stats::quantile(bm, probs = 0.025, na.rm = TRUE),
               ci_high = stats::quantile(bm, probs = 0.975, na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  out[order(out$observed_sum_lcs, decreasing = TRUE), , drop = FALSE]
}

#' Cell-Label Permutation Null for Cell-Type Communication
#'
#' @param ct_comm Optional observed object. If NULL, it is computed from inputs.
#' @param reo_mat Binary REO matrix or LogicCommREOResult.
#' @param cell_labels Named cell-type labels. Defaults to \code{ct_comm$cell_labels}
#'   when \code{ct_comm} is supplied. If \code{ct_comm} was built after filtering
#'   unlabeled cells, \code{reo_mat} and named \code{knn_mat} inputs are aligned to
#'   those retained cells automatically.
#' @param knn_mat Optional KNN/SNN graph. Named graphs are subset to retained cells
#'   when necessary.
#' @param lr_db LR database. Defaults to \code{ct_comm$lr_db} when available.
#' @param n_perm Number of permutations.
#' @param metric Pair-level metric to test.
#' @param seed Optional random seed.
#' @param adaptive If TRUE, refine top preliminary candidates after the first
#'   \code{n_perm} permutations.
#' @param adaptive_top_n Number of preliminary candidates to refine.
#' @param adaptive_n_perm Target total permutations for adaptive candidates. For
#'   example, \code{n_perm = 100} and \code{adaptive_n_perm = 1000} runs 900
#'   additional full permutations for the top candidates.
#' @param min_n_perm_publication Recommended minimum permutation count.
#' @param verbose Print progress, including adaptive-phase elapsed time.
#' @param ... Passed to \code{summarize_celltype_communication()}.
#' @return Pair-level null summary with empirical p-values.
#' @export
permute_celltype_communication <- function(ct_comm = NULL,
                                           reo_mat,
                                           cell_labels = NULL,
                                           knn_mat = NULL,
                                           lr_db = NULL,
                                           n_perm = 100,
                                           metric = "sum_lcs",
                                           seed = NULL,
                                           adaptive = FALSE,
                                           adaptive_top_n = 25,
                                           adaptive_n_perm = 1000,
                                           min_n_perm_publication = 1000,
                                           verbose = TRUE,
                                           ...) {
  if (!is.null(seed)) set.seed(seed)
  dots <- list(...)
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_mat <- reo_mat$logic
  if (is.null(ct_comm)) {
    ct_comm <- do.call(
      summarize_celltype_communication,
      c(list(reo_mat = reo_mat, cell_labels = cell_labels, knn_mat = knn_mat,
             lr_db = lr_db, verbose = FALSE), dots)
    )
  }
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  inherited_params <- intersect(
    names(ct_comm$params),
    c("mode", "graph_symmetrize", "edge_weight_mode", "remove_self_edges",
      "include_self", "lcs_threshold", "min_edges", "min_active_edges")
  )
  for (param in inherited_params) {
    if (is.null(dots[[param]])) {
      dots[[param]] <- ct_comm$params[[param]]
    } else if (!identical(dots[[param]], ct_comm$params[[param]])) {
      warning(
        "Permutation parameter '", param, "' differs from ct_comm$params$", param,
        "; observed and null summaries may use different scoring definitions.",
        call. = FALSE
      )
    }
  }
  if (is.null(cell_labels)) cell_labels <- ct_comm$cell_labels
  label_names <- names(cell_labels)
  cell_labels <- as.character(cell_labels)
  if (!is.null(label_names)) names(cell_labels) <- label_names
  if (is.null(names(cell_labels))) {
    if (length(cell_labels) != ncol(reo_mat)) stop("Unnamed cell_labels must have length ncol(reo_mat).")
    names(cell_labels) <- colnames(reo_mat)
  }
  missing_cells <- setdiff(colnames(reo_mat), names(cell_labels))
  if (length(missing_cells) > 0) {
    shared_cells <- intersect(colnames(reo_mat), names(cell_labels))
    if (!length(shared_cells)) stop("cell_labels and reo_mat do not share any cells.")
    if (isTRUE(verbose)) {
      message(
        "Filtered ", length(missing_cells),
        " cells from permutation inputs because they are not present in cell_labels; first affected cells: ",
        paste(utils::head(missing_cells, 5), collapse = ", "),
        "."
      )
    }
    reo_mat <- reo_mat[, shared_cells, drop = FALSE]
    cell_labels <- cell_labels[shared_cells]
  } else {
    cell_labels <- cell_labels[colnames(reo_mat)]
  }
  invalid_labels <- is.na(cell_labels) | !nzchar(cell_labels)
  if (any(invalid_labels)) {
    n_filtered <- sum(invalid_labels)
    filtered_cells <- names(cell_labels)[invalid_labels]
    if (n_filtered == length(cell_labels)) {
      stop("All cell_labels are missing or empty; assign valid cell type labels before permutation.")
    }
    if (isTRUE(verbose)) {
      message(
        "Filtered ", n_filtered,
        " cells with missing or empty cell type labels before permutation; first affected cells: ",
        paste(utils::head(filtered_cells, 5), collapse = ", "),
        "."
      )
    }
    cell_labels <- cell_labels[!invalid_labels]
    reo_mat <- reo_mat[, names(cell_labels), drop = FALSE]
  }
  if (!is.null(knn_mat) && !is.null(rownames(knn_mat)) && !is.null(colnames(knn_mat)) &&
      all(names(cell_labels) %in% rownames(knn_mat)) && all(names(cell_labels) %in% colnames(knn_mat))) {
    knn_mat <- knn_mat[names(cell_labels), names(cell_labels), drop = FALSE]
  }
  if (is.null(lr_db)) lr_db <- ct_comm$lr_db
  observed <- ct_comm$pair_summary
  if (!metric %in% names(observed)) stop("metric not found in pair_summary: ", metric)
  feature <- paste(observed$sender_type, observed$receiver_type, sep = "|")
  null_mat <- matrix(NA_real_, nrow = nrow(observed), ncol = n_perm,
                     dimnames = list(feature, paste0("perm", seq_len(n_perm))))
  run_perm <- function() {
    perm_labels <- sample(as.character(cell_labels))
    names(perm_labels) <- names(cell_labels)
    perm <- do.call(
      summarize_celltype_communication,
      c(list(reo_mat = reo_mat, cell_labels = perm_labels, knn_mat = knn_mat,
             lr_db = lr_db, verbose = FALSE), dots)
    )
    key_perm <- paste(perm$pair_summary$sender_type, perm$pair_summary$receiver_type, sep = "|")
    vals <- stats::setNames(perm$pair_summary[[metric]], key_perm)
    vals[feature]
  }
  null_start <- proc.time()[["elapsed"]]
  null_step <- max(1, floor(n_perm / 10))
  for (b in seq_len(n_perm)) {
    null_mat[, b] <- run_perm()
    if (isTRUE(verbose) && (b == 1 || b %% null_step == 0 || b == n_perm)) {
      elapsed_min <- (proc.time()[["elapsed"]] - null_start) / 60
      message(sprintf("[CellTypeComm null] %d/%d permutations | elapsed %.1f min", b, n_perm, elapsed_min))
    }
  }
  obs <- observed[[metric]]
  if (isTRUE(adaptive) && is.finite(adaptive_n_perm) && adaptive_n_perm > n_perm) {
    prelim_p <- vapply(seq_along(obs), function(i) {
      if (!is.finite(obs[i])) return(NA_real_)
      (1 + sum(null_mat[i, ] >= obs[i], na.rm = TRUE)) / (1 + sum(!is.na(null_mat[i, ])))
    }, numeric(1))
    ord <- order(prelim_p, -obs, na.last = NA)
    top_idx <- ord[seq_len(min(adaptive_top_n, length(ord)))]
    extra_n <- as.integer(adaptive_n_perm - n_perm)
    if (length(top_idx) > 0 && extra_n > 0) {
      if (isTRUE(verbose)) {
        message(
          sprintf(
            "[CellTypeComm null adaptive] Starting adaptive refinement: %d additional full permutations for top %d candidates; target n_perm=%d. This is about %.1fx the preliminary run.",
            extra_n, length(top_idx), adaptive_n_perm, extra_n / max(1, n_perm)
          )
        )
      }
      extra_mat <- matrix(NA_real_, nrow = nrow(observed), ncol = extra_n,
                          dimnames = list(feature, paste0("adaptive", seq_len(extra_n))))
      adaptive_start <- proc.time()[["elapsed"]]
      adaptive_step <- max(1, floor(extra_n / 20))
      for (b in seq_len(extra_n)) {
        vals <- run_perm()
        extra_mat[top_idx, b] <- vals[top_idx]
        if (isTRUE(verbose) && (b == 1 || b %% adaptive_step == 0 || b == extra_n)) {
          elapsed_min <- (proc.time()[["elapsed"]] - adaptive_start) / 60
          message(sprintf("[CellTypeComm null adaptive] %d/%d additional permutations for top %d candidates | elapsed %.1f min", b, extra_n, length(top_idx), elapsed_min))
        }
      }
      null_mat <- cbind(null_mat, extra_mat)
    }
  }
  null_mean <- rowMeans(null_mat, na.rm = TRUE)
  null_sd <- apply(null_mat, 1, stats::sd, na.rm = TRUE)
  n_null_nonmissing <- rowSums(!is.na(null_mat))
  min_possible_p <- ifelse(n_null_nonmissing > 0, 1 / (1 + n_null_nonmissing), NA_real_)
  empirical_p <- vapply(seq_along(obs), function(i) {
    if (!is.finite(obs[i])) return(NA_real_)
    (1 + sum(null_mat[i, ] >= obs[i], na.rm = TRUE)) / (1 + sum(!is.na(null_mat[i, ])))
  }, numeric(1))
  z_score <- ifelse(null_sd > 0, (obs - null_mean) / null_sd,
                    ifelse(is.finite(obs) & is.finite(null_mean) & obs > null_mean, Inf,
                           ifelse(is.finite(obs) & is.finite(null_mean) & obs < null_mean, -Inf, NA_real_)))
  degenerate_null <- is.finite(null_sd) & null_sd == 0 & is.finite(obs) & is.finite(null_mean)
  degenerate_positive_null <- degenerate_null & obs > null_mean
  degenerate_negative_null <- degenerate_null & obs < null_mean
  fdr <- stats::p.adjust(empirical_p, method = "BH")
  min_possible_bh_fdr <- stats::p.adjust(min_possible_p, method = "BH")
  null_interpretation_flag <- ifelse(degenerate_positive_null, "degenerate_positive_null",
                                     ifelse(n_null_nonmissing < min_n_perm_publication, "low_resolution_null", "resolved_null"))
  null_biological_interpretation <- ifelse(
    degenerate_positive_null,
    "Permutation null has zero variance below the observation; treat as structural-null evidence requiring biological validation.",
    ifelse(n_null_nonmissing < min_n_perm_publication,
           "Permutation count limits p-value resolution; increase n_perm for publication claims.",
           "Permutation resolution is adequate for screening-level interpretation.")
  )

  data.frame(sender_type = observed$sender_type,
             receiver_type = observed$receiver_type,
             observed = obs,
             null_mean = null_mean,
             null_sd = null_sd,
             z_score = z_score,
             empirical_p = empirical_p,
             fdr = fdr,
             n_perm = as.integer(n_null_nonmissing),
             n_null_nonmissing = as.integer(n_null_nonmissing),
             min_possible_p = min_possible_p,
             min_possible_bh_fdr = min_possible_bh_fdr,
             degenerate_null = degenerate_null,
             degenerate_positive_null = degenerate_positive_null,
             degenerate_negative_null = degenerate_negative_null,
             null_interpretation_flag = null_interpretation_flag,
             null_biological_interpretation = null_biological_interpretation,
             metric = metric,
             stringsAsFactors = FALSE)
}

#' Diagnose Permutation Null Resolution
#'
#' @param null_pair Output from \code{permute_celltype_communication()}.
#' @return A \code{LogicCommPermutationDiagnostic} list.
#' @export
diagnose_permutation_resolution <- function(null_pair) {
  stopifnot(is.data.frame(null_pair))
  out <- list(
    n_tests = nrow(null_pair),
    n_low_resolution = sum(null_pair$null_interpretation_flag == "low_resolution_null", na.rm = TRUE),
    n_degenerate_null = sum(null_pair$degenerate_null %in% TRUE, na.rm = TRUE),
    n_degenerate_positive_null = sum(null_pair$degenerate_positive_null %in% TRUE, na.rm = TRUE),
    min_n_perm = if ("n_null_nonmissing" %in% names(null_pair)) min(null_pair$n_null_nonmissing, na.rm = TRUE) else NA_integer_,
    recommendations = c(
      "Use permutation p-values as structural evidence, not proof of causality.",
      "Increase n_perm for low-resolution nulls before publication-level claims.",
      "Inspect edge support and specificity for degenerate positive nulls."
    ),
    input = null_pair
  )
  class(out) <- "LogicCommPermutationDiagnostic"
  out
}

#' @export
print.LogicCommPermutationDiagnostic <- function(x, ...) {
  cat(sprintf("LogicCommPermutationDiagnostic | %d tests | %d degenerate positive nulls\n",
              x$n_tests, x$n_degenerate_positive_null))
  if (length(x$recommendations)) {
    cat("Recommendations:\n")
    cat(paste0("- ", x$recommendations, collapse = "\n"), "\n")
  }
  invisible(x)
}

#' Sensitivity Analysis Across REO Rank Thresholds
#'
#' @param expr_mat Seurat object or matrix-like input accepted by \code{calc_REO_matrix()}.
#' @param rank_threshold_grid Numeric thresholds to evaluate.
#' @param lr_db LR database.
#' @param cell_labels Cell-type labels.
#' @param knn_mat Optional graph.
#' @param layer Expression layer/slot.
#' @param level Feature level passed to \code{celltype_comm_to_lcs()}.
#' @param metric Metric passed to \code{celltype_comm_to_lcs()}.
#' @param verbose Print progress.
#' @param ... Passed to \code{summarize_celltype_communication()}.
#' @return A list with threshold-specific results and a stability table.
#' @export
sensitivity_REO_threshold <- function(expr_mat,
                                      rank_threshold_grid = c(0.4, 0.5, 0.6),
                                      lr_db = lr_pairs_human,
                                      cell_labels,
                                      knn_mat = NULL,
                                      layer = "counts",
                                      level = "celltype_lr",
                                      metric = "lcs",
                                      verbose = TRUE,
                                      ...) {
  lr_genes <- all_lr_genes(lr_db)
  results <- vector("list", length(rank_threshold_grid))
  names(results) <- paste0("thr_", rank_threshold_grid)
  vecs <- vector("list", length(rank_threshold_grid))
  names(vecs) <- names(results)
  for (i in seq_along(rank_threshold_grid)) {
    thr <- rank_threshold_grid[i]
    if (isTRUE(verbose)) message(sprintf("[Sensitivity] rank_threshold=%.3f", thr))
    reo <- calc_REO_matrix(expr_mat, lr_genes = lr_genes, rank_threshold = thr, layer = layer, verbose = FALSE)
    ct <- summarize_celltype_communication(reo, cell_labels = cell_labels, knn_mat = knn_mat,
                                           lr_db = lr_db, verbose = FALSE, ...)
    results[[i]] <- ct
    vecs[[i]] <- celltype_comm_to_lcs(ct, level = level, metric = metric, active_only = FALSE)
  }
  features <- sort(unique(unlist(lapply(vecs, names), use.names = FALSE)))
  mat <- matrix(0, nrow = length(features), ncol = length(vecs), dimnames = list(features, names(vecs)))
  for (i in seq_along(vecs)) mat[names(vecs[[i]]), i] <- vecs[[i]]
  stability <- data.frame(
    feature = rownames(mat),
    active_fraction = rowMeans(mat > 0, na.rm = TRUE),
    mean_lcs = rowMeans(mat, na.rm = TRUE),
    sd_lcs = apply(mat, 1, stats::sd, na.rm = TRUE),
    n_thresholds = ncol(mat),
    stringsAsFactors = FALSE
  )
  stability <- stability[order(stability$active_fraction, stability$mean_lcs, decreasing = TRUE), , drop = FALSE]
  out <- list(stability = stability, thresholds = rank_threshold_grid, results = results, matrix = mat)
  class(out) <- "LogicCommSensitivity"
  out
}
