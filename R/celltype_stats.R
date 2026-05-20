# R/celltype_stats.R

#' Bootstrap LogicComm Cell-Type Communication Scores
#'
#' Performs a lightweight binomial edge bootstrap from the observed LCS and edge
#' counts. This estimates score uncertainty without rerunning REO scanning.
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param n_boot Number of bootstrap replicates.
#' @param level \code{"celltype_lr"} or \code{"celltype_pair"}.
#' @param seed Optional random seed.
#' @return Data frame with observed score, bootstrap mean/sd/95 percent interval.
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
  boot_mat <- matrix(NA_real_, nrow = nrow(df), ncol = n_boot)
  for (b in seq_len(n_boot)) {
    boot_mat[, b] <- stats::rbinom(nrow(df), size = pmax(0, as.integer(round(df$n_edges))), prob = probs) / df$n_edges
  }
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
#' Recomputes cell-type communication after shuffling cell labels within a sample.
#' This tests whether observed sender-receiver structure exceeds what is expected
#' from the same graph and REO states but random cell-type assignments.
#'
#' @param ct_comm Optional observed object. If NULL, it is computed from inputs.
#' @param reo_mat Binary REO matrix or LogicCommREOResult.
#' @param cell_labels Named cell-type labels.
#' @param knn_mat Optional KNN/SNN graph.
#' @param lr_db LR database. Defaults to \code{ct_comm$lr_db} when available.
#' @param n_perm Number of permutations.
#' @param metric Pair-level metric to test.
#' @param seed Optional random seed.
#' @param adaptive If TRUE, refine the top preliminary candidates with additional permutations.
#' @param adaptive_top_n Number of preliminary candidates to refine when adaptive is TRUE.
#' @param adaptive_n_perm Target number of permutations for adaptive candidates.
#' @param min_n_perm_publication Recommended minimum permutation count for publication-level p-value resolution.
#' @param verbose Print progress.
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
  if (is.list(reo_mat) && !is.null(reo_mat$logic)) reo_mat <- reo_mat$logic
  if (is.null(ct_comm)) {
    ct_comm <- summarize_celltype_communication(reo_mat, cell_labels = cell_labels,
                                                knn_mat = knn_mat, lr_db = lr_db,
                                                verbose = FALSE, ...)
  }
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  if (is.null(cell_labels)) cell_labels <- ct_comm$cell_labels
  if (is.null(lr_db)) lr_db <- ct_comm$lr_db
  observed <- ct_comm$pair_summary
  if (!metric %in% names(observed)) stop("metric not found in pair_summary: ", metric)
  feature <- paste(observed$sender_type, observed$receiver_type, sep = "|")
  null_mat <- matrix(NA_real_, nrow = nrow(observed), ncol = n_perm,
                     dimnames = list(feature, paste0("perm", seq_len(n_perm))))
  for (b in seq_len(n_perm)) {
    perm_labels <- sample(as.character(cell_labels))
    names(perm_labels) <- names(cell_labels)
    perm <- summarize_celltype_communication(reo_mat, cell_labels = perm_labels,
                                             knn_mat = knn_mat, lr_db = lr_db,
                                             verbose = FALSE, ...)
    key_perm <- paste(perm$pair_summary$sender_type, perm$pair_summary$receiver_type, sep = "|")
    vals <- stats::setNames(perm$pair_summary[[metric]], key_perm)
    null_mat[, b] <- vals[feature]
    if (isTRUE(verbose) && (b %% max(1, floor(n_perm / 10)) == 0)) {
      message(sprintf("[CellTypeComm null] %d/%d permutations", b, n_perm))
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
      extra_mat <- matrix(NA_real_, nrow = nrow(observed), ncol = extra_n,
                          dimnames = list(feature, paste0("adaptive", seq_len(extra_n))))
      for (b in seq_len(extra_n)) {
        perm_labels <- sample(as.character(cell_labels))
        names(perm_labels) <- names(cell_labels)
        perm <- summarize_celltype_communication(reo_mat, cell_labels = perm_labels,
                                                 knn_mat = knn_mat, lr_db = lr_db,
                                                 verbose = FALSE, ...)
        key_perm <- paste(perm$pair_summary$sender_type, perm$pair_summary$receiver_type, sep = "|")
        vals <- stats::setNames(perm$pair_summary[[metric]], key_perm)
        extra_mat[top_idx, b] <- vals[feature[top_idx]]
        if (isTRUE(verbose) && (b %% max(1, floor(extra_n / 10)) == 0)) {
          message(sprintf("[CellTypeComm null adaptive] %d/%d additional permutations for top %d candidates", b, extra_n, length(top_idx)))
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
  
  out <- data.frame(sender_type = observed$sender_type,
                    receiver_type = observed$receiver_type,
                    observed = obs,
                    null_mean = null_mean,
                    null_sd = null_sd,
                    z_score = z_score,
                    empirical_p = empirical_p,
                    fdr = fdr,
                    n_perm = as.integer(n_null_nonmissing),
                    stringsAsFactors = FALSE)
  out
}
