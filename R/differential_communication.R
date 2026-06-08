# R/differential_communication.R — Subgroup-resolved multi-sample differential communication

#' Differential Cell-Type-Resolved Communication Between Sample Groups
#'
#' Compares cell-type-resolved ligand-receptor communication between two sample
#' groups and organizes the result by the sender -> receiver subgroup on which
#' each difference occurs. This makes the practical question "which differential
#' L-R pairs change, and between which cell-type pairs?" directly answerable.
#'
#' @details
#' The biological sample is the unit of comparison. For each sample, the
#' cell-type-resolved feature vector (sender|receiver|L-R) is extracted with
#' \code{\link{celltype_comm_to_lcs}} and compared across groups with
#' \code{\link{CompareLogicGroups}}. The composite feature keys are then split
#' into explicit \code{sender_type}, \code{receiver_type}, and \code{lr_pair}
#' columns, and a per-subgroup summary counts how many L-R pairs change in each
#' direction on each sender -> receiver cell-type pair.
#'
#' \strong{On FDR and sample size.} Fisher's exact p-values are bounded by the
#' number of samples per group; with only a few replicates per group the
#' smallest attainable p-value is large, so Benjamini-Hochberg FDR saturates
#' near 1 even for fully separating effects. In that regime, prioritize by
#' effect size (\code{asymmetry}, \code{log2fc_lcs}) and corroborate with the
#' edge-count model \code{\link{fit_celltype_comm_glm}} and additional
#' biological replicates rather than reading the FDR as evidence of absence.
#'
#' @param sample_ct_list Named list of per-sample \code{LogicCommCellTypeComm}
#'   objects (one per biological sample), or a named list of pre-computed
#'   cell-type-resolved LCS vectors with \code{sender|receiver|lr_pair} names.
#' @param group_info Named vector mapping sample names to group labels.
#' @param case_label,ctrl_label Group labels. Defaults: \code{"Case"}, \code{"Ctrl"}.
#' @param metric Cell-type LR metric to compare. Default: \code{"lcs"}.
#' @param lcs_threshold Threshold for calling a sample positive in
#'   \code{\link{CompareLogicGroups}}. Default: \code{0.01}.
#' @param min_samples_per_group Minimum non-missing samples per group. Default: \code{1}.
#' @param min_asymmetry Minimum absolute Case-Ctrl frequency asymmetry for the
#'   per-subgroup change counts. Default: \code{0.2}.
#' @param lr_db Optional LR database used to attach pathway annotation. Defaults
#'   to the database stored in the supplied \code{ct_comm} objects.
#' @param verbose Print progress messages. Default: \code{TRUE}.
#'
#' @return A list of class \code{LogicCommDifferential} with \code{lr} (per
#'   sender/receiver/L-R differential table), \code{subgroup} (per sender ->
#'   receiver summary), and the underlying \code{comparison}.
#' @seealso \code{\link{plot_differential_communication_summary}},
#'   \code{\link{plot_differential_celltype_heatmap}}, \code{\link{fit_celltype_comm_glm}}
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
#' fun <- get("differential_celltype_communication")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
differential_celltype_communication <- function(sample_ct_list,
                                                group_info,
                                                case_label = "Case",
                                                ctrl_label = "Ctrl",
                                                metric = "lcs",
                                                lcs_threshold = 0.01,
                                                min_samples_per_group = 1,
                                                min_asymmetry = 0.2,
                                                lr_db = NULL,
                                                verbose = TRUE) {
  stopifnot(is.list(sample_ct_list), length(sample_ct_list) >= 2)
  is_ct <- vapply(sample_ct_list, function(x) inherits(x, "LogicCommCellTypeComm"), logical(1))
  if (any(is_ct) && !all(is_ct)) {
    stop("sample_ct_list must be either all LogicCommCellTypeComm objects or all named LCS vectors.")
  }
  if (all(is_ct)) {
    if (is.null(lr_db)) lr_db <- sample_ct_list[[which(is_ct)[1]]]$lr_db
    lcs_list <- lapply(sample_ct_list, celltype_comm_to_lcs, level = "celltype_lr", metric = metric)
  } else {
    lcs_list <- lapply(sample_ct_list, function(v) {
      if (!is.numeric(v)) stop("Non-ct_comm entries must be numeric named LCS vectors.")
      v
    })
  }

  cmp <- CompareLogicGroups(lcs_list, group_info = group_info, case_label = case_label,
                            ctrl_label = ctrl_label, lcs_threshold = lcs_threshold,
                            min_samples_per_group = min_samples_per_group, lr_db = NULL,
                            verbose = FALSE)

  keys <- strsplit(as.character(cmp$lr_pair), "|", fixed = TRUE)
  sender <- vapply(keys, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
  receiver <- vapply(keys, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  lrp <- vapply(keys, function(x) if (length(x) >= 3) paste(x[3:length(x)], collapse = "|") else NA_character_, character(1))

  lr <- data.frame(sender_type = sender, receiver_type = receiver, lr_pair = lrp,
                   stringsAsFactors = FALSE)
  carry <- intersect(c("case_freq", "ctrl_freq", "asymmetry", "case_mean_lcs", "ctrl_mean_lcs",
                       "log2fc_lcs", "p_fisher", "fdr_fisher", "p_wilcox", "fdr_wilcox",
                       "n_case_avail", "n_ctrl_avail"), names(cmp))
  for (cc in carry) lr[[cc]] <- cmp[[cc]]

  if (!is.null(lr_db) && all(c("lr_pair", "pathway") %in% names(lr_db))) {
    lr$pathway <- as.character(lr_db$pathway)[match(lr$lr_pair, lr_db$lr_pair)]
  } else {
    lr$pathway <- NA_character_
  }
  lr$direction <- ifelse(is.na(lr$asymmetry), NA_character_,
                  ifelse(lr$asymmetry > 0, paste0(case_label, "-up"),
                  ifelse(lr$asymmetry < 0, paste0(ctrl_label, "-up"), "unchanged")))
  lr <- lr[order(-abs(lr$asymmetry), -abs(lr$log2fc_lcs), na.last = TRUE), , drop = FALSE]
  rownames(lr) <- NULL

  subgroup <- .differential_subgroup_summary(lr, min_asymmetry, case_label, ctrl_label)

  res <- list(
    lr = lr,
    subgroup = subgroup,
    comparison = cmp,
    params = list(case_label = case_label, ctrl_label = ctrl_label, metric = metric,
                  min_asymmetry = min_asymmetry,
                  n_case = sum(group_info == case_label, na.rm = TRUE),
                  n_ctrl = sum(group_info == ctrl_label, na.rm = TRUE))
  )
  class(res) <- "LogicCommDifferential"
  if (isTRUE(verbose)) {
    message(sprintf("[Differential] %d cell-type-resolved L-R features across %d sender->receiver subgroups.",
                    nrow(lr), nrow(subgroup)))
  }
  res
}

#' Per sender -> receiver subgroup differential summary
#' @keywords internal
.differential_subgroup_summary <- function(lr, min_asymmetry, case_label, ctrl_label) {
  if (is.null(lr) || !nrow(lr)) return(data.frame())
  key <- paste(lr$sender_type, lr$receiver_type, sep = "|||")
  sp <- split(seq_len(nrow(lr)), key)
  out <- lapply(sp, function(ii) {
    sub <- lr[ii, , drop = FALSE]
    big <- !is.na(sub$asymmetry) & abs(sub$asymmetry) >= min_asymmetry
    up <- sum(big & sub$asymmetry > 0)
    down <- sum(big & sub$asymmetry < 0)
    top_i <- if (any(!is.na(sub$asymmetry))) which.max(abs(sub$asymmetry)) else integer(0)
    data.frame(
      sender_type = sub$sender_type[1],
      receiver_type = sub$receiver_type[1],
      n_changed = up + down,
      n_case_up = up,
      n_ctrl_up = down,
      mean_abs_asymmetry = mean(abs(sub$asymmetry), na.rm = TRUE),
      max_abs_asymmetry = if (any(!is.na(sub$asymmetry))) max(abs(sub$asymmetry), na.rm = TRUE) else NA_real_,
      top_lr_pair = if (length(top_i)) sub$lr_pair[top_i] else NA_character_,
      dominant_direction = if (up > down) paste0(case_label, "-up") else if (down > up) paste0(ctrl_label, "-up") else "balanced",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[order(-out$n_changed, -out$mean_abs_asymmetry), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Print method for LogicCommDifferential
#' @param x LogicCommDifferential object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}, called for side effects.
#' @examples
#' diff <- structure(
#'   list(
#'     lr = data.frame(lr_pair = "L_R", fdr_fisher = 0.5, asymmetry = 0.2),
#'     subgroup = data.frame(sender_type = "T", receiver_type = "B",
#'                           n_changed = 1, n_case_up = 1, n_ctrl_up = 0,
#'                           top_lr_pair = "L_R", dominant_direction = "case"),
#'     params = list(case_label = "Case", ctrl_label = "Ctrl", n_case = 2, n_ctrl = 2)
#'   ),
#'   class = "LogicCommDifferential"
#' )
#' print(diff)
#' @export
print.LogicCommDifferential <- function(x, ...) {
  cat(sprintf("LogicCommDifferential | %d L-R features | %d sender->receiver subgroups | %s vs %s\n",
              nrow(x$lr), nrow(x$subgroup), x$params$case_label, x$params$ctrl_label))
  if (nrow(x$subgroup)) {
    cat("Top differential subgroups (by number of changed L-R pairs):\n")
    show <- intersect(c("sender_type", "receiver_type", "n_changed", "n_case_up",
                        "n_ctrl_up", "top_lr_pair", "dominant_direction"), names(x$subgroup))
    print(utils::head(x$subgroup[, show, drop = FALSE], 6), row.names = FALSE)
  }
  saturated <- x$params$n_case < 3 || x$params$n_ctrl < 3 ||
    (("fdr_fisher" %in% names(x$lr)) && all(x$lr$fdr_fisher >= 0.99, na.rm = TRUE))
  if (isTRUE(saturated)) {
    cat(sprintf(
      "\nNote: with %d vs %d samples, Fisher FDR saturates near 1; rank by effect size (asymmetry / log2fc_lcs) and corroborate with fit_celltype_comm_glm() and more replicates.\n",
      x$params$n_case, x$params$n_ctrl))
  }
  invisible(x)
}

#' Plot Differential L-R Pairs Per Sender -> Receiver Subgroup
#'
#' Shows, for each sender -> receiver cell-type subgroup, how many ligand-receptor
#' pairs change in each direction, so the cell-type pairs that carry the
#' differential communication are immediately visible.
#'
#' @param diff A \code{LogicCommDifferential} object.
#' @param top_n Number of top subgroups (by number of changed pairs). Default: \code{20}.
#' @param title Optional plot title.
#' @return A ggplot2 object.
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
#' fun <- get("plot_differential_communication_summary")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
plot_differential_communication_summary <- function(diff, top_n = 20, title = NULL) {
  stopifnot(inherits(diff, "LogicCommDifferential"))
  sg <- diff$subgroup
  sg <- sg[sg$n_changed > 0, , drop = FALSE]
  if (!nrow(sg)) stop("No differential subgroups at min_asymmetry = ", diff$params$min_asymmetry, ".")
  sg <- utils::head(sg, top_n)
  sg$sender_receiver <- paste(sg$sender_type, "\u2192", sg$receiver_type)
  case_lab <- paste0(diff$params$case_label, "-up")
  ctrl_lab <- paste0(diff$params$ctrl_label, "-up")
  long <- rbind(
    data.frame(sender_receiver = sg$sender_receiver, count = sg$n_case_up,
               diff_direction = case_lab, stringsAsFactors = FALSE),
    data.frame(sender_receiver = sg$sender_receiver, count = sg$n_ctrl_up,
               diff_direction = ctrl_lab, stringsAsFactors = FALSE)
  )
  long$sender_receiver <- factor(long$sender_receiver, levels = rev(sg$sender_receiver))
  if (is.null(title)) title <- "Differential L-R pairs per sender -> receiver subgroup"
  ggplot2::ggplot(long, ggplot2::aes(x = count, y = sender_receiver, fill = diff_direction)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(c("#b2182b", "#2166ac"), c(case_lab, ctrl_lab)), name = NULL) +
    ggplot2::labs(title = title, x = "Number of differential L-R pairs", y = NULL) +
    theme_logiccomm()
}
