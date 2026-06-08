# R/api_validation.R

.logic_require_ct_comm <- function(ct_comm) {
  if (!inherits(ct_comm, "LogicCommCellTypeComm")) {
    stop("ct_comm must be a LogicCommCellTypeComm object.")
  }
  if (is.null(ct_comm$lr_table) || !is.data.frame(ct_comm$lr_table)) {
    stop("ct_comm must contain an lr_table data frame.")
  }
  invisible(TRUE)
}

.logic_top_names <- function(names, scores, n = 5) {
  if (!length(names)) return(NA_character_)
  names <- as.character(names)
  scores <- as.numeric(scores)
  ok <- nzchar(names) & !is.na(names)
  names <- names[ok]
  scores <- scores[ok]
  if (!length(names)) return(NA_character_)
  ord <- order(scores, decreasing = TRUE, na.last = TRUE)
  paste(unique(names[ord])[seq_len(min(n, length(unique(names))))], collapse = ";")
}

.logic_mean_or_na <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

.logic_max_or_na <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

.logic_min_or_na <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

#' Validate receiver read-receipt response for LogicComm events
#'
#' @param ct_comm A \code{LogicCommCellTypeComm} object.
#' @param reo_mat REO matrix or \code{LogicCommREOResult}.
#' @param response_db Data frame with \code{lr_pair} and \code{response_genes}.
#' @param response_mode Response gene logic mode.
#' @param attach If \code{TRUE}, return the updated \code{ct_comm}; otherwise return
#'   the validation table.
#' @return Updated \code{LogicCommCellTypeComm} or validation data frame.
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
#' fun <- get("logic_validate_receiver_response")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
logic_validate_receiver_response <- function(ct_comm,
                                             reo_mat,
                                             response_db,
                                             response_mode = c("any", "all"),
                                             attach = TRUE) {
  .logic_require_ct_comm(ct_comm)
  response_mode <- match.arg(response_mode)
  out <- add_receiver_response_score(ct_comm, reo_mat, response_db, response_mode = response_mode)
  lr <- out$lr_table

  map <- c(
    response_gene_count = "validation_response_gene_count",
    response_active_frac = "validation_response_active_frac",
    receiver_response_score = "validation_response_score",
    response_integrated_score = "validation_response_integrated_lcs"
  )
  for (old in names(map)) {
    if (old %in% names(lr)) lr[[map[[old]]]] <- lr[[old]]
  }
  out$lr_table <- lr

  validation_cols <- intersect(
    c("sender_type", "receiver_type", "lr_pair", unname(map)),
    names(lr)
  )
  if (is.null(out$validation)) out$validation <- list()
  out$validation$receiver_response <- lr[, validation_cols, drop = FALSE]

  if (isTRUE(attach)) out else out$validation$receiver_response
}

#' Attach NicheNet-style ligand target evidence to LogicComm events
#'
#' @param ct_comm A \code{LogicCommCellTypeComm} object.
#' @param target_evidence Data frame with target evidence. Required columns are
#'   \code{target_gene}, \code{target_score}, and either \code{lr_pair} or
#'   \code{ligand}, depending on \code{by}. Optional \code{sender_type} and
#'   \code{receiver_type} columns restrict matches to specific cell-type pairs.
#' @param by Match target evidence by \code{"lr_pair"} or \code{"ligand"}.
#' @param attach If \code{TRUE}, return updated \code{ct_comm}; otherwise return
#'   the validation table.
#' @return Updated \code{LogicCommCellTypeComm} or validation data frame.
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
#' fun <- get("logic_validate_targets")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
logic_validate_targets <- function(ct_comm,
                                   target_evidence,
                                   by = c("lr_pair", "ligand"),
                                   attach = TRUE) {
  .logic_require_ct_comm(ct_comm)
  by <- match.arg(by)
  target_evidence <- as.data.frame(target_evidence, stringsAsFactors = FALSE)
  required <- c(by, "target_gene", "target_score")
  missing <- setdiff(required, names(target_evidence))
  if (length(missing) > 0) {
    stop("target_evidence is missing required columns: ", paste(missing, collapse = ", "))
  }
  target_evidence$target_score <- as.numeric(target_evidence$target_score)

  lr <- ct_comm$lr_table
  if (!by %in% names(lr)) stop("ct_comm$lr_table does not contain column '", by, "'.")

  summary <- data.frame(
    sender_type = lr$sender_type,
    receiver_type = lr$receiver_type,
    lr_pair = lr$lr_pair,
    validation_target_n = integer(nrow(lr)),
    validation_target_score_mean = NA_real_,
    validation_target_score_max = NA_real_,
    validation_target_active_frac = NA_real_,
    validation_target_top_genes = NA_character_,
    validation_target_source = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(lr))) {
    idx <- target_evidence[[by]] == lr[[by]][i]
    if ("sender_type" %in% names(target_evidence)) idx <- idx & target_evidence$sender_type == lr$sender_type[i]
    if ("receiver_type" %in% names(target_evidence)) idx <- idx & target_evidence$receiver_type == lr$receiver_type[i]
    ev <- target_evidence[idx %in% TRUE, , drop = FALSE]
    if (!nrow(ev)) next

    summary$validation_target_n[i] <- length(unique(ev$target_gene[nzchar(ev$target_gene) & !is.na(ev$target_gene)]))
    summary$validation_target_score_mean[i] <- .logic_mean_or_na(ev$target_score)
    summary$validation_target_score_max[i] <- .logic_max_or_na(ev$target_score)
    if ("target_active_frac" %in% names(ev)) {
      summary$validation_target_active_frac[i] <- .logic_mean_or_na(ev$target_active_frac)
    }
    summary$validation_target_top_genes[i] <- .logic_top_names(ev$target_gene, ev$target_score)
    if ("evidence_source" %in% names(ev)) {
      src <- unique(as.character(ev$evidence_source[nzchar(ev$evidence_source) & !is.na(ev$evidence_source)]))
      if (length(src)) summary$validation_target_source[i] <- paste(src, collapse = ";")
    }
  }

  if (!isTRUE(attach)) return(summary)
  for (col in setdiff(names(summary), c("sender_type", "receiver_type", "lr_pair"))) {
    ct_comm$lr_table[[col]] <- summary[[col]]
  }
  if (is.null(ct_comm$validation)) ct_comm$validation <- list()
  ct_comm$validation$ligand_targets <- summary
  ct_comm
}

#' Attach TF-switch validation evidence to LogicComm events
#'
#' @param ct_comm A \code{LogicCommCellTypeComm} object.
#' @param tf_activity Data frame with receiver TF activity. Required columns are
#'   \code{receiver_type}, \code{tf}, and \code{tf_activity_score}.
#' @param lr_tf_map Data frame mapping LR events, pathways, or ligands to TFs.
#'   Must contain \code{tf} and at least one of \code{lr_pair}, \code{pathway},
#'   or \code{ligand}.
#' @param attach If \code{TRUE}, return updated \code{ct_comm}; otherwise return
#'   the validation table.
#' @return Updated \code{LogicCommCellTypeComm} or validation data frame.
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
#' fun <- get("logic_validate_tf_switch")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
logic_validate_tf_switch <- function(ct_comm,
                                     tf_activity,
                                     lr_tf_map,
                                     attach = TRUE) {
  .logic_require_ct_comm(ct_comm)
  tf_activity <- as.data.frame(tf_activity, stringsAsFactors = FALSE)
  lr_tf_map <- as.data.frame(lr_tf_map, stringsAsFactors = FALSE)

  missing_activity <- setdiff(c("receiver_type", "tf", "tf_activity_score"), names(tf_activity))
  if (length(missing_activity) > 0) {
    stop("tf_activity is missing required columns: ", paste(missing_activity, collapse = ", "))
  }
  if (!"tf" %in% names(lr_tf_map)) stop("lr_tf_map must contain a 'tf' column.")

  lr <- ct_comm$lr_table
  key <- intersect(c("lr_pair", "pathway", "ligand"), intersect(names(lr), names(lr_tf_map)))[1]
  if (is.na(key)) stop("lr_tf_map must share at least one key with lr_table: lr_pair, pathway, or ligand.")
  tf_activity$tf_activity_score <- as.numeric(tf_activity$tf_activity_score)

  summary <- data.frame(
    sender_type = lr$sender_type,
    receiver_type = lr$receiver_type,
    lr_pair = lr$lr_pair,
    validation_tf_n = integer(nrow(lr)),
    validation_tf_score_mean = NA_real_,
    validation_tf_score_max = NA_real_,
    validation_tf_switch_delta = NA_real_,
    validation_tf_switch_fdr_min = NA_real_,
    validation_tf_top = NA_character_,
    validation_tf_source = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(lr))) {
    tf_set <- unique(as.character(lr_tf_map$tf[lr_tf_map[[key]] == lr[[key]][i]]))
    tf_set <- tf_set[nzchar(tf_set) & !is.na(tf_set)]
    if (!length(tf_set)) next
    ev <- tf_activity[tf_activity$receiver_type == lr$receiver_type[i] & tf_activity$tf %in% tf_set, , drop = FALSE]
    if (!nrow(ev)) next

    summary$validation_tf_n[i] <- length(unique(ev$tf))
    summary$validation_tf_score_mean[i] <- .logic_mean_or_na(ev$tf_activity_score)
    summary$validation_tf_score_max[i] <- .logic_max_or_na(ev$tf_activity_score)
    if ("tf_activity_delta" %in% names(ev)) {
      summary$validation_tf_switch_delta[i] <- .logic_mean_or_na(ev$tf_activity_delta)
    }
    if ("tf_activity_fdr" %in% names(ev)) {
      summary$validation_tf_switch_fdr_min[i] <- .logic_min_or_na(ev$tf_activity_fdr)
    }
    summary$validation_tf_top[i] <- .logic_top_names(ev$tf, ev$tf_activity_score)
    if ("evidence_source" %in% names(ev)) {
      src <- unique(as.character(ev$evidence_source[nzchar(ev$evidence_source) & !is.na(ev$evidence_source)]))
      if (length(src)) summary$validation_tf_source[i] <- paste(src, collapse = ";")
    }
  }

  if (!isTRUE(attach)) return(summary)
  for (col in setdiff(names(summary), c("sender_type", "receiver_type", "lr_pair"))) {
    ct_comm$lr_table[[col]] <- summary[[col]]
  }
  if (is.null(ct_comm$validation)) ct_comm$validation <- list()
  ct_comm$validation$tf_switches <- summary
  ct_comm
}

#' Attach external validation evidence to LogicComm events
#'
#' @param ct_comm A \code{LogicCommCellTypeComm} object.
#' @param evidence_table Data frame with optional identity columns such as
#'   \code{sender_type}, \code{receiver_type}, \code{lr_pair}, \code{ligand},
#'   \code{receptor}, or \code{pathway}, and optional \code{evidence_score},
#'   \code{evidence_label}, \code{evidence_source}, and \code{evidence_type}.
#' @param evidence_type Optional evidence type to assign when not present in
#'   \code{evidence_table}.
#' @param attach If \code{TRUE}, return updated \code{ct_comm}; otherwise return
#'   the validation table.
#' @return Updated \code{LogicCommCellTypeComm} or validation data frame.
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
#' fun <- get("logic_add_external_evidence")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
logic_add_external_evidence <- function(ct_comm,
                                        evidence_table,
                                        evidence_type = NULL,
                                        attach = TRUE) {
  .logic_require_ct_comm(ct_comm)
  evidence_table <- as.data.frame(evidence_table, stringsAsFactors = FALSE)
  if (!is.null(evidence_type) && !"evidence_type" %in% names(evidence_table)) {
    evidence_table$evidence_type <- evidence_type
  }
  if ("evidence_score" %in% names(evidence_table)) {
    evidence_table$evidence_score <- as.numeric(evidence_table$evidence_score)
  }

  lr <- ct_comm$lr_table
  key_cols <- intersect(c("sender_type", "receiver_type", "lr_pair", "ligand", "receptor", "pathway"),
                        intersect(names(lr), names(evidence_table)))
  if (!length(key_cols)) stop("evidence_table must share at least one identity column with ct_comm$lr_table.")

  summary <- data.frame(
    sender_type = lr$sender_type,
    receiver_type = lr$receiver_type,
    lr_pair = lr$lr_pair,
    validation_external_n = integer(nrow(lr)),
    validation_external_score_mean = NA_real_,
    validation_external_label = NA_character_,
    validation_external_source = NA_character_,
    validation_external_type = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(lr))) {
    idx <- rep(TRUE, nrow(evidence_table))
    for (key in key_cols) idx <- idx & evidence_table[[key]] == lr[[key]][i]
    ev <- evidence_table[idx %in% TRUE, , drop = FALSE]
    if (!nrow(ev)) next

    summary$validation_external_n[i] <- nrow(ev)
    if ("evidence_score" %in% names(ev)) {
      summary$validation_external_score_mean[i] <- .logic_mean_or_na(ev$evidence_score)
    }
    if ("evidence_label" %in% names(ev)) {
      labs <- unique(as.character(ev$evidence_label[nzchar(ev$evidence_label) & !is.na(ev$evidence_label)]))
      if (length(labs)) summary$validation_external_label[i] <- paste(labs, collapse = ";")
    }
    if ("evidence_source" %in% names(ev)) {
      src <- unique(as.character(ev$evidence_source[nzchar(ev$evidence_source) & !is.na(ev$evidence_source)]))
      if (length(src)) summary$validation_external_source[i] <- paste(src, collapse = ";")
    }
    if ("evidence_type" %in% names(ev)) {
      typ <- unique(as.character(ev$evidence_type[nzchar(ev$evidence_type) & !is.na(ev$evidence_type)]))
      if (length(typ)) summary$validation_external_type[i] <- paste(typ, collapse = ";")
    }
  }

  if (!isTRUE(attach)) return(summary)
  for (col in setdiff(names(summary), c("sender_type", "receiver_type", "lr_pair"))) {
    ct_comm$lr_table[[col]] <- summary[[col]]
  }
  if (is.null(ct_comm$validation)) ct_comm$validation <- list()
  ct_comm$validation$external <- summary
  ct_comm
}

#' Grade LogicComm communication evidence
#'
#' @param ct_comm A \code{LogicCommCellTypeComm} object.
#' @param use_validation Include attached \code{validation_*} evidence columns when
#'   assigning grades.
#' @param min_edges Optional minimum edge threshold for caution text.
#' @param min_active_edges Optional minimum active-edge threshold for caution text.
#' @return Updated \code{LogicCommCellTypeComm} with evidence-grade columns added
#'   to \code{lr_table}.
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
#' fun <- get("logic_grade_evidence")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
logic_grade_evidence <- function(ct_comm,
                                 use_validation = TRUE,
                                 min_edges = NULL,
                                 min_active_edges = NULL) {
  .logic_require_ct_comm(ct_comm)
  lr <- ct_comm$lr_table
  n <- nrow(lr)

  has_positive <- function(col) {
    col %in% names(lr) & !is.na(lr[[col]]) & as.numeric(lr[[col]]) > 0
  }

  local <- if ("local_active" %in% names(lr)) lr$local_active %in% TRUE else rep(FALSE, n)
  global <- if ("global_candidate_active" %in% names(lr)) lr$global_candidate_active %in% TRUE else rep(FALSE, n)
  distal <- if ("distal_candidate" %in% names(lr)) lr$distal_candidate %in% TRUE else rep(FALSE, n)
  active <- if ("active" %in% names(lr)) lr$active %in% TRUE else local | global

  validation <- rep(FALSE, n)
  if (isTRUE(use_validation)) {
    for (col in c("validation_response_score", "validation_target_n", "validation_tf_n", "validation_external_n")) {
      if (col %in% names(lr)) validation <- validation | has_positive(col)
    }
  }

  grade <- rep(NA_character_, n)
  tier <- rep(NA_integer_, n)
  grade[local & validation] <- "A"
  tier[local & validation] <- 5L
  grade[local & is.na(grade)] <- "B"
  tier[local & is.na(tier)] <- 4L
  grade[(distal | global) & is.na(grade)] <- "C"
  tier[(distal | global) & is.na(tier)] <- 3L
  grade[active & is.na(grade)] <- "D"
  tier[active & is.na(tier)] <- 2L

  components <- character(n)
  reason <- character(n)
  caution <- character(n)
  recommended <- character(n)

  for (i in seq_len(n)) {
    comp <- character(0)
    if (local[i]) comp <- c(comp, "local_neighborhood")
    if (global[i]) comp <- c(comp, "global_candidate")
    if (distal[i]) comp <- c(comp, "distal_candidate")
    if (validation[i]) comp <- c(comp, "validation")
    if (!length(comp)) comp <- "insufficient_evidence"
    components[i] <- paste(comp, collapse = ";")

    if (identical(grade[i], "A")) {
      reason[i] <- "Local graph-supported LR logic with downstream or external validation evidence."
      caution[i] <- "Strong candidate; still requires biological context and orthogonal validation for causal claims."
      recommended[i] <- "Prioritize for mechanistic follow-up."
    } else if (identical(grade[i], "B")) {
      reason[i] <- "Local graph-supported LR logic with sufficient core LogicComm evidence."
      caution[i] <- "Communication candidate; add receiver target, TF, protein, spatial, or perturbation evidence when possible."
      recommended[i] <- "Validate receiver response or external support."
    } else if (identical(grade[i], "C")) {
      reason[i] <- "Global or distal candidate potential without strong local evidence."
      caution[i] <- "Interpret as distal/global potential rather than local communication proof."
      recommended[i] <- "Check spatial distance, ligand class, receiver response, and sample recurrence."
    } else if (identical(grade[i], "D")) {
      reason[i] <- "Weak or limited communication evidence."
      caution[i] <- "Low-priority candidate unless supported by external biology."
      recommended[i] <- "Inspect edge counts and repeat under sensitivity analyses."
    } else {
      reason[i] <- "Insufficient evidence for an active communication call."
      caution[i] <- "Do not interpret as active communication."
      recommended[i] <- "No follow-up unless independently motivated."
    }
  }

  if (!is.null(min_edges) && "n_edges" %in% names(lr)) {
    low <- !is.na(lr$n_edges) & lr$n_edges < min_edges
    caution[low] <- paste(caution[low], "Edge opportunity count is below the requested threshold.")
  }
  if (!is.null(min_active_edges) && "n_active_edges" %in% names(lr)) {
    low <- !is.na(lr$n_active_edges) & lr$n_active_edges < min_active_edges
    caution[low] <- paste(caution[low], "Active edge support is below the requested threshold.")
  }

  ct_comm$lr_table$evidence_grade <- grade
  ct_comm$lr_table$evidence_tier <- tier
  ct_comm$lr_table$evidence_components <- components
  ct_comm$lr_table$evidence_reason <- reason
  ct_comm$lr_table$interpretation_caution <- caution
  ct_comm$lr_table$recommended_validation <- recommended
  ct_comm
}
