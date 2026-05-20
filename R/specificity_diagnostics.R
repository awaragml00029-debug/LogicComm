# R/specificity_diagnostics.R

.specificity_entropy <- function(values, levels) {
  x <- setNames(rep(0, length(levels)), levels)
  if (length(values) > 0) {
    vals <- tapply(values$value, values$key, sum, na.rm = TRUE)
    common <- intersect(names(vals), names(x))
    x[common] <- as.numeric(vals[common])
  }
  x[!is.finite(x) | x < 0] <- 0
  if (sum(x) <= 0 || length(x) <= 1) return(NA_real_)
  p <- x / sum(x)
  p <- p[p > 0]
  if (!length(p)) return(NA_real_)
  -sum(p * log(p)) / log(length(x))
}

.specificity_score <- function(values, levels) {
  h <- .specificity_entropy(values, levels)
  ifelse(is.na(h), NA_real_, 1 - h)
}

.first_nonempty <- function(x, default = NA_character_) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[1] else default
}

.weighted_mean_safe <- function(x, w = NULL) {
  ok <- is.finite(x)
  if (!is.null(w)) ok <- ok & is.finite(w) & w >= 0
  if (!any(ok)) return(NA_real_)
  if (is.null(w)) return(mean(x[ok], na.rm = TRUE))
  if (sum(w[ok], na.rm = TRUE) <= 0) return(mean(x[ok], na.rm = TRUE))
  stats::weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

.parse_celltype_feature_names <- function(feature) {
  parts <- strsplit(as.character(feature), "\\|", fixed = FALSE)
  out <- lapply(parts, function(x) {
    x <- c(x, rep(NA_character_, max(0, 3 - length(x))))
    data.frame(sender_type = x[1], receiver_type = x[2], feature = x[3], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Score Communication Specificity and Ubiquity
#'
#' Adds specificity annotations to a cell-type communication object. This is
#' useful for distinguishing strong but broad axes from more specific candidates.
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @param active_only If TRUE, compute specificity using active L-R rows only.
#' @param broad_pair_specificity_cutoff Pair specificity below this threshold is
#'   flagged as broad/ubiquitous.
#' @param ubiquitous_pair_fraction Fraction of possible sender-receiver pairs
#'   above which an LR pair is flagged as ubiquitous.
#' @param identity_pathway_patterns Regular expressions for identity pathways.
#' @param identity_lr_patterns Regular expressions for identity L-R names.
#' @param verbose Print a short summary.
#' @return Updated \code{LogicCommCellTypeComm} object.
#' @export
score_communication_specificity <- function(ct_comm,
                                            active_only = TRUE,
                                            broad_pair_specificity_cutoff = 0.25,
                                            ubiquitous_pair_fraction = 0.25,
                                            identity_pathway_patterns = c("^MHC", "ANTIGEN", "BCR", "TCR"),
                                            identity_lr_patterns = c("^HLA", "_CD4$", "_CD8", "B2M"),
                                            verbose = TRUE) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  lr <- ct_comm$lr_table
  if (!nrow(lr)) return(ct_comm)
  
  cell_types <- ct_comm$role_summary$cell_type
  if (is.null(cell_types)) cell_types <- sort(unique(c(lr$sender_type, lr$receiver_type)))
  
  pair_levels <- paste(ct_comm$pair_summary$sender_type, ct_comm$pair_summary$receiver_type, sep = "|||")
  
  use <- !is.na(lr$lcs) & is.finite(lr$lcs)
  if (isTRUE(active_only) && "active" %in% names(lr)) use <- use & lr$active %in% TRUE
  dat <- lr[use, , drop = FALSE]
  if (!nrow(dat)) {
    ct_comm$specificity_summary <- data.frame()
    return(ct_comm)
  }
  dat$pair_key <- paste(dat$sender_type, dat$receiver_type, sep = "|||")

  sp <- split(seq_len(nrow(dat)), dat$lr_pair)
  out <- lapply(sp, function(ii) {
    sub <- dat[ii, , drop = FALSE]
    sender_spec <- .specificity_score(data.frame(key = sub$sender_type, value = sub$lcs), cell_types)
    receiver_spec <- .specificity_score(data.frame(key = sub$receiver_type, value = sub$lcs), cell_types)
    pair_spec <- .specificity_score(data.frame(key = sub$pair_key, value = sub$lcs), pair_levels)
    
    active_pair_n <- length(unique(sub$pair_key))
    pair_fraction <- active_pair_n / max(1, length(pair_levels))
    
    pathway <- .first_nonempty(sub$pathway, "Unknown")
    lr_pair <- sub$lr_pair[1]
    
    identity_flag <- any(grepl(paste(identity_pathway_patterns, collapse = "|"), pathway, ignore.case = TRUE)) ||
      any(grepl(paste(identity_lr_patterns, collapse = "|"), lr_pair, ignore.case = TRUE))
    
    ubiquitous_flag <- isTRUE(pair_fraction >= ubiquitous_pair_fraction) ||
      (!is.na(pair_spec) && pair_spec <= broad_pair_specificity_cutoff)
    
    class <- if (identity_flag && ubiquitous_flag) {
      "Identity-associated broad axis"
    } else if (identity_flag) {
      "Identity-associated axis"
    } else if (ubiquitous_flag) {
      "Broad/ubiquitous axis"
    } else if (!is.na(pair_spec) && pair_spec >= 0.6) {
      "Pair-specific candidate"
    } else {
      "Intermediate specificity"
    }
    
    data.frame(
      lr_pair = lr_pair,
      pathway = pathway,
      active_celltype_pair_fraction = pair_fraction,
      total_lcs = sum(sub$lcs, na.rm = TRUE),
      pair_specificity = pair_spec,
      ubiquitous_interaction_flag = ubiquitous_flag,
      identity_associated_flag = identity_flag,
      specificity_class = class,
      stringsAsFactors = FALSE
    )
  })
  
  summary <- do.call(rbind, out)
  ct_comm$specificity_summary <- summary
  
  # Augment lr_table
  idx <- match(lr$lr_pair, summary$lr_pair)
  lr$specificity_class <- summary$specificity_class[idx]
  lr$pair_specificity <- summary$pair_specificity[idx]
  lr$identity_associated_flag <- summary$identity_associated_flag[idx]
  lr$ubiquitous_interaction_flag <- summary$ubiquitous_interaction_flag[idx]
  ct_comm$lr_table <- lr
  
  if (verbose) message(sprintf("[Specificity] Annotated %d LR pairs.", nrow(summary)))
  ct_comm
}

#' Rank Publication-Priority Communication Axes
#'
#' @param ct_comm Cell-type communication object.
#' @param null_pair Optional output from \code{permute_celltype_communication()}.
#' @param sens Optional output from \code{sensitivity_REO_threshold()}.
#' @param top_n Number of rows to return.
#' @return Data frame ranked by publication priority.
#' @export
rank_communication_axes <- function(ct_comm,
                                    null_pair = NULL,
                                    sens = NULL,
                                    top_n = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  if (is.null(ct_comm$specificity_summary)) {
    ct_comm <- score_communication_specificity(ct_comm, verbose = FALSE)
  }
  
  df <- ct_comm$lr_table
  df <- df[df$active %in% TRUE & is.finite(df$lcs), , drop = FALSE]
  if (!nrow(df)) return(df)
  
  df$publication_priority_score <- .rescale01(df$lcs) * 0.4 + 
    .rescale01(ifelse(is.na(df$pair_specificity), 0, df$pair_specificity)) * 0.4
    
  if (!is.null(null_pair)) {
    # Add permutation component
  }
  
  df <- df[order(df$publication_priority_score, decreasing = TRUE), , drop = FALSE]
  if (!is.null(top_n)) df <- head(df, top_n)
  df
}

#' Interpret Cell-Type Communication Roles
#'
#' @param ct_comm Output from \code{summarize_celltype_communication()}.
#' @return Data frame with role interpretation and caution labels.
#' @export
interpret_celltype_roles <- function(ct_comm) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  roles <- ct_comm$role_summary
  if (is.null(roles) || !nrow(roles)) return(data.frame())

  reliability <- roles$role_reliability_label
  if (is.null(reliability)) reliability <- roles$role_confidence_label
  if (is.null(reliability)) reliability <- rep(NA_character_, nrow(roles))

  separation <- roles$role_separation_label
  if (is.null(separation)) separation <- rep(NA_character_, nrow(roles))

  low_communication <- if ("low_communication" %in% names(roles)) roles$low_communication %in% TRUE else rep(FALSE, nrow(roles))
  evidence <- if ("communication_evidence_label" %in% names(roles)) roles$communication_evidence_label else rep(NA_character_, nrow(roles))
  dominant_role <- if ("dominant_role" %in% names(roles)) roles$dominant_role else roles$dominant_role_strict
  strict_role <- if ("dominant_role_strict" %in% names(roles)) roles$dominant_role_strict else dominant_role
  biological <- if ("role_biological_interpretation" %in% names(roles)) roles$role_biological_interpretation else rep(NA_character_, nrow(roles))

  caution <- ifelse(low_communication | reliability %in% c("Low", "Low-communication"),
                    "Low communication evidence; interpret role assignment cautiously.",
                    ifelse(separation %in% c("Ambiguous", "Mixed"),
                           "Role scores are not well separated; consider secondary roles.",
                           "Role assignment is supported for hypothesis generation."))

  data.frame(
    cell_type = roles$cell_type,
    dominant_role = dominant_role,
    dominant_role_strict = strict_role,
    secondary_role = if ("secondary_role" %in% names(roles)) roles$secondary_role else NA_character_,
    communication_evidence_label = evidence,
    role_reliability_label = reliability,
    role_separation_label = separation,
    role_interpretation = biological,
    role_interpretation_caution = caution,
    stringsAsFactors = FALSE
  )
}

#' Diagnose Cell-Type Communication Results
#'
#' @param ct_comm Cell-type communication object.
#' @param null_pair Optional permutation-null result.
#' @param sens Optional sensitivity-analysis result.
#' @return LogicCommDiagnostic list.
#' @export
diagnose_celltype_communication <- function(ct_comm,
                                            null_pair = NULL,
                                            sens = NULL) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  lr <- ct_comm$lr_table
  pair <- ct_comm$pair_summary
  roles <- ct_comm$role_summary
  out <- list(
    n_cells = length(ct_comm$cell_labels),
    n_cell_types = length(unique(as.character(ct_comm$cell_labels))),
    n_lr_rows = nrow(lr),
    n_celltype_pairs = nrow(pair),
    n_active_events = sum(lr$active %in% TRUE, na.rm = TRUE),
    n_distal_candidates = if ("distal_candidate" %in% names(lr)) sum(lr$distal_candidate %in% TRUE, na.rm = TRUE) else NA_integer_,
    n_low_reliability_roles = if (!is.null(roles) && "role_reliability_label" %in% names(roles)) sum(roles$role_reliability_label %in% c("Low", "Low-communication"), na.rm = TRUE) else NA_integer_,
    role_interpretation = interpret_celltype_roles(ct_comm),
    null_diagnostic = if (!is.null(null_pair)) diagnose_permutation_resolution(null_pair) else NULL,
    sensitivity = sens
  )
  class(out) <- "LogicCommDiagnostic"
  out
}

#' @export
print.LogicCommDiagnostic <- function(x, ...) {
  cat(sprintf("LogicCommDiagnostic | %d cells | %d active LR events\n",
              x$n_cells, x$n_active_events))
  if (!is.null(x$n_distal_candidates) && is.finite(x$n_distal_candidates)) {
    cat(sprintf("Distal/global candidates: %d\n", x$n_distal_candidates))
  }
  if (!is.null(x$n_low_reliability_roles) && is.finite(x$n_low_reliability_roles)) {
    cat(sprintf("Low-reliability roles: %d\n", x$n_low_reliability_roles))
  }
  invisible(x)
}
