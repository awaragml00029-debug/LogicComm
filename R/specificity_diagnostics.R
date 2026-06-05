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
  
  pair_levels <- unique(paste(ct_comm$pair_summary$sender_type, ct_comm$pair_summary$receiver_type, sep = "|||"))
  
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

#' Rank Communication Axes by Integrated Discovery Evidence
#'
#' Produces a single prioritized list of candidate sender -> receiver -> L-R
#' communication axes by integrating the evidence layers LogicComm computes:
#' communication strength, cell-type-pair specificity, REO threshold stability,
#' cell-label permutation support, and (when available) receiver downstream
#' response. This is the end of the discovery workflow: it turns the separate
#' diagnostics into one evidence-tiered ranking.
#'
#' @details
#' Each component is scaled to \code{[0, 1]} and combined into a
#' \code{discovery_score} using weights over whichever components are available
#' (strength 0.35, specificity 0.25, threshold stability 0.15, permutation
#' support 0.15, receiver response 0.10, renormalized). Stability is joined from
#' \code{sens} by the \code{sender|receiver|lr_pair} feature key; permutation
#' support (\code{1 - empirical_p}) is joined from \code{null_pair} by the
#' sender/receiver subgroup. An \code{evidence_tier} summarizes the axes from
#' Tier 1 (strong, specific, and supported) to Tier 4 (weak or context
#' dependent); when stability or permutation evidence is not supplied, the tier
#' degrades gracefully to strength and specificity. Tiers are an interpretation
#' aid for hypothesis prioritization, not hypothesis-test p-values. Broad or
#' ubiquitous axes flagged by \code{\link{score_communication_specificity}} are
#' down-ranked (see \code{demote_broad}) so that cell-type-pair-specific
#' candidates surface above ubiquitous interactions such as MHC-I -> CD8A.
#'
#' @param ct_comm Cell-type communication object.
#' @param null_pair Optional output from \code{\link{permute_celltype_communication}}.
#' @param sens Optional output from \code{\link{sensitivity_REO_threshold}}.
#' @param top_n Number of rows to return.
#' @param demote_broad If \code{TRUE} (default), broad/ubiquitous axes (those with
#'   \code{ubiquitous_interaction_flag == TRUE}, e.g. MHC-I -> CD8A or CD99 - CD99)
#'   are down-ranked: their \code{discovery_score} is multiplied by
#'   \code{broad_penalty} and their \code{evidence_tier} is capped at Tier 3,
#'   because by definition they are not cell-type-pair specific. This changes
#'   ranking only -- the \code{active} flag and LCS are untouched. Set to
#'   \code{FALSE} to restore the previous behaviour.
#' @param broad_penalty Multiplicative penalty in \code{[0, 1]} applied to the
#'   \code{discovery_score} of broad/ubiquitous axes when
#'   \code{demote_broad = TRUE} (default \code{0.5}).
#' @return The active \code{lr_table} rows augmented with \code{strength_score},
#'   \code{specificity_score}, \code{threshold_stability}, \code{null_support},
#'   \code{discovery_score}, \code{broad_axis_flag}, and \code{evidence_tier},
#'   sorted by \code{discovery_score}. \code{publication_priority_score} is kept
#'   as a backward-compatible alias.
#' @seealso \code{\link{plot_communication_discovery}}
#' @export
rank_communication_axes <- function(ct_comm,
                                    null_pair = NULL,
                                    sens = NULL,
                                    top_n = NULL,
                                    demote_broad = TRUE,
                                    broad_penalty = 0.5) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  stopifnot(is.numeric(broad_penalty), length(broad_penalty) == 1,
            broad_penalty >= 0, broad_penalty <= 1)
  if (is.null(ct_comm$specificity_summary)) {
    ct_comm <- score_communication_specificity(ct_comm, verbose = FALSE)
  }

  df <- ct_comm$lr_table
  df <- df[df$active %in% TRUE & is.finite(df$lcs), , drop = FALSE]
  if (!nrow(df)) return(df)

  strength <- .rescale01(df$lcs)
  specificity <- .rescale01(ifelse(is.na(df$pair_specificity), 0, df$pair_specificity))

  # Threshold stability: fraction of REO thresholds at which the axis stays
  # active, joined from sensitivity_REO_threshold() by the feature key.
  feature_key <- paste(df$sender_type, df$receiver_type, df$lr_pair, sep = "|")
  stability <- rep(NA_real_, nrow(df))
  if (!is.null(sens) && !is.null(sens$stability) && "feature" %in% names(sens$stability)) {
    stability <- sens$stability$active_fraction[match(feature_key, sens$stability$feature)]
  }
  df$threshold_stability <- stability

  # Permutation support: 1 - empirical_p from the cell-label shuffle null,
  # joined by sender/receiver subgroup (the null is pair-level).
  support <- rep(NA_real_, nrow(df))
  if (!is.null(null_pair) && all(c("sender_type", "receiver_type", "empirical_p") %in% names(null_pair))) {
    nk <- paste(null_pair$sender_type, null_pair$receiver_type)
    ep <- null_pair$empirical_p[match(paste(df$sender_type, df$receiver_type), nk)]
    support <- 1 - ep
    df$permutation_empirical_p <- ep
  }
  df$null_support <- support

  response <- rep(NA_real_, nrow(df))
  resp_col <- intersect(c("response_integrated_score", "receiver_response_score"), names(df))
  if (length(resp_col)) {
    response <- .rescale01(ifelse(is.na(df[[resp_col[1]]]), 0, df[[resp_col[1]]]))
  }

  comp <- list(strength = strength, specificity = specificity,
               stability = stability, support = support, response = response)
  w0 <- c(strength = 0.35, specificity = 0.25, stability = 0.15, support = 0.15, response = 0.10)
  avail <- vapply(comp, function(x) any(is.finite(x)), logical(1))
  w <- w0[avail]
  w <- w / sum(w)
  score_mat <- vapply(names(w), function(nm) {
    v <- comp[[nm]]; v[!is.finite(v)] <- 0; v
  }, numeric(nrow(df)))
  if (is.null(dim(score_mat))) score_mat <- matrix(score_mat, nrow = nrow(df))
  df$discovery_score <- as.numeric(score_mat %*% w)
  df$strength_score <- strength
  df$specificity_score <- specificity

  has_stab <- any(is.finite(stability))
  has_stat <- any(is.finite(support))
  stab_ok <- if (has_stab) is.finite(stability) & stability >= 0.5 else rep(TRUE, nrow(df))
  stat_ok <- if (has_stat) is.finite(support) & support >= 0.6 else rep(TRUE, nrow(df))
  spec_hi <- specificity >= 0.6
  str_hi <- strength >= 0.5
  tier_rank <- ifelse(str_hi & spec_hi & stab_ok & stat_ok, 1L,
                ifelse(str_hi & (spec_hi | (has_stat & stat_ok)), 2L,
                ifelse(strength >= 0.25 | spec_hi, 3L, 4L)))

  # Broad/ubiquitous axes (e.g. MHC-I -> CD8A, CD99 - CD99) are by definition not
  # cell-type-pair specific. score_communication_specificity() already flags them
  # (ubiquitous_interaction_flag), but the flag used to be annotation only, so a
  # strong, null-supported broad axis still reached Tier 1/2 -- the tier rule
  # never consulted it. When demote_broad = TRUE we multiply discovery_score by
  # broad_penalty and cap the tier at 3, so a broad axis cannot outrank a
  # genuinely pair-specific candidate. Ranking layer only: active and lcs are
  # untouched, and demote_broad = FALSE restores the previous behaviour.
  broad <- if ("ubiquitous_interaction_flag" %in% names(df)) {
    df$ubiquitous_interaction_flag %in% TRUE
  } else rep(FALSE, nrow(df))
  df$broad_axis_flag <- broad
  demoted <- isTRUE(demote_broad) & broad & tier_rank < 3L
  if (isTRUE(demote_broad)) {
    df$discovery_score <- df$discovery_score * ifelse(broad, broad_penalty, 1)
    tier_rank[demoted] <- 3L
  }

  tier_labels <- c("Tier 1: strong, specific, supported",
                   "Tier 2: strong candidate",
                   "Tier 3: emerging candidate",
                   "Tier 4: weak / context-dependent")
  df$evidence_tier <- tier_labels[tier_rank]
  df$evidence_tier[demoted] <- "Tier 3: broad / non-specific"

  df$publication_priority_score <- df$discovery_score
  df <- df[order(df$discovery_score, decreasing = TRUE), , drop = FALSE]
  rownames(df) <- NULL
  if (!is.null(top_n)) df <- utils::head(df, top_n)
  df
}

#' Plot the Communication Discovery Landscape
#'
#' Scatter of candidate communication axes by scaled communication strength
#' versus cell-type-pair specificity, sized by permutation support (or
#' discovery score) and coloured by evidence tier. Strong, specific, supported
#' axes sit in the upper right; broad or weak axes fall to the lower left.
#'
#' @param ranked Output of \code{\link{rank_communication_axes}}.
#' @param top_n_label Number of top axes (by discovery score) to label.
#' @param title Optional plot title.
#' @return A ggplot2 object.
#' @export
plot_communication_discovery <- function(ranked, top_n_label = 15, title = NULL) {
  stopifnot(is.data.frame(ranked))
  needed <- c("strength_score", "specificity_score", "discovery_score",
              "sender_type", "receiver_type", "lr_pair")
  if (!all(needed %in% names(ranked))) {
    stop("ranked must be the output of rank_communication_axes().")
  }
  if (!nrow(ranked)) stop("No active communication axes to plot.")
  df <- ranked
  df$axis <- .compact_feature_label(paste(df$sender_type, df$receiver_type, df$lr_pair, sep = "|"))
  df$tier <- if ("evidence_tier" %in% names(df)) df$evidence_tier else "Unranked"
  use_support <- "null_support" %in% names(df) && any(is.finite(df$null_support))
  size_raw <- if (use_support) df$null_support else df$discovery_score
  df$disc_size <- ifelse(is.finite(size_raw), size_raw, 0)
  df$to_label <- FALSE
  df$to_label[seq_len(min(top_n_label, nrow(df)))] <- TRUE
  if (is.null(title)) title <- "Communication discovery landscape"
  ggplot2::ggplot(df, ggplot2::aes(x = strength_score, y = specificity_score)) +
    ggplot2::geom_point(ggplot2::aes(size = disc_size, color = tier), alpha = 0.8) +
    ggrepel::geom_text_repel(data = df[df$to_label, , drop = FALSE],
                             ggplot2::aes(label = axis), size = 3, max.overlaps = 12,
                             min.segment.length = 0) +
    ggplot2::scale_size_continuous(range = c(2, 7),
                                   name = if (use_support) "Null support" else "Discovery score") +
    ggplot2::labs(title = title, x = "Communication strength (scaled)",
                  y = "Pair specificity (scaled)", color = "Evidence tier") +
    theme_logiccomm()
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
