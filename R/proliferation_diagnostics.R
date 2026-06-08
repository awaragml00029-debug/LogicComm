# R/proliferation_diagnostics.R

# Tirosh et al. (2016) S and G2M cell-cycle gene sets (the lists shipped with
# Seurat as cc.genes). Used as defaults to score proliferation; any genes absent
# from the expression matrix are dropped, so newer symbol aliases simply shrink
# the effective set rather than erroring.
.logiccomm_s_genes <- function() c(
  "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
  "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
  "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
  "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
  "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
  "BRIP1", "E2F8"
)

.logiccomm_g2m_genes <- function() c(
  "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
  "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
  "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
  "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
  "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
  "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
  "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
)

.col_sums <- function(x) if (methods::is(x, "Matrix")) Matrix::colSums(x) else base::colSums(x)
.col_means <- function(x) if (methods::is(x, "Matrix")) Matrix::colMeans(x) else base::colMeans(x)

# Mean expression of a gene set per cell (a simple, interpretable module score;
# the goal is to identify cycling cell types, not to phase individual cells).
.module_score <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  if (!length(genes)) return(rep(NA_real_, ncol(expr)))
  .col_means(expr[genes, , drop = FALSE])
}

.zscore <- function(x) {
  ok <- is.finite(x)
  if (sum(ok) < 2) return(rep(NA_real_, length(x)))
  mu <- mean(x[ok]); sdv <- stats::sd(x[ok])
  if (!is.finite(sdv) || sdv == 0) return(ifelse(is.finite(x), 0, NA_real_))
  (x - mu) / sdv
}

#' Diagnose the Proliferation / Transcriptional-Breadth Confound
#'
#' Cycling cells have unusually broad transcriptomes (many detected genes), so in
#' an REO/logic framework they call a large number of genes "active" and form
#' active ligand-receptor edges with almost every partner. The result is a
#' proliferation "hub" -- e.g. a \code{Cycling_Lymphocytes} cluster appearing in
#' the large majority of top-ranked axes -- that reflects transcriptional breadth
#' rather than specific biology.
#'
#' This is a \emph{diagnostic only}. It scores each cell's S and G2M signatures
#' and transcriptional breadth, aggregates them per cell type, flags cell types
#' that are proliferation/breadth outliers, and annotates each communication axis
#' with whether its sender or receiver is such a hub. It does \strong{not} change
#' the \code{active} flag, the LCS, or the \code{discovery_score}: use the
#' returned flags to filter or down-weight axes yourself, or pass the object on
#' to \code{\link{rank_communication_axes}} (the flags ride along in
#' \code{lr_table}).
#'
#' @param ct_comm Output of \code{\link{summarize_celltype_communication}}.
#' @param expr Gene x cell expression matrix (raw counts or log-normalized).
#'   A base matrix or a \code{Matrix} sparse matrix is accepted. If both
#'   \code{colnames(expr)} and \code{names(ct_comm$cell_labels)} are present the
#'   columns are matched by cell name; otherwise \code{ncol(expr)} must equal
#'   \code{length(ct_comm$cell_labels)} and be in the same cell order.
#' @param s_genes,g2m_genes Character vectors of S-phase and G2M-phase marker
#'   genes. Default to the Tirosh et al. (2016) sets; genes absent from
#'   \code{expr} are dropped.
#' @param hub_z_cutoff A cell type is flagged as a proliferation/breadth hub when
#'   its cycling-score z-score OR its breadth z-score across cell types is at
#'   least this value (default \code{1.5}). z-scores need several cell types to
#'   be meaningful.
#' @param verbose Print a short summary.
#' @return \code{ct_comm} with a new \code{proliferation_summary} data frame (one
#'   row per cell type: \code{n_cells}, \code{mean_s_score}, \code{mean_g2m_score},
#'   \code{mean_cycling_score}, \code{cycling_z}, \code{mean_breadth},
#'   \code{breadth_z}, \code{proliferation_hub_flag}) and four columns added to
#'   \code{lr_table}: \code{sender_proliferation_hub},
#'   \code{receiver_proliferation_hub}, \code{proliferation_confound_flag} (either
#'   end is a hub), and \code{proliferation_hub_role}
#'   ("sender"/"receiver"/"both"/"none").
#' @seealso \code{\link{rank_communication_axes}},
#'   \code{\link{score_communication_specificity}}
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
#' fun <- get("diagnose_proliferation_confound")
#' toy_args <- toy_args[intersect(names(toy_args), names(formals(fun)))]
#' try(do.call(fun, toy_args), silent = TRUE)
#' @export
diagnose_proliferation_confound <- function(ct_comm,
                                            expr,
                                            s_genes = NULL,
                                            g2m_genes = NULL,
                                            hub_z_cutoff = 1.5,
                                            verbose = TRUE) {
  stopifnot(inherits(ct_comm, "LogicCommCellTypeComm"))
  stopifnot(is.numeric(hub_z_cutoff), length(hub_z_cutoff) == 1, is.finite(hub_z_cutoff))
  labels <- ct_comm$cell_labels
  if (is.null(labels)) {
    stop("ct_comm has no cell_labels; cannot map cells to cell types.")
  }
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (is.null(dim(expr))) stop("'expr' must be a gene x cell matrix.")

  # Align expr columns to the cells in ct_comm.
  if (!is.null(names(labels)) && !is.null(colnames(expr)) &&
      length(intersect(names(labels), colnames(expr)))) {
    miss <- setdiff(names(labels), colnames(expr))
    if (length(miss)) {
      stop(sprintf("expr is missing %d cell(s) present in ct_comm$cell_labels (e.g. %s).",
                   length(miss), paste(utils::head(miss, 3), collapse = ", ")))
    }
    expr <- expr[, names(labels), drop = FALSE]
  } else if (ncol(expr) != length(labels)) {
    stop(sprintf(
      "ncol(expr) (%d) must equal length(ct_comm$cell_labels) (%d) and be in the same cell order.",
      ncol(expr), length(labels)))
  }

  if (is.null(s_genes)) s_genes <- .logiccomm_s_genes()
  if (is.null(g2m_genes)) g2m_genes <- .logiccomm_g2m_genes()
  n_s <- length(intersect(s_genes, rownames(expr)))
  n_g2m <- length(intersect(g2m_genes, rownames(expr)))
  if (n_s + n_g2m == 0) {
    warning("None of the S/G2M genes were found in rownames(expr); cycling scores are NA. ",
            "Check that 'expr' uses the same gene symbols as the marker sets.")
  }

  s_score <- .module_score(expr, s_genes)
  g2m_score <- .module_score(expr, g2m_genes)
  cycling_score <- pmax(s_score, g2m_score, na.rm = TRUE)
  cycling_score[!is.finite(cycling_score)] <- NA_real_
  breadth <- as.numeric(.col_sums(expr > 0))   # detected genes per cell

  labels <- as.character(labels)
  cell_types <- ct_comm$role_summary$cell_type
  if (is.null(cell_types)) cell_types <- sort(unique(labels))

  agg <- function(x) vapply(cell_types, function(ct) {
    v <- x[labels == ct]; v <- v[is.finite(v)]
    if (!length(v)) NA_real_ else mean(v)
  }, numeric(1))

  n_cells <- vapply(cell_types, function(ct) sum(labels == ct), integer(1))
  mean_cycling <- agg(cycling_score)
  mean_breadth <- agg(breadth)
  cycling_z <- .zscore(mean_cycling)
  breadth_z <- .zscore(mean_breadth)
  hub_flag <- (is.finite(cycling_z) & cycling_z >= hub_z_cutoff) |
              (is.finite(breadth_z) & breadth_z >= hub_z_cutoff)

  prolif <- data.frame(
    cell_type = cell_types,
    n_cells = as.integer(n_cells),
    mean_s_score = as.numeric(agg(s_score)),
    mean_g2m_score = as.numeric(agg(g2m_score)),
    mean_cycling_score = as.numeric(mean_cycling),
    cycling_z = as.numeric(cycling_z),
    mean_breadth = as.numeric(mean_breadth),
    breadth_z = as.numeric(breadth_z),
    proliferation_hub_flag = as.logical(hub_flag),
    stringsAsFactors = FALSE
  )
  prolif <- prolif[order(prolif$mean_cycling_score, decreasing = TRUE), , drop = FALSE]
  rownames(prolif) <- NULL
  ct_comm$proliferation_summary <- prolif

  # Annotate the lr_table (flags ride along into rank_communication_axes()).
  hub_types <- prolif$cell_type[prolif$proliferation_hub_flag %in% TRUE]
  lr <- ct_comm$lr_table
  if (!is.null(lr) && nrow(lr)) {
    s_hub <- lr$sender_type %in% hub_types
    r_hub <- lr$receiver_type %in% hub_types
    lr$sender_proliferation_hub <- s_hub
    lr$receiver_proliferation_hub <- r_hub
    lr$proliferation_confound_flag <- s_hub | r_hub
    lr$proliferation_hub_role <- ifelse(s_hub & r_hub, "both",
                                 ifelse(s_hub, "sender",
                                 ifelse(r_hub, "receiver", "none")))
    ct_comm$lr_table <- lr
  }

  if (isTRUE(verbose)) {
    frac <- if (!is.null(lr) && nrow(lr) && "active" %in% names(lr)) {
      a <- lr$active %in% TRUE
      if (any(a)) mean(lr$proliferation_confound_flag[a]) else NA_real_
    } else NA_real_
    message(sprintf(
      "[Proliferation] %d/%d cell types flagged as proliferation/breadth hubs%s.%s",
      length(hub_types), nrow(prolif),
      if (length(hub_types)) sprintf(" (%s)", paste(hub_types, collapse = ", ")) else "",
      if (is.finite(frac)) sprintf(" %.0f%% of active axes involve a hub.", 100 * frac) else ""
    ))
  }
  ct_comm
}
