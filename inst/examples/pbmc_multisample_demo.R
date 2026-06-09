# PBMC3K four-sample LogicComm demo
#
# This script derives four sample-level PBMC matrices from SeuratData::pbmc3k,
# applies a controlled Case perturbation to pbmc1 and pbmc2, and validates that
# LogicComm recovers the simulated ligand-receptor signal at the sample level.
# It is intended as an executable companion to inst/tutorials/PBMC_multisample_demo.Rmd.

load_pbmc3k_counts <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The PBMC demo requires the Seurat package.", call. = FALSE)
  }
  if (!requireNamespace("pbmc3k.SeuratData", quietly = TRUE)) {
    stop(
      "The PBMC demo requires the pbmc3k.SeuratData package. Install it with SeuratData::InstallData('pbmc3k').",
      call. = FALSE
    )
  }

  data("pbmc3k", package = "pbmc3k.SeuratData", envir = environment())
  pbmc <- Seurat::UpdateSeuratObject(pbmc3k)
  counts <- tryCatch(
    Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts"),
    error = function(e) Seurat::GetAssayData(pbmc, assay = "RNA", slot = "counts")
  )
  counts
}

split_pbmc3k_samples <- function(counts, n_samples = 4, seed = 20260526) {
  stopifnot(n_samples >= 2, !is.null(colnames(counts)))
  set.seed(seed)
  shuffled_cells <- sample(colnames(counts))
  sample_ids <- rep(paste0("pbmc", seq_len(n_samples)), length.out = length(shuffled_cells))
  cell_groups <- split(shuffled_cells, sample_ids)
  out <- lapply(cell_groups, function(cells) counts[, cells, drop = FALSE])
  out[paste0("pbmc", seq_len(n_samples))]
}

boost_pbmc_lr_axis <- function(mat,
                               genes = c("MIF", "CD74", "LGALS3", "CD44", "HLA-A", "CD8A"),
                               lambda = 50,
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  genes <- intersect(genes, rownames(mat))
  for (gene in genes) {
    mat[gene, ] <- mat[gene, ] + Matrix::Matrix(stats::rpois(ncol(mat), lambda), nrow = 1, sparse = TRUE)
  }
  mat
}

run_pbmc_multisample_demo <- function(seed = 20260526,
                                      perturb_lambda = 50,
                                      lcs_threshold = 0.2,
                                      rank_threshold_grid = c(0.4, 0.5, 0.6, 0.7, 0.8),
                                      run_rank_evidence = TRUE,
                                      verbose = TRUE) {
  if (!requireNamespace("LogicComm", quietly = TRUE)) {
    stop("The PBMC demo requires LogicComm to be installed.", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The PBMC demo requires Matrix.", call. = FALSE)
  }

  counts <- load_pbmc3k_counts()
  pbmc_list <- split_pbmc3k_samples(counts, n_samples = 4, seed = seed)

  pbmc_list$pbmc1 <- boost_pbmc_lr_axis(pbmc_list$pbmc1, lambda = perturb_lambda, seed = seed + 1L)
  pbmc_list$pbmc2 <- boost_pbmc_lr_axis(pbmc_list$pbmc2, lambda = perturb_lambda, seed = seed + 2L)

  group_info <- c(pbmc1 = "Case", pbmc2 = "Case", pbmc3 = "Ctrl", pbmc4 = "Ctrl")

  result <- LogicComm::logic_run(
    pbmc_list,
    group_info = group_info,
    lr_db = LogicComm::lr_pairs_human,
    rank_threshold = 0.5,
    lcs_threshold = lcs_threshold,
    case_label = "Case",
    ctrl_label = "Ctrl",
    chunk_size = 2000,
    mc_cores = 1,
    min_samples_per_group = 2,
    verbose = verbose
  )

  targets <- c("MIF_CD74", "MIF_CD44", "LGALS3_CD44", "HLA-A_CD8A", "HLA-B_CD8A")
  target_rows <- result$comparison[result$comparison$lr_pair %in% targets, , drop = FALSE]
  target_rows <- target_rows[order(-target_rows$asymmetry, -target_rows$log2fc_lcs), , drop = FALSE]

  rank_results <- rank_comparison <- rank_target_rows <- NULL
  if (isTRUE(run_rank_evidence)) {
    lr_genes <- LogicComm::logic_get_lr_genes(LogicComm::lr_pairs_human)
    rank_results <- lapply(names(pbmc_list), function(sname) {
      if (verbose) message(sprintf("[PBMC rank demo] Scoring rank-aware REO evidence for %s...", sname))
      reo <- LogicComm::logic_prepare(
        pbmc_list[[sname]],
        genes = lr_genes,
        lr_db = NULL,
        rank_threshold = 0.5,
        return_rank = TRUE,
        rank_output = "percentile",
        chunk_size = 2000,
        verbose = FALSE
      )
      LogicComm::IdentifyRankLogicConsensus(
        reo,
        lr_db = LogicComm::lr_pairs_human,
        threshold_grid = rank_threshold_grid,
        lcs_threshold = 0.01,
        verbose = FALSE
      )
    })
    names(rank_results) <- names(pbmc_list)
    rank_comparison <- LogicComm::CompareRankLogicGroups(
      rank_results,
      group_info = group_info,
      lr_db = LogicComm::lr_pairs_human,
      verbose = FALSE
    )
    rank_target_rows <- rank_comparison[rank_comparison$lr_pair %in% targets, , drop = FALSE]
    rank_target_rows <- rank_target_rows[order(-rank_target_rows$delta_rank_score), , drop = FALSE]
  }

  list(
    pbmc_list = pbmc_list,
    group_info = group_info,
    result = result,
    target_rows = target_rows,
    rank_results = rank_results,
    rank_comparison = rank_comparison,
    rank_target_rows = rank_target_rows,
    parameters = list(
      seed = seed,
      perturb_lambda = perturb_lambda,
      lcs_threshold = lcs_threshold,
      rank_threshold_grid = rank_threshold_grid,
      run_rank_evidence = run_rank_evidence
    )
  )
}

if (sys.nframe() == 0L) {
  demo <- run_pbmc_multisample_demo(verbose = TRUE)
  display_cols <- intersect(
    c("lr_pair", "pathway", "case_freq", "ctrl_freq", "asymmetry", "case_mean_lcs",
      "ctrl_mean_lcs", "log2fc_lcs", "fdr_fisher", "n_case_avail", "n_ctrl_avail"),
    names(demo$target_rows)
  )
  print(demo$target_rows[, display_cols, drop = FALSE], row.names = FALSE)
  if (!is.null(demo$rank_target_rows)) {
    rank_cols <- intersect(
      c("lr_pair", "pathway", "case_mean_rank_score", "ctrl_mean_rank_score",
        "delta_rank_score", "log2fc_rank_score", "case_mean_binary_lcs", "ctrl_mean_binary_lcs"),
      names(demo$rank_target_rows)
    )
    cat("\nRank-aware 0.7 evidence for targeted LR pairs:\n")
    print(demo$rank_target_rows[, rank_cols, drop = FALSE], row.names = FALSE)
  }
  cat("\nDemo validation: ", sum(demo$target_rows$asymmetry == 1, na.rm = TRUE),
      " targeted LR pairs are active in both Case PBMC samples and inactive in both Ctrl PBMC samples at the selected threshold.\n",
      sep = "")
  if (!is.null(demo$rank_target_rows)) {
    cat("Rank validation: ", sum(demo$rank_target_rows$delta_rank_score > 0, na.rm = TRUE),
        " targeted LR pairs have stronger rank-aware evidence in Case than Ctrl.\n",
        sep = "")
  }
}
