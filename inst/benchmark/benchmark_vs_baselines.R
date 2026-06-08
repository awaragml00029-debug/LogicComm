# inst/benchmark/benchmark_vs_baselines.R
#
# Stage 6 benchmark harness (sandbox-runnable).
#
# Goal: quantify whether LogicComm recovers TRUE cell-type-resolved ligand-receptor
# communication while down-ranking ubiquitous/housekeeping confounds, versus
# baseline cell-cell-communication scores. We simulate data with a known
# ground truth (specific axes, ubiquitous confounds, and null pairs), score every
# sender -> receiver -> L-R axis with each method, and compare AUROC / AUPRC.
#
# Baselines reimplemented here (so the comparison runs without installing other
# packages): a CellPhoneDB/CellChat-style mean-expression product, and a naive
# raw-count co-detection product. Adapters for the real CellChat / CellPhoneDB /
# LIANA packages are stubbed at the bottom -- run those where the packages are
# installed; the ground-truth scoring helpers here are reusable for them.
#
# Usage:
#   devtools::load_all("."); source("inst/benchmark/benchmark_vs_baselines.R")
#   res <- run_benchmark(seed = 1); print(res$summary)

suppressPackageStartupMessages({
  library(Matrix)
})

# ---- 1. Simulation with ground truth ----------------------------------------

simulate_communication_dataset <- function(n_types = 6,
                                            cells_per_type = 80,
                                            n_background = 150,
                                            n_true = 10,
                                            n_ubiquitous = 5,
                                            n_null = 10,
                                            signal = 8,
                                            seed = 1) {
  set.seed(seed)
  types <- paste0("T", seq_len(n_types))
  n_cells <- n_types * cells_per_type
  labels <- rep(types, each = cells_per_type)
  cell_ids <- paste0("C", seq_len(n_cells))
  names(labels) <- cell_ids

  # gene pools
  bg_genes <- paste0("BG", seq_len(n_background))
  true_lig <- paste0("TL", seq_len(n_true));  true_rec <- paste0("TR", seq_len(n_true))
  ubi_lig  <- paste0("UL", seq_len(n_ubiquitous)); ubi_rec <- paste0("UR", seq_len(n_ubiquitous))
  null_lig <- paste0("NL", seq_len(n_null));  null_rec <- paste0("NR", seq_len(n_null))
  all_genes <- c(bg_genes, true_lig, true_rec, ubi_lig, ubi_rec, null_lig, null_rec)

  counts <- matrix(stats::rpois(length(all_genes) * n_cells, lambda = 1),
                   nrow = length(all_genes), ncol = n_cells,
                   dimnames = list(all_genes, cell_ids))

  # TRUE axes: ligand high only in a designated sender type, receptor high only
  # in a designated receiver type (sender != receiver), so the signal is
  # genuinely cell-type-pair specific.
  truth <- data.frame(lr_pair = character(), sender_type = character(),
                       receiver_type = character(), stringsAsFactors = FALSE)
  for (k in seq_len(n_true)) {
    s <- types[((k - 1) %% n_types) + 1]
    r <- types[(k %% n_types) + 1]            # different type
    counts[true_lig[k], labels == s] <- counts[true_lig[k], labels == s] + signal
    counts[true_rec[k], labels == r] <- counts[true_rec[k], labels == r] + signal
    truth <- rbind(truth, data.frame(lr_pair = paste0(true_lig[k], "_", true_rec[k]),
                                     sender_type = s, receiver_type = r,
                                     stringsAsFactors = FALSE))
  }
  # UBIQUITOUS confounds: ligand and receptor high in EVERY cell (housekeeping).
  for (k in seq_len(n_ubiquitous)) {
    counts[ubi_lig[k], ] <- counts[ubi_lig[k], ] + signal
    counts[ubi_rec[k], ] <- counts[ubi_rec[k], ] + signal
  }
  # NULL pairs: left at background.

  lr_db <- data.frame(
    lr_pair  = c(paste0(true_lig, "_", true_rec),
                 paste0(ubi_lig,  "_", ubi_rec),
                 paste0(null_lig, "_", null_rec)),
    ligand   = c(true_lig, ubi_lig, null_lig),
    receptor = c(true_rec, ubi_rec, null_rec),
    pathway  = c(rep("true", n_true), rep("ubiquitous", n_ubiquitous), rep("null", n_null)),
    stringsAsFactors = FALSE
  )
  lr_db$ligand_genes   <- as.list(lr_db$ligand)
  lr_db$receptor_genes <- as.list(lr_db$receptor)

  list(counts = Matrix::Matrix(counts, sparse = TRUE), labels = labels,
       lr_db = lr_db, truth = truth, types = types)
}

# ---- 2. Ground-truth-scored universe -----------------------------------------

# Build the universe of all (sender_type, receiver_type, lr_pair) and mark the
# planted TRUE axes as positives.
make_universe <- function(sim) {
  grid <- expand.grid(sender_type = sim$types, receiver_type = sim$types,
                      lr_pair = sim$lr_db$lr_pair, stringsAsFactors = FALSE)
  key <- function(s, r, p) paste(s, r, p, sep = "|")
  pos <- key(sim$truth$sender_type, sim$truth$receiver_type, sim$truth$lr_pair)
  grid$key <- key(grid$sender_type, grid$receiver_type, grid$lr_pair)
  grid$is_true <- grid$key %in% pos
  grid
}

# ---- 3. Methods: each returns a score per universe key -----------------------

score_logiccomm <- function(sim, universe) {
  reo <- calc_REO_matrix(sim$counts, lr_genes = rownames(sim$counts), verbose = FALSE)
  ct <- summarize_celltype_communication(reo, cell_labels = sim$labels,
                                         lr_db = sim$lr_db, include_self = FALSE,
                                         verbose = FALSE)
  r <- rank_communication_axes(ct)   # discovery_score integrates specificity + broad demotion
  k <- paste(r$sender_type, r$receiver_type, r$lr_pair, sep = "|")
  # raw co-expression LCS and the specificity-aware discovery score
  lcs <- stats::setNames(r$lcs, k)
  disc <- stats::setNames(r$discovery_score, k)
  list(
    `LogicComm (LCS)`       = unname(ifelse(universe$key %in% names(lcs),  lcs[universe$key],  0)),
    `LogicComm (discovery)` = unname(ifelse(universe$key %in% names(disc), disc[universe$key], 0))
  )
}

# CellPhoneDB / CellChat core: product of mean expression of ligand in sender
# type and receptor in receiver type.
score_mean_product <- function(sim, universe) {
  expr <- as.matrix(sim$counts)
  type_mean <- sapply(sim$types, function(t) rowMeans(expr[, sim$labels == t, drop = FALSE]))
  lig <- sim$lr_db$ligand; rec <- sim$lr_db$receptor
  names(lig) <- names(rec) <- sim$lr_db$lr_pair
  vapply(seq_len(nrow(universe)), function(i) {
    type_mean[lig[universe$lr_pair[i]], universe$sender_type[i]] *
      type_mean[rec[universe$lr_pair[i]], universe$receiver_type[i]]
  }, numeric(1))
}

# Naive raw co-detection: frac(ligand detected in sender) x frac(receptor detected in receiver).
score_naive_codetect <- function(sim, universe) {
  expr <- as.matrix(sim$counts) > 0
  type_frac <- sapply(sim$types, function(t) rowMeans(expr[, sim$labels == t, drop = FALSE]))
  lig <- sim$lr_db$ligand; rec <- sim$lr_db$receptor
  names(lig) <- names(rec) <- sim$lr_db$lr_pair
  vapply(seq_len(nrow(universe)), function(i) {
    type_frac[lig[universe$lr_pair[i]], universe$sender_type[i]] *
      type_frac[rec[universe$lr_pair[i]], universe$receiver_type[i]]
  }, numeric(1))
}

# ---- 4. Metrics: AUROC (Mann-Whitney) and AUPRC (step) -----------------------

auroc <- function(score, label) {
  pos <- score[label]; neg <- score[!label]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

auprc <- function(score, label) {
  # Break ties at random (not by index) so saturating scores -- e.g. detection
  # fractions pinned at 1 for both true and ubiquitous axes -- are not silently
  # ordered in the positives' favour.
  o <- order(score, stats::runif(length(score)), decreasing = TRUE)
  lab <- label[o]
  tp <- cumsum(lab); fp <- cumsum(!lab)
  prec <- tp / (tp + fp); rec <- tp / sum(label)
  rec0 <- c(0, rec); prec0 <- c(prec[1], prec)
  sum(diff(rec0) * prec0[-1])
}

# ---- 5. Driver ---------------------------------------------------------------

run_benchmark <- function(seeds = 1:5, ...) {
  per_seed <- lapply(seeds, function(s) {
    sim <- simulate_communication_dataset(seed = s, ...)
    universe <- make_universe(sim)
    scores <- c(score_logiccomm(sim, universe),
                list(`Mean-product (CellPhoneDB-like)` = score_mean_product(sim, universe),
                     `Naive co-detection` = score_naive_codetect(sim, universe)))
    set.seed(s)  # reproducible random tie-breaking in auprc
    do.call(rbind, lapply(names(scores), function(m) {
      data.frame(method = m,
                 AUROC = auroc(scores[[m]], universe$is_true),
                 AUPRC = auprc(scores[[m]], universe$is_true),
                 stringsAsFactors = FALSE)
    }))
  })
  all <- do.call(rbind, per_seed)
  summary <- stats::aggregate(cbind(AUROC, AUPRC) ~ method, data = all, FUN = mean)
  summary$AUROC <- round(summary$AUROC, 3)
  summary$AUPRC <- round(summary$AUPRC, 3)
  summary <- summary[order(-summary$AUPRC, -summary$AUROC), ]
  rownames(summary) <- NULL
  ex <- simulate_communication_dataset(seed = seeds[1], ...)
  list(summary = summary, n_seeds = length(seeds),
       n_positive = nrow(ex$truth),
       n_universe = length(ex$types)^2 * nrow(ex$lr_db))
}

# ---- 6. Adapters for real packages (run where installed) ---------------------
# These intentionally do not run in the sandbox. Convert `sim` to each tool's
# input, obtain its per-(sender,receiver,lr) score, align to `universe$key`, then
# reuse auroc()/auprc() above for an apples-to-apples comparison.
#
# score_cellchat <- function(sim, universe) { ... }   # CellChat::computeCommunProb
# score_liana    <- function(sim, universe) { ... }   # liana::liana_wrap
# score_cellphonedb <- function(sim, universe) { ... } # statistical_analysis

if (identical(environment(), globalenv()) && !interactive()) {
  res <- run_benchmark(seeds = 1:5)
  cat(sprintf("Universe: %d axes, %d true positives | mean over %d sims\n",
              res$n_universe, res$n_positive, res$n_seeds))
  print(res$summary)
}
