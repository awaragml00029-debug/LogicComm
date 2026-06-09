# inst/benchmark/benchmark_vs_baselines.R
#
# Stage 6 benchmark harness (sandbox-runnable core + real-package adapters).
#
# Quantifies whether LogicComm recovers TRUE cell-type-resolved ligand-receptor
# communication while suppressing the confounds that dominate naive scores, vs
# baseline methods. We simulate data with a known ground truth and several
# realistic confounds, score every sender -> receiver -> L-R axis with each
# method, and compare AUROC / AUPRC.
#
# Confounds simulated (the ones that break mean-expression methods in practice):
#   * ubiquitous / housekeeping pairs (high in every cell)            -> ubiquitous
#   * a broad, moderate-abundance pair active across many types       -> abundance
#   * a CYCLING / transcriptional-breadth hub cell type that co-expresses many
#     L-R genes and so forms spurious axes with everything            -> breadth hub
#   * uneven cell-type sizes incl. a rare type                        -> abundance/size
#   * TRUE axes planted at three signal strengths (strong/med/weak)   -> sensitivity
#
# In-sandbox methods: LogicComm (binary), LogicComm (rank), LogicComm + the
# proliferation filter, and re-implemented mean-product / naive baselines. Real
# CellChat / LIANA adapters are at the bottom (run where those packages exist;
# they return scores aligned to the same universe, then reuse auroc()/auprc()).
#
# Usage:
#   devtools::load_all("."); source("inst/benchmark/benchmark_vs_baselines.R")
#   res <- run_benchmark(seeds = 1:5); print(res$summary)

suppressPackageStartupMessages(library(Matrix))

# ---- 1. Simulation with ground truth + confounds ----------------------------

simulate_communication_dataset <- function(n_background = 200,
                                            n_true = 12,        # 3 strengths x 4
                                            n_ubiquitous = 5,
                                            n_abundance = 4,
                                            n_null = 12,
                                            seed = 1) {
  set.seed(seed)
  # Uneven sizes incl. a rare type; "Cycling" is the breadth hub.
  sizes <- c(T1 = 120, T2 = 100, T3 = 80, T4 = 45, T5 = 20, Cycling = 70)
  types <- names(sizes)
  noncyc <- setdiff(types, "Cycling")
  labels <- rep(types, times = sizes)
  cell_ids <- paste0("C", seq_along(labels))
  names(labels) <- cell_ids
  nc <- length(labels)

  bg   <- paste0("BG", seq_len(n_background))
  tL <- paste0("TL", seq_len(n_true));  tR <- paste0("TR", seq_len(n_true))
  uL <- paste0("UL", seq_len(n_ubiquitous)); uR <- paste0("UR", seq_len(n_ubiquitous))
  aL <- paste0("AL", seq_len(n_abundance));  aR <- paste0("AR", seq_len(n_abundance))
  nL <- paste0("NL", seq_len(n_null));  nR <- paste0("NR", seq_len(n_null))
  all_genes <- c(bg, tL, tR, uL, uR, aL, aR, nL, nR)
  counts <- matrix(stats::rpois(length(all_genes) * nc, lambda = 1),
                   nrow = length(all_genes), ncol = nc,
                   dimnames = list(all_genes, cell_ids))

  hi <- function(idx, lam) { counts[g, idx] <<- counts[g, idx] + stats::rpois(sum(idx), lam) }

  # TRUE axes between non-cycling types at three signal strengths.
  strengths <- rep(c(10, 6, 3), length.out = n_true)   # strong / medium / weak
  truth <- data.frame(lr_pair = character(), sender_type = character(),
                      receiver_type = character(), strength = numeric(),
                      stringsAsFactors = FALSE)
  for (k in seq_len(n_true)) {
    s <- noncyc[((k - 1) %% length(noncyc)) + 1]
    r <- noncyc[(k %% length(noncyc)) + 1]
    counts[tL[k], labels == s] <- counts[tL[k], labels == s] + stats::rpois(sum(labels == s), strengths[k])
    counts[tR[k], labels == r] <- counts[tR[k], labels == r] + stats::rpois(sum(labels == r), strengths[k])
    truth <- rbind(truth, data.frame(lr_pair = paste0(tL[k], "_", tR[k]),
                                     sender_type = s, receiver_type = r,
                                     strength = strengths[k], stringsAsFactors = FALSE))
  }
  # Ubiquitous: high everywhere.
  for (k in seq_len(n_ubiquitous)) {
    counts[uL[k], ] <- counts[uL[k], ] + stats::rpois(nc, 9)
    counts[uR[k], ] <- counts[uR[k], ] + stats::rpois(nc, 9)
  }
  # Abundance/breadth: moderate, in a broad random subset of types (not specific).
  for (k in seq_len(n_abundance)) {
    tset <- sample(noncyc, 3)
    sel <- labels %in% tset
    counts[aL[k], sel] <- counts[aL[k], sel] + stats::rpois(sum(sel), 5)
    counts[aR[k], sel] <- counts[aR[k], sel] + stats::rpois(sum(sel), 5)
  }
  # Cycling breadth hub: elevated across many genes (broad transcriptome) and
  # co-expresses a random half of ALL L-R genes, so it forms spurious axes with
  # everything. This is what the proliferation/breadth diagnostic must catch.
  cyc <- labels == "Cycling"
  broad_genes <- sample(rownames(counts), floor(0.6 * nrow(counts)))
  counts[broad_genes, cyc] <- counts[broad_genes, cyc] + stats::rpois(length(broad_genes) * sum(cyc), 4)
  lr_genes_all <- c(tL, tR, uL, uR, aL, aR, nL, nR)
  cyc_lr <- sample(lr_genes_all, floor(0.5 * length(lr_genes_all)))
  counts[cyc_lr, cyc] <- counts[cyc_lr, cyc] + stats::rpois(length(cyc_lr) * sum(cyc), 8)

  lr_db <- data.frame(
    lr_pair  = c(paste0(tL, "_", tR), paste0(uL, "_", uR),
                 paste0(aL, "_", aR), paste0(nL, "_", nR)),
    ligand   = c(tL, uL, aL, nL),
    receptor = c(tR, uR, aR, nR),
    pathway  = c(rep("true", n_true), rep("ubiquitous", n_ubiquitous),
                 rep("abundance", n_abundance), rep("null", n_null)),
    stringsAsFactors = FALSE
  )
  lr_db$ligand_genes   <- as.list(lr_db$ligand)
  lr_db$receptor_genes <- as.list(lr_db$receptor)

  list(counts = Matrix::Matrix(counts, sparse = TRUE), labels = labels,
       lr_db = lr_db, truth = truth, types = types)
}

make_universe <- function(sim) {
  grid <- expand.grid(sender_type = sim$types, receiver_type = sim$types,
                      lr_pair = sim$lr_db$lr_pair, stringsAsFactors = FALSE)
  k <- function(s, r, p) paste(s, r, p, sep = "|")
  grid$key <- k(grid$sender_type, grid$receiver_type, grid$lr_pair)
  grid$is_true <- grid$key %in% k(sim$truth$sender_type, sim$truth$receiver_type, sim$truth$lr_pair)
  grid
}

.align <- function(named_scores, universe) {
  unname(ifelse(universe$key %in% names(named_scores), named_scores[universe$key], 0))
}

# ---- 2. LogicComm scorers ----------------------------------------------------

.lc_rank_discovery <- function(ct) {
  r <- rank_communication_axes(ct)
  list(r = r, key = paste(r$sender_type, r$receiver_type, r$lr_pair, sep = "|"))
}

score_logiccomm <- function(sim, universe) {
  reo  <- calc_REO_matrix(sim$counts, lr_genes = rownames(sim$counts),
                          return_rank = TRUE, verbose = FALSE)
  ctb  <- summarize_celltype_communication(reo, cell_labels = sim$labels,
                                           lr_db = sim$lr_db, include_self = FALSE, verbose = FALSE)
  ctb  <- score_communication_specificity(ctb, verbose = FALSE)
  ctr  <- summarize_celltype_communication(reo, cell_labels = sim$labels, lr_db = sim$lr_db,
                                           include_self = FALSE, lcs_weighting = "rank", verbose = FALSE)
  ctr  <- score_communication_specificity(ctr, verbose = FALSE)
  ctp  <- diagnose_proliferation_confound(ctb, expr = as.matrix(sim$counts), verbose = FALSE)

  db <- .lc_rank_discovery(ctb); dr <- .lc_rank_discovery(ctr); dp <- .lc_rank_discovery(ctp)
  # proliferation-filtered score: zero out axes whose sender/receiver is a
  # breadth hub (i.e. the filter would drop them).
  pflag <- if ("proliferation_confound_flag" %in% names(dp$r)) dp$r$proliferation_confound_flag %in% TRUE
           else rep(FALSE, nrow(dp$r))
  list(
    `LogicComm (binary)`      = .align(stats::setNames(db$r$discovery_score, db$key), universe),
    `LogicComm (rank)`        = .align(stats::setNames(dr$r$discovery_score, dr$key), universe),
    `LogicComm (+prolif filt)`= .align(stats::setNames(ifelse(pflag, 0, dp$r$discovery_score), dp$key), universe)
  )
}

# ---- 3. Baselines (re-implemented, sandbox-runnable) -------------------------

score_mean_product <- function(sim, universe) {
  expr <- as.matrix(sim$counts)
  tm <- sapply(sim$types, function(t) rowMeans(expr[, sim$labels == t, drop = FALSE]))
  lig <- stats::setNames(sim$lr_db$ligand, sim$lr_db$lr_pair)
  rec <- stats::setNames(sim$lr_db$receptor, sim$lr_db$lr_pair)
  vapply(seq_len(nrow(universe)), function(i)
    tm[lig[universe$lr_pair[i]], universe$sender_type[i]] *
      tm[rec[universe$lr_pair[i]], universe$receiver_type[i]], numeric(1))
}

score_naive_codetect <- function(sim, universe) {
  det <- as.matrix(sim$counts) > 0
  tf <- sapply(sim$types, function(t) rowMeans(det[, sim$labels == t, drop = FALSE]))
  lig <- stats::setNames(sim$lr_db$ligand, sim$lr_db$lr_pair)
  rec <- stats::setNames(sim$lr_db$receptor, sim$lr_db$lr_pair)
  vapply(seq_len(nrow(universe)), function(i)
    tf[lig[universe$lr_pair[i]], universe$sender_type[i]] *
      tf[rec[universe$lr_pair[i]], universe$receiver_type[i]], numeric(1))
}

# ---- 4. Metrics --------------------------------------------------------------

auroc <- function(score, label) {
  pos <- score[label]; neg <- score[!label]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

auprc <- function(score, label) {
  o <- order(score, stats::runif(length(score)), decreasing = TRUE)  # random tie-break
  lab <- label[o]; tp <- cumsum(lab); fp <- cumsum(!lab)
  prec <- tp / (tp + fp); rec <- tp / sum(label)
  sum(diff(c(0, rec)) * c(prec[1], prec)[-1])
}

# sensitivity at a fixed precision-oriented cutoff: recall of true axes within the
# top-K ranked, K = number of true positives (precision-at-k == recall-at-k here).
sens_at_k <- function(score, label) {
  k <- sum(label)
  o <- order(score, stats::runif(length(score)), decreasing = TRUE)
  mean(label[o][seq_len(k)])
}

# ---- 5. Driver ---------------------------------------------------------------

run_benchmark <- function(seeds = 1:5, ...) {
  per_seed <- lapply(seeds, function(s) {
    sim <- simulate_communication_dataset(seed = s, ...)
    u <- make_universe(sim)
    scores <- c(score_logiccomm(sim, u),
                list(`Mean-product (CellPhoneDB-like)` = score_mean_product(sim, u),
                     `Naive co-detection` = score_naive_codetect(sim, u)))
    set.seed(s)
    do.call(rbind, lapply(names(scores), function(m)
      data.frame(method = m, AUROC = auroc(scores[[m]], u$is_true),
                 AUPRC = auprc(scores[[m]], u$is_true),
                 sens_at_k = sens_at_k(scores[[m]], u$is_true), stringsAsFactors = FALSE)))
  })
  all <- do.call(rbind, per_seed)
  summary <- stats::aggregate(cbind(AUROC, AUPRC, sens_at_k) ~ method, all, mean)
  summary[, -1] <- round(summary[, -1], 3)
  summary <- summary[order(-summary$AUPRC, -summary$AUROC), ]
  rownames(summary) <- NULL
  ex <- simulate_communication_dataset(seed = seeds[1], ...)
  list(summary = summary, n_seeds = length(seeds), n_positive = nrow(ex$truth),
       n_universe = length(ex$types)^2 * nrow(ex$lr_db))
}

# ---- 6. Real-package adapters (run where the packages are installed) ---------
# Each returns a named numeric vector keyed by "sender|receiver|lr_pair", aligned
# to make_universe()$key via .align(), then scored with auroc()/auprc(). They do
# not run in the sandbox (packages absent); they are reference implementations.

# CellChat (R). Maps each CellChat interaction back to the simulated lr_pair by
# matching its ligand/receptor symbols.
score_cellchat <- function(sim, universe, ...) {
  if (!requireNamespace("CellChat", quietly = TRUE)) stop("CellChat not installed.")
  data.input <- as.matrix(sim$counts)
  data.input <- log1p(sweep(data.input, 2, pmax(colSums(data.input), 1), "/") * 1e4)  # simple lognorm
  meta <- data.frame(labels = sim$labels, row.names = colnames(data.input))
  cc <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "labels")
  # Build a custom DB from the simulated lr_db so identifiers match the universe.
  inter <- data.frame(
    interaction_name = sim$lr_db$lr_pair,
    ligand = sim$lr_db$ligand, receptor = sim$lr_db$receptor,
    pathway_name = sim$lr_db$pathway, annotation = "sim",
    stringsAsFactors = FALSE)
  rownames(inter) <- inter$interaction_name
  cc@DB <- list(interaction = inter,
                complex = data.frame(), cofactor = data.frame(),
                geneInfo = data.frame(Symbol = rownames(data.input), stringsAsFactors = FALSE))
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc)
  df <- CellChat::subsetCommunication(cc)   # source, target, interaction_name, prob
  key <- paste(df$source, df$target, df$interaction_name, sep = "|")
  .align(stats::setNames(df$prob, key), universe)
}

# LIANA (R). liana_wrap aggregates several CCC methods (incl. a CellPhoneDB-style
# and NATMI/SCA scores); we use the aggregate rank (lower = stronger -> negate).
score_liana <- function(sim, universe, method_col = "aggregate_rank", ...) {
  if (!requireNamespace("liana", quietly = TRUE) ||
      !requireNamespace("SeuratObject", quietly = TRUE)) stop("liana/SeuratObject not installed.")
  so <- SeuratObject::CreateSeuratObject(counts = sim$counts)
  so <- SeuratObject::NormalizeData(so, verbose = FALSE)
  SeuratObject::Idents(so) <- factor(sim$labels[colnames(so)])
  res <- liana::liana_wrap(so, resource = "custom",
                           external_resource = data.frame(
                             source_genesymbol = sim$lr_db$ligand,
                             target_genesymbol = sim$lr_db$receptor))
  agg <- liana::liana_aggregate(res)
  lrp <- paste0(agg$ligand.complex, "_", agg$receptor.complex)
  key <- paste(agg$source, agg$target, lrp, sep = "|")
  score <- -as.numeric(agg[[method_col]])   # rank: smaller is better
  .align(stats::setNames(score, key), universe)
}
# CellPhoneDB: available as a method inside LIANA (resource/method = "cellphonedb")
# for an apples-to-apples R comparison; the standalone Python CLI can also be run
# and its means/pvalues mapped to universe$key the same way.

if (identical(environment(), globalenv()) && !interactive()) {
  res <- run_benchmark(seeds = 1:5)
  cat(sprintf("Universe: %d axes, %d true positives | mean over %d sims\n",
              res$n_universe, res$n_positive, res$n_seeds))
  print(res$summary)
}
