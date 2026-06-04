library(testthat)
library(LogicComm)

# Regression tests for the cell-type expressing-fraction gate (v0.9.2).
#
# Background: neighborhood LCS is an edge fraction. A handful of high-degree
# (hub) receiver cells that carry residual ambient receptor signal can make a
# whole sender -> receiver L-R axis look active even though the receptor is
# expressed in a negligible fraction of that cell type (e.g. ambient CD8A in a
# couple of Treg cells that happen to be KNN hubs). The gene-background guard
# removes most ambient calls, but the few that survive can still dominate the
# edge fraction. The min_expr_frac gate requires the ligand/receptor to be
# active in a minimum fraction of sender/receiver cells, which is
# graph-independent and therefore immune to the hub effect.

make_hub_fixture <- function(seed = 7, hub_count = 30L) {
  set.seed(seed)
  nCD8 <- 60; nTreg <- 60
  cells <- c(paste0("CD8T", seq_len(nCD8)), paste0("Treg", seq_len(nTreg)))
  genes <- c("CD8A", "HLAE", paste0("BG", 1:100))
  m <- matrix(rpois(length(genes) * length(cells), 0.4), length(genes), length(cells),
              dimnames = list(genes, cells))
  for (hk in paste0("BG", 1:6)) m[hk, ] <- rpois(ncol(m), 25)
  CD8 <- seq_len(nCD8); TREG <- (nCD8 + 1):(nCD8 + nTreg)
  m["CD8A", ] <- 0L
  m["CD8A", CD8] <- rpois(nCD8, 10)            # genuine CD8A in CD8 T cells
  hub_treg <- TREG[1:2]                        # 2/60 = 3.3% of Treg, contaminated
  m["CD8A", hub_treg] <- 30L                   # high counts -> survive the bg gate
  m["HLAE", ] <- rpois(ncol(m), 5)             # broadly expressed ligand
  labs <- stats::setNames(c(rep("CD8T", nCD8), rep("Treg", nTreg)), cells)

  # KNN graph where the 2 contaminated Treg cells are hubs for many cells.
  knn <- matrix(0, ncol(m), ncol(m), dimnames = list(cells, cells))
  for (j in seq_len(ncol(m))) {
    nbrs <- unique(c(hub_treg, sample(seq_len(ncol(m)), 8)))
    knn[nbrs, j] <- 1
  }
  lr_db <- data.frame(lr_pair = "HLAE_CD8A", ligand = "HLAE", receptor = "CD8A",
                      pathway = "MHC", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("HLAE"); lr_db$receptor_genes <- list("CD8A")
  list(m = m, labs = labs, knn = knn, lr_db = lr_db, TREG = TREG, CD8 = CD8)
}

axis_row <- function(ct, s, r) {
  t <- ct$lr_table
  t[t$sender_type == s & t$receiver_type == r, , drop = FALSE]
}

test_that("min_expr_frac gates a hub-inflated spurious axis but keeps genuine signal", {
  d <- make_hub_fixture()
  reo <- calc_REO_matrix(d$m, lr_genes = all_lr_genes(d$lr_db),
                         gene_background = "quantile", verbose = FALSE)

  # The background gate leaves a small but nonzero CD8A fraction in Treg
  # (the high-count hub cells survive it): this is exactly the residue the
  # min_expr_frac gate must catch.
  treg_frac <- mean(as.numeric(reo["CD8A", ])[d$TREG])
  expect_gt(treg_frac, 0)
  expect_lt(treg_frac, 0.1)
  expect_gt(mean(as.numeric(reo["CD8A", ])[d$CD8]), 0.3)

  off <- summarize_celltype_communication(reo, cell_labels = d$labs, knn_mat = d$knn,
           lr_db = d$lr_db, mode = "neighborhood", lcs_threshold = 0.01,
           min_edges = 5, min_expr_frac = 0, verbose = FALSE)
  on <- summarize_celltype_communication(reo, cell_labels = d$labs, knn_mat = d$knn,
           lr_db = d$lr_db, mode = "neighborhood", lcs_threshold = 0.01,
           min_edges = 5, min_expr_frac = 0.1, verbose = FALSE)

  spur_off <- axis_row(off, "CD8T", "Treg")
  spur_on  <- axis_row(on,  "CD8T", "Treg")
  # Without the gate the spurious axis is called active; with it, it is not.
  expect_true(spur_off$active)
  expect_false(spur_on$active)
  # The receptor fraction column is reported and is what drives the gate.
  expect_lt(spur_on$receptor_active_frac_receiver, 0.1)

  # The genuine homotypic CD8 T axis (CD8A widely expressed) is preserved.
  gen_on <- axis_row(on, "CD8T", "CD8T")
  expect_true(gen_on$active)
  expect_gte(gen_on$receptor_active_frac_receiver, 0.1)
})

test_that("min_expr_frac is validated and recorded in params", {
  d <- make_hub_fixture()
  reo <- calc_REO_matrix(d$m, lr_genes = all_lr_genes(d$lr_db),
                         gene_background = "quantile", verbose = FALSE)
  expect_error(
    summarize_celltype_communication(reo, cell_labels = d$labs, lr_db = d$lr_db,
      mode = "global", min_expr_frac = 1.5, verbose = FALSE)
  )
  expect_error(
    summarize_celltype_communication(reo, cell_labels = d$labs, lr_db = d$lr_db,
      mode = "global", min_expr_frac = -0.1, verbose = FALSE)
  )
  ct <- summarize_celltype_communication(reo, cell_labels = d$labs, lr_db = d$lr_db,
    mode = "global", min_expr_frac = 0.2, verbose = FALSE)
  expect_equal(ct$params$min_expr_frac, 0.2)
})

test_that("permute_celltype_communication inherits min_expr_frac from ct_comm", {
  d <- make_hub_fixture()
  reo <- calc_REO_matrix(d$m, lr_genes = all_lr_genes(d$lr_db),
                         gene_background = "quantile", verbose = FALSE)
  ct <- summarize_celltype_communication(reo, cell_labels = d$labs, knn_mat = d$knn,
          lr_db = d$lr_db, mode = "neighborhood", lcs_threshold = 0.01,
          min_edges = 5, min_expr_frac = 0.1, verbose = FALSE)
  # No warning about mismatched scoring parameters: the null inherits the gate.
  expect_no_warning(
    null_pair <- permute_celltype_communication(
      ct_comm = ct, reo_mat = reo, knn_mat = d$knn,
      n_perm = 3, metric = "sum_lcs", seed = 1, verbose = FALSE
    )
  )
})
