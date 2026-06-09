library(testthat)
library(LogicComm)

test_that("discover_celltype_communication runs end to end with an axis-level null", {
  set.seed(1)
  types <- paste0("T", 1:5)
  per <- 40
  nc <- length(types) * per
  labs <- stats::setNames(rep(types, each = per), paste0("C", seq_len(nc)))
  genes <- c("L1", "R1", "UL", "UR", paste0("BG", seq_len(16)))
  m <- matrix(stats::rpois(length(genes) * nc, 2), nrow = length(genes), ncol = nc,
              dimnames = list(genes, names(labs)))
  # Specific axis: L1 only in T1, R1 only in T2 (clearly above background).
  m["L1", ] <- 0; m["L1", labs == "T1"] <- stats::rpois(per, 12)
  m["R1", ] <- 0; m["R1", labs == "T2"] <- stats::rpois(per, 12)
  # Ubiquitous confound everywhere.
  m["UL", ] <- stats::rpois(nc, 12); m["UR", ] <- stats::rpois(nc, 12)
  reo <- calc_REO_matrix(m, lr_genes = genes, verbose = FALSE)
  lr <- data.frame(lr_pair = c("L1_R1", "UL_UR"), ligand = c("L1", "UL"),
                   receptor = c("R1", "UR"), pathway = c("specific", "ubi"),
                   stringsAsFactors = FALSE)
  lr$ligand_genes <- list("L1", "UL")
  lr$receptor_genes <- list("R1", "UR")

  res <- discover_celltype_communication(reo, cell_labels = labs, lr_db = lr,
                                         n_perm = 99, n_cores = 1L, seed = 1,
                                         min_edges = 1, verbose = FALSE)

  expect_true(all(c("ct_comm", "null", "ranked", "view", "shortlist") %in% names(res)))
  # Axis-level null by default (carries lr_pair) with per-axis p and FDR.
  expect_true("lr_pair" %in% names(res$null))
  expect_true(all(c("permutation_empirical_p", "permutation_fdr") %in% names(res$ranked)))

  # The specificity layer separates the specific axis from the ubiquitous one:
  # the ubiquitous pair is flagged broad, the specific T1 -> T2 axis is not, and
  # it survives the confound-filtered discovery view.
  ub <- res$ranked[res$ranked$lr_pair == "UL_UR", ]
  sp <- res$ranked[res$ranked$lr_pair == "L1_R1" &
                     res$ranked$sender_type == "T1" & res$ranked$receiver_type == "T2", ]
  expect_true(all(ub$broad_axis_flag))
  expect_false(isTRUE(sp$broad_axis_flag))
  expect_true("L1_R1" %in% res$view$lr_pair)
})
