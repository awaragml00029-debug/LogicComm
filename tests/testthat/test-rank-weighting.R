library(testthat)
library(LogicComm)

test_that("lcs_weighting = 'rank' re-scores by REO intensity with a coherent null", {
  set.seed(1)
  types <- paste0("T", 1:4)
  per <- 40
  nc <- length(types) * per
  labs <- stats::setNames(rep(types, each = per), paste0("C", seq_len(nc)))
  genes <- c("L1", "R1", paste0("BG", seq_len(18)))
  m <- matrix(stats::rpois(length(genes) * nc, 2), nrow = length(genes), ncol = nc,
              dimnames = list(genes, names(labs)))
  m["L1", ] <- 0; m["L1", labs == "T1"] <- stats::rpois(per, 12)
  m["R1", ] <- 0; m["R1", labs == "T2"] <- stats::rpois(per, 12)
  lr <- data.frame(lr_pair = "L1_R1", ligand = "L1", receptor = "R1",
                   pathway = "p", stringsAsFactors = FALSE)
  lr$ligand_genes <- list("L1"); lr$receptor_genes <- list("R1")

  reo <- calc_REO_matrix(m, lr_genes = genes, return_rank = TRUE, verbose = FALSE)

  ct_b <- summarize_celltype_communication(reo, cell_labels = labs, lr_db = lr,
                                           min_edges = 1, verbose = FALSE)
  ct_r <- summarize_celltype_communication(reo, cell_labels = labs, lr_db = lr,
                                           min_edges = 1, lcs_weighting = "rank",
                                           verbose = FALSE)

  expect_equal(ct_b$params$lcs_weighting, "binary")
  expect_equal(ct_r$params$lcs_weighting, "rank")
  # The active call is unchanged; only the reported lcs differs.
  expect_identical(ct_b$lr_table$active, ct_r$lr_table$active)
  expect_true(all(ct_r$lr_table$lcs >= 0 & ct_r$lr_table$lcs <= 1))
  row <- function(d) d[d$sender_type == "T1" & d$receiver_type == "T2" & d$lr_pair == "L1_R1", ]
  expect_true(row(ct_r$lr_table)$active)
  expect_gt(row(ct_r$lr_table)$lcs, 0)

  # 'rank' requires the rank matrix.
  expect_error(
    summarize_celltype_communication(reo$logic, cell_labels = labs, lr_db = lr,
                                     lcs_weighting = "rank", verbose = FALSE),
    "rank matrix|return_rank"
  )

  # The permutation null inherits the rank weighting from ct_comm (coherent
  # observed-vs-null comparison) and stays axis-level.
  np <- permute_celltype_communication(ct_r, reo_mat = reo, n_perm = 19, seed = 1,
                                       verbose = FALSE)
  expect_true("lr_pair" %in% names(np))
  expect_true(all(c("empirical_p", "fdr") %in% names(np)))
})
