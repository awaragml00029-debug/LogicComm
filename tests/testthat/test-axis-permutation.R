library(testthat)
library(LogicComm)

test_that("permute_celltype_communication builds an axis-level null with per-axis p and FDR", {
  set.seed(1)
  genes <- c("L1", "R1", "L2", "R2", paste0("BG", seq_len(10)))
  nc <- 60
  m <- matrix(stats::rpois(length(genes) * nc, 1), nrow = length(genes), ncol = nc,
              dimnames = list(genes, paste0("C", seq_len(nc))))
  # Cell-type-specific axis L1 (A) -> R1 (B); broad axis L2/R2 active everywhere.
  m["L1", 1:30] <- 8
  m["R1", 31:60] <- 8
  m["L2", ] <- 8
  m["R2", ] <- 8
  reo <- calc_REO_matrix(m, lr_genes = genes, verbose = FALSE)
  labs <- stats::setNames(rep(c("A", "B"), each = 30), paste0("C", seq_len(nc)))
  lr <- data.frame(lr_pair = c("L1_R1", "L2_R2"), ligand = c("L1", "L2"),
                   receptor = c("R1", "R2"), pathway = c("P", "P"),
                   stringsAsFactors = FALSE)
  lr$ligand_genes <- list("L1", "L2")
  lr$receptor_genes <- list("R1", "R2")

  ct <- summarize_celltype_communication(reo, cell_labels = labs, lr_db = lr,
                                         min_edges = 1, verbose = FALSE)
  np <- permute_celltype_communication(ct, reo_mat = reo, n_perm = 99, seed = 1,
                                       verbose = FALSE)

  # Axis-level null: carries lr_pair, an empirical p, and a BH FDR per axis.
  expect_true(all(c("sender_type", "receiver_type", "lr_pair",
                    "empirical_p", "fdr") %in% names(np)))
  # The same cell-type pair carries multiple L-R axis rows (not a single shared p).
  expect_true(any(duplicated(paste(np$sender_type, np$receiver_type))))

  # rank_communication_axes joins the per-axis p and FDR onto each L-R axis.
  r <- rank_communication_axes(ct, null_pair = np)
  expect_true(all(c("permutation_empirical_p", "permutation_fdr") %in% names(r)))
})
