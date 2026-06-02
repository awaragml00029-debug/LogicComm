library(testthat)
library(LogicComm)

# Tests for the integrated discovery ranking (phase 3/3 of 0.8.0).

make_discovery_ct <- function() {
  reo <- Matrix::Matrix(matrix(c(
    1, 1, 0, 0, 0, 0,   # L1 active in A cells
    0, 0, 1, 1, 0, 0,   # R1 active in B cells
    1, 0, 0, 0, 1, 0,   # L2 spread across A and C
    0, 0, 1, 0, 0, 1),  # R2 spread across B and C
    nrow = 4, byrow = TRUE,
    dimnames = list(c("L1", "R1", "L2", "R2"), paste0("C", 1:6))), sparse = TRUE)
  labels <- stats::setNames(rep(c("A", "B", "Cc"), each = 2), paste0("C", 1:6))
  lr_db <- data.frame(lr_pair = c("L1_R1", "L2_R2"), ligand = c("L1", "L2"),
                      receptor = c("R1", "R2"), pathway = c("PWa", "PWb"),
                      stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L1", "L2")
  lr_db$receptor_genes <- list("R1", "R2")
  summarize_celltype_communication(reo, cell_labels = labels, lr_db = lr_db, mode = "global",
                                   lcs_threshold = 0.1, min_edges = 1, verbose = FALSE)
}

test_that("rank_communication_axes returns an ordered discovery ranking", {
  ct <- make_discovery_ct()
  ranked <- rank_communication_axes(ct)
  expect_true(all(c("discovery_score", "evidence_tier", "strength_score",
                    "specificity_score", "publication_priority_score") %in% names(ranked)))
  expect_false(is.unsorted(rev(ranked$discovery_score)))  # non-increasing
  # Backward-compatible alias.
  expect_equal(ranked$publication_priority_score, ranked$discovery_score)
})

test_that("rank_communication_axes integrates permutation and stability evidence", {
  ct <- make_discovery_ct()
  base <- rank_communication_axes(ct)
  feats <- paste(base$sender_type, base$receiver_type, base$lr_pair, sep = "|")
  null_pair <- data.frame(sender_type = c("A", "A"), receiver_type = c("B", "A"),
                          empirical_p = c(0.01, 0.9), stringsAsFactors = FALSE)
  sens <- list(stability = data.frame(
    feature = feats, active_fraction = seq(1, 0.2, length.out = length(feats)),
    stringsAsFactors = FALSE))
  ranked <- rank_communication_axes(ct, null_pair = null_pair, sens = sens)

  expect_true(any(is.finite(ranked$null_support)))
  expect_true(any(is.finite(ranked$threshold_stability)))
  # The A -> B : L1_R1 axis is strong, specific, stable, and supported -> Tier 1 and top.
  top <- ranked[1, ]
  expect_equal(top$sender_type, "A")
  expect_equal(top$receiver_type, "B")
  expect_equal(top$lr_pair, "L1_R1")
  expect_match(top$evidence_tier, "Tier 1")
  # Permutation support joins by sender->receiver subgroup (1 - empirical_p).
  ab <- ranked[ranked$sender_type == "A" & ranked$receiver_type == "B", ]
  expect_true(all(abs(ab$null_support - 0.99) < 1e-8))
})

test_that("plot_communication_discovery returns a ggplot", {
  ct <- make_discovery_ct()
  ranked <- rank_communication_axes(ct)
  expect_s3_class(plot_communication_discovery(ranked), "ggplot")
})
