library(testthat)
library(LogicComm)

# Regression tests for correctness fixes from the v0.7.x communication audit.

test_that("mediator betweenness ranks bridge cell types above endpoints", {
  # A -> B -> C -> D with strong edges plus a weak A -> D shortcut.
  # B and C are the bridges and must score higher than the endpoints A and D.
  adj_strength <- matrix(0, 4, 4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  adj_strength["A", "B"] <- 5
  adj_strength["B", "C"] <- 5
  adj_strength["C", "D"] <- 5
  adj_strength["A", "D"] <- 0.5
  adj_count <- (adj_strength > 0) * 1

  rs <- .communication_role_summary(adj_strength, adj_count, min_role_event_count = 0)
  med <- stats::setNames(rs$mediator_role_score, rs$cell_type)

  expect_gt(med[["B"]], med[["A"]])
  expect_gt(med[["C"]], med[["D"]])
  expect_true(all(is.finite(rs$mediator_role_score)))
  expect_true(all(is.finite(rs$influencer_role_score)))
})

test_that("single cell type role summary does not produce NaN or NA roles", {
  adj <- matrix(2, 1, 1, dimnames = list("A", "A"))
  rs <- .communication_role_summary(adj, (adj > 0) * 1, min_role_event_count = 0)
  expect_equal(nrow(rs), 1)
  expect_false(is.nan(rs$mediator_role_score))
  expect_false(is.na(rs$dominant_role_strict))
  expect_true(is.finite(rs$influencer_role_score))
})

test_that("built-in lr_pairs_human has unique lr_pair identifiers", {
  data(lr_pairs_human)
  expect_equal(sum(duplicated(lr_pairs_human$lr_pair)), 0L)
  # The previously duplicated POSTN/ITGAV interactions are now distinguished by
  # their receptor subunit.
  expect_true(all(c("POSTN_ITGAV_ITGB1", "POSTN_ITGAV_ITGB3") %in% lr_pairs_human$lr_pair))
})

test_that("rank_comm_cells aligns unnamed cell_labels to ranked cells", {
  scores <- list(
    comm_score     = stats::setNames(c(0.1, 0.9, 0.2), c("c1", "c2", "c3")),
    sender_score   = stats::setNames(c(0.1, 0.9, 0.2), c("c1", "c2", "c3")),
    receiver_score = stats::setNames(c(0.1, 0.9, 0.2), c("c1", "c2", "c3"))
  )
  # Unnamed labels follow the original (c1, c2, c3) scoring order.
  ranked <- rank_comm_cells(scores, n = 3, cell_labels = c("A", "B", "A"))
  expect_equal(ranked$cell[1], "c2")
  expect_equal(ranked$cluster[1], "B")
  expect_equal(ranked$cluster[ranked$cell == "c1"], "A")

  # Named labels (possibly in a different order) must agree with the unnamed path.
  ranked_named <- rank_comm_cells(
    scores, n = 3,
    cell_labels = stats::setNames(c("A", "A", "B"), c("c3", "c1", "c2"))
  )
  expect_equal(ranked_named$cluster[ranked_named$cell == "c2"], "B")
  expect_equal(ranked_named$cluster[ranked_named$cell == "c1"], "A")
})

test_that("GLM count table keeps success within the n_edges universe for distal candidates", {
  # A -> B is a distal candidate: one local edge that is not active, but the
  # ligand/receptor are active globally. The old code set success to the global
  # cell-count product and clamped it to 100% against the neighborhood n_edges.
  reo <- Matrix::Matrix(matrix(c(
    0, 1, 0, 0,
    0, 0, 1, 1), nrow = 2, byrow = TRUE,
    dimnames = list(c("L", "R"), c("A1", "A2", "B1", "B2"))), sparse = TRUE)
  labels <- stats::setNames(c("A", "A", "B", "B"), c("A1", "A2", "B1", "B2"))
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      pathway = "p", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")
  knn <- Matrix::Matrix(0, 4, 4, sparse = TRUE)
  rownames(knn) <- colnames(knn) <- c("A1", "A2", "B1", "B2")
  knn["A1", "B1"] <- 1

  ct <- summarize_celltype_communication(
    reo, cell_labels = labels, knn_mat = knn, lr_db = lr_db,
    lcs_threshold = 0.4, min_edges = 1, min_active_edges = 1, verbose = FALSE
  )
  ctab <- .extract_comm_count_table(ct, level = "celltype_lr")
  expect_true(all(ctab$success <= ctab$total))
  ab <- ctab[ctab$feature == "A|B|L_R", , drop = FALSE]
  expect_equal(nrow(ab), 1)
  # The single local A->B edge is inactive, so the binomial success is 0, not a
  # clamped 100%.
  expect_equal(ab$success, 0)
  expect_equal(ab$total, 1)
})
