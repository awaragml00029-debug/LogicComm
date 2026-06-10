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

test_that("PageRank influencer score is non-degenerate on an acyclic graph", {
  # On a feed-forward graph X -> Y -> Z, eigen_centrality(directed = TRUE)
  # returned all zeros. PageRank must give finite, differentiated scores.
  adj <- matrix(0, 3, 3, dimnames = list(c("X", "Y", "Z"), c("X", "Y", "Z")))
  adj["X", "Y"] <- 3
  adj["Y", "Z"] <- 3
  rs <- .communication_role_summary(adj, (adj > 0) * 1, min_role_event_count = 0)
  inf <- stats::setNames(rs$influencer_role_score, rs$cell_type)
  expect_true(all(is.finite(inf)))
  expect_gt(max(inf), 0)
  expect_gt(length(unique(round(inf, 6))), 1)
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

test_that("GLM count table keeps binomial success within the total universe", {
  # Cell-type co-expression: L is active in 1 of 2 A cells and R in both B cells,
  # so the A->B opportunity universe is 2 x 2 = 4 with 1 x 2 = 2 co-expressing.
  # The binomial success must never exceed the total opportunity universe.
  reo <- Matrix::Matrix(matrix(c(
    0, 1, 0, 0,
    0, 0, 1, 1), nrow = 2, byrow = TRUE,
    dimnames = list(c("L", "R"), c("A1", "A2", "B1", "B2"))), sparse = TRUE)
  labels <- stats::setNames(c("A", "A", "B", "B"), c("A1", "A2", "B1", "B2"))
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      pathway = "p", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")

  ct <- summarize_celltype_communication(
    reo, cell_labels = labels, lr_db = lr_db,
    lcs_threshold = 0.4, min_edges = 1, min_active_edges = 1, verbose = FALSE
  )
  ctab <- .extract_comm_count_table(ct, level = "celltype_lr")
  expect_true(all(ctab$success <= ctab$total))
  ab <- ctab[ctab$feature == "A|B|L_R", , drop = FALSE]
  expect_equal(nrow(ab), 1)
  expect_true(ab$success <= ab$total)
  expect_true(ab$total > 0)
})

test_that("dominant role can be Mediator or Influencer, not only Sender/Receiver", {
  # The four role scores live on incomparable scales (summed-LCS strength is
  # effectively unbounded; mediator betweenness is normalized to [0, 1]; the
  # PageRank influencer score sums to 1 over cell types, ~1/n each). dominant_role
  # must therefore be derived from per-column rescaled scores, otherwise raw
  # strength always wins and Mediator/Influencer are mathematically unreachable.
  # On the canonical bridge A -> B -> C -> D, B and C are the mediators.
  adj <- matrix(0, 4, 4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  adj["A", "B"] <- 5; adj["B", "C"] <- 5; adj["C", "D"] <- 5; adj["A", "D"] <- 0.5
  rs <- .communication_role_summary(adj, (adj > 0) * 1, min_role_event_count = 0)
  dom <- stats::setNames(rs$dominant_role_strict, rs$cell_type)
  expect_equal(unname(dom[c("B", "C")]), c("Mediator", "Mediator"))
  # Raw role-score columns remain the interpretable, un-rescaled values.
  expect_equal(rs$sender_role_score[rs$cell_type == "A"], 5.5)
  expect_equal(rs$receiver_role_score[rs$cell_type == "D"], 5.5)
})

test_that("run_multisample warns on and ignores legacy neighborhood arguments", {
  mk <- function(seed) {
    set.seed(seed)
    matrix(stats::rpois(2 * 30, 5), nrow = 2,
           dimnames = list(c("L", "R"), paste0("c", seq_len(30))))
  }
  samples <- list(s1 = mk(1), s2 = mk(2), s3 = mk(3), s4 = mk(4))
  groups <- c(s1 = "Case", s2 = "Case", s3 = "Ctrl", s4 = "Ctrl")
  lr_db <- data.frame(lr_pair = "L_R", ligand = "L", receptor = "R",
                      pathway = "P", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("L")
  lr_db$receptor_genes <- list("R")
  expect_warning(
    run_multisample(samples, group_info = groups, lr_db = lr_db,
                    min_samples_per_group = 2, verbose = FALSE,
                    graph_symmetrize = "or"),
    "deprecated and ignored"
  )
})

test_that("calc_rank_shift ranks within the full transcriptome, not the L-R panel", {
  # MIF dominates the whole transcriptome in Case but is background-level in Ctrl,
  # surrounded by ~300 non-panel genes. Ranking within the 2-gene L-R panel (the
  # former behaviour) would put MIF at a mid normalized rank in Case (~0.5) and
  # miss the transcriptome-wide shift; full-transcriptome ranking puts it near 0.
  mk <- function(mif_level, seed) {
    set.seed(seed)
    G <- 300; Ce <- 25
    m <- matrix(stats::rpois(G * Ce, 2) + 1, nrow = G,
                dimnames = list(c("MIF", "CD74", paste0("bg", seq_len(G - 2))),
                                paste0("c", seq_len(Ce))))
    m["MIF", ] <- mif_level
    m
  }
  samples <- list(c1 = mk(80, 1), c2 = mk(80, 2), t1 = mk(3, 3), t2 = mk(3, 4))
  groups <- c(c1 = "Case", c2 = "Case", t1 = "Ctrl", t2 = "Ctrl")
  lr_db <- data.frame(lr_pair = "MIF_CD74", ligand = "MIF", receptor = "CD74",
                      pathway = "P", stringsAsFactors = FALSE)
  lr_db$ligand_genes <- list("MIF")
  lr_db$receptor_genes <- list("CD74")

  rs <- calc_rank_shift(samples, group_info = groups, lr_db = lr_db,
                        min_detection_frac = 0, verbose = FALSE)
  mif <- rs[rs$gene == "MIF", , drop = FALSE]
  expect_equal(nrow(mif), 1)
  expect_lt(mif$median_rank_case, 0.05)   # transcriptome-dominant in Case
  expect_gt(mif$shift_score, 0)           # more prominent (higher rank) in Case
})
